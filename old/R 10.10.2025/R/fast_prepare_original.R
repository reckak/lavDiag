#' Prepare original data + model-based and empirical item curves (continuous & ordinal)
#'
#' Returns .rid, .gid, .group, original indicators, factor scores (+ SEs for
#' continuous-only models), and both model-based yhat CIs (m_*) from augment()
#' and empirical expected-score CIs (e_*) from mgcv::gam() per item.
#'
#' @param fit   Fitted lavaan object.
#' @param data  Optional lavPredict data/list (will be computed if NULL).
#' @param info  Optional model_info(fit) (will be computed if NULL).
#' @param level Confidence level for CIs (default 0.95).
#' @param fam_cont mgcv family for continuous items (default gaussian()).
#' @param fam_ord  mgcv family for ordinal items (default betar(link="logit")).
#' @param ...   Extra args forwarded to mgcv::gam().
#' @return tibble with base columns, m_* columns from augment(), and e_* columns from mgcv.
#' @export

prepare_original <- function(fit,
                             data     = NULL,
                             info     = NULL,
                             level    = 0.95,
                             fam_cont = stats::gaussian(),
                             fam_ord  = mgcv::betar(link = "logit"),
                             gam_args_cont = list(method = "REML"),
                             gam_args_ord  = list(method = "REML"),
                             # --- new parallel args ---
                             plan     = c("auto","multisession","multicore",
                                          "sequential","cluster","none"),
                             workers  = NULL,
                             cluster  = NULL,
                             progress = FALSE) {

  plan <- match.arg(plan)

  .assert_lavaan_fit(fit, require_converged = TRUE, require_latent = TRUE, forbid_multilevel = TRUE)
  if (is.null(info)) info <- model_info(fit)

  ov_all   <- info$observed_variables
  ov_cont  <- info$ov_continuous
  ov_ord   <- info$ov_ordinal
  lats     <- info$latent_variables
  n_groups <- info$n_groups %||% 1L

  build_gam_call <- function(fml, dat, fam, user_args) {
    # Remove any user-supplied fields we lock
    user_args$formula <- NULL
    user_args$data    <- NULL
    user_args$family  <- NULL
    # Compose final argument list (user args take effect for everything else)
    args <- c(list(formula = fml, data = dat, family = fam), user_args)
    # Call mgcv::gam()
    do.call(mgcv::gam, args)
  }

  # --- 1) model-based via augment() (beze změny) -----------------------------
  out <- augment(
    fit            = fit,
    data           = data,
    info           = info,
    yhat           = TRUE, ci = TRUE, resid = FALSE,
    se_yhat        = FALSE, ystar = FALSE, pr = FALSE,
    se_fs          = TRUE,
    prefix_yhat    = "m_est_yhat_",
    prefix_ci      = c("m_lwr_yhat_", "m_upr_yhat_"),
    col_layout     = "by_type"
  )

  # --- helpers (beze změny) --------------------------------------------------
  rescale_01 <- function(x) { r <- range(x, na.rm = TRUE, finite = TRUE)
  if (!is.finite(r[1]) || !is.finite(r[2]) || r[2] == r[1]) return(rep(0.5, length(x)))
  (x - r[1]) / (r[2] - r[1]) }
  sv01 <- function(y01) { n <- sum(!is.na(y01)); if (n <= 1L) return(y01); ((y01 * (n - 1)) + 0.5) / n }
  linkinv_fun <- function(fam) if (!is.null(fam$linkinv)) fam$linkinv else stop("Family lacks linkinv.")

  PE_raw <- lavaan::parameterEstimates(fit, standardized = FALSE)
  PE_sel <- PE_raw[PE_raw$op == "=~" & PE_raw$rhs %in% ov_all, , drop = FALSE]
  keep_cols <- intersect(c("group","lhs","rhs","est"), names(PE_sel))
  PE <- PE_sel[, keep_cols, drop = FALSE]
  if (!"group" %in% names(PE)) PE$group <- 1L

  tol <- 1e-10
  loads_by_item <- tapply(seq_len(nrow(PE)), PE$rhs, function(idx) {
    lhs <- PE$lhs[idx]; est <- PE$est[idx]
    unique(lhs[is.finite(est) & abs(est) > tol])
  })
  loads_by_item <- loads_by_item[ov_all]

  if (!requireNamespace("mgcv", quietly = TRUE)) stop("Package 'mgcv' is required for empirical curves.")
  zcrit <- stats::qnorm(1 - (1 - level) / 2)
  fam_for_item <- function(j) if (j %in% ov_cont) fam_cont else fam_ord
  inv_for_item <- function(j) linkinv_fun(fam_for_item(j))

  # Use full ordinal scale 1..K if a factor; otherwise numeric range
  rng_by_item <- lapply(ov_all, function(j) {
    if (!j %in% names(out)) return(c(0,1))
    x <- out[[j]]
    if (is.factor(x)) c(1L, nlevels(x)) else range(as.numeric(x), na.rm = TRUE, finite = TRUE)
  })
  names(rng_by_item) <- ov_all

  # --- set up future::plan safely --------------------------------------------
  reset_plan <- .set_future_plan(plan = plan, workers = workers, cluster = cluster)
  on.exit(reset_plan(), add = TRUE)

  # --- 2) empirical GAM curves, parallel over items within group -------------
  for (g in seq_len(n_groups)) {
    rows_g <- if (".gid" %in% names(out)) which(out$.gid == g) else seq_len(nrow(out))
    if (!length(rows_g)) next

    eta_cols <- intersect(lats, names(out))
    if (!length(eta_cols)) stop("Factor scores missing in data; cannot build empirical curves.")

    # Freeze per-group data once to keep futures light
    out_g <- out[rows_g, c(ov_all, eta_cols), drop = FALSE]

    item_results <- furrr::future_map(
      ov_all,
      function(j) {
        # -- pick relevant factors for item j
        rel_lats <- loads_by_item[[j]]
        if (is.null(rel_lats) || !length(rel_lats)) rel_lats <- eta_cols
        rel_lats <- rel_lats[rel_lats %in% eta_cols]
        if (!length(rel_lats) || !j %in% names(out_g)) return(NULL)

        df_g <- out_g[, c(j, rel_lats), drop = FALSE]
        fam_j   <- fam_for_item(j)
        invlink <- inv_for_item(j)
        rj      <- rng_by_item[[j]]

        # Prepare response
        if (j %in% ov_cont) {
          y_resp <- df_g[[j]]
          back_to_original <- function(mu_link) invlink(mu_link)
        } else {
          y_resp <- sv01(rescale_01(df_g[[j]]))
          to_original_scale <- function(p01) { d <- rj[2] - rj[1]; p01 * d + rj[1] }
          back_to_original  <- function(mu_link) to_original_scale(invlink(mu_link))
        }

        sm_terms <- paste0("s(", rel_lats, ")", collapse = " + ")
        fml <- stats::as.formula(paste0("y_resp ~ ", sm_terms))

        ok    <- stats::complete.cases(df_g[, c(j, rel_lats), drop = FALSE])
        df_ok <- df_g[ok, , drop = FALSE]
        if (!nrow(df_ok)) {
          n <- nrow(df_g)
          return(list(
            nm_est = paste0("e_est_yhat_", j),
            nm_lwr = paste0("e_lwr_yhat_", j),
            nm_upr = paste0("e_upr_yhat_", j),
            mu  = rep(NA_real_, n),
            lwr = rep(NA_real_, n),
            upr = rep(NA_real_, n)
          ))
        }

        y_fit <- if (j %in% ov_cont) df_ok[[j]] else sv01(rescale_01(df_ok[[j]]))
        X_fit <- df_ok[, rel_lats, drop = FALSE]

        if (nrow(X_fit) >= length(rel_lats) + 2L) {
          # Choose per-item arg set: continuous vs ordinal
          user_args <- if (j %in% ov_cont) gam_args_cont else gam_args_ord

          fit_gam <- build_gam_call(
            fml,
            dat = cbind(y_resp = y_fit, X_fit),
            fam = fam_j,
            user_args = user_args
          )

          pr <- predict(fit_gam,
                        newdata = df_g[, rel_lats, drop = FALSE],
                        type = "link",
                        se.fit = TRUE)
          mu  <- back_to_original(pr$fit)
          lwr <- back_to_original(pr$fit - zcrit * pr$se.fit)
          upr <- back_to_original(pr$fit + zcrit * pr$se.fit)
        } else {
          n <- nrow(df_g)
          mu <- lwr <- upr <- rep(NA_real_, n)
        }

        list(
          nm_est = paste0("e_est_yhat_", j),
          nm_lwr = paste0("e_lwr_yhat_", j),
          nm_upr = paste0("e_upr_yhat_", j),
          mu  = as.numeric(mu),
          lwr = as.numeric(lwr),
          upr = as.numeric(upr)
        )
      },
      .options  = furrr::furrr_options(
        seed     = TRUE,
        # Explicit packages so workers have what they need
        packages = c("mgcv", "stats")
      ),
      .progress = progress
    )

    # Write results back (sequential, small overhead)
    for (res in item_results) {
      if (is.null(res)) next
      out[rows_g, res$nm_est] <- res$mu
      out[rows_g, res$nm_lwr] <- res$lwr
      out[rows_g, res$nm_upr] <- res$upr
    }
  }

  # Ensure numeric columns are double (keep .rid/.gid integer if present)
  is_num <- vapply(out, is.numeric, logical(1L))
  keep_int <- names(out) %in% c(".rid", ".gid")
  to_double <- is_num & !keep_int
  for (nm in names(out)[to_double]) out[[nm]] <- as.double(out[[nm]])

  tibble::as_tibble(out)
}
