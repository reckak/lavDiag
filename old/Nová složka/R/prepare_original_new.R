#' Prepare original data + model-based and empirical item curves (continuous & ordinal)
#'
#' Returns `.rid`, `.gid`, `.group`, original indicators, factor scores
#' (including SEs for continuous-only models), and two kinds of item curves:
#' (i) model-based expected scores + CIs from `augment()` with prefixes `m_*`,
#' and (ii) empirical expected-score curves + CIs from `mgcv::gam()` with prefixes `e_*`.
#'
#' @param fit        Fitted `lavaan` object; must be converged and contain latent variables.
#' @param data       Optional `lavPredict` data/list; computed if `NULL`.
#' @param info       Optional `model_info(fit)`; computed if `NULL`.
#' @param level      Confidence level for empirical CIs (default `0.95`).
#' @param fam_cont   `mgcv` family for continuous items (default `stats::gaussian()`).
#' @param fam_ord    `mgcv` family for ordinal items (default `mgcv::betar(link = "logit")`).
#' @param gam_args_cont List of additional arguments passed to `mgcv::gam()` **for continuous items**.
#'   The function locks `formula`, `data`, and `family`; all other fields (e.g., `method`,
#'   `select`, `discrete`, `gamma`, `knots`, `sp`, `optimizer`, `control`, `scale`,
#'   `weights`, `offset`, ...) can be set here. Default `list(method = "REML")`.
#' @param gam_args_ord  List of additional arguments passed to `mgcv::gam()` **for ordinal items**.
#'   Same rules as above. Default `list(method = "REML")`.
#' @param plan       Parallel plan used via `future`/`furrr`. One of
#'   `c("auto","multisession","multicore","sequential","cluster","none")`.
#'   `"auto"` chooses a cross-platform safe default (multisession).
#' @param workers    Integer number of workers (defaults to `parallel::detectCores()-1` inside
#'   `.set_future_plan()` when `NULL`).
#' @param cluster    Optional cluster object/host spec if `plan = "cluster"`.
#' @param progress   Logical; show `furrr` progress bar (default `FALSE`).
#'
#' @details
#' Empirical curves are fit per item and per group using `mgcv::gam()` on the link scale,
#' then transformed back using the family's inverse link. For ordinal outcomes the response
#' is rescaled to `(0,1)` with the Smithson–Verkuilen adjustment and modeled with `betar`,
#' after which predictions are mapped back to the original item scale.
#'
#' Parallelization: groups are processed sequentially while items within a group are
#' processed in parallel via `furrr::future_map()`. The helper `.set_future_plan()` is used
#' to safely set and restore `future::plan()`.
#'
#' Factor loadings: relevant latent predictors for each item can vary by group. The function
#' first tries group-specific loadings; if missing, it falls back to overall loadings and then
#' to all available factor-score columns. For speed, if loadings **do not** vary across groups,
#' group-specific lookup is skipped entirely.
#'
#' @return A tibble containing base columns, model-based columns (`m_*`) and empirical columns (`e_*`).
#'
#' @examples
#' # Typical call with stronger penalization/selection for ordinal items:
#' # prepare_original(
#' #   fit,
#' #   gam_args_cont = list(method="REML"),
#' #   gam_args_ord  = list(method="REML", select=TRUE, discrete=TRUE)
#' # )
#'
#' @export
prepare_original_old <- function(fit,
                             data     = NULL,
                             info     = NULL,
                             level    = 0.95,
                             fam_cont = stats::gaussian(),
                             fam_ord  = mgcv::betar(link = "logit"),
                             gam_args_cont = list(method = "REML"),
                             gam_args_ord  = list(method = "REML"),
                             plan     = c("auto","multisession","multicore",
                                          "sequential","cluster","none"),
                             workers  = NULL,
                             cluster  = NULL,
                             progress = FALSE) {

  plan <- match.arg(plan)

  # -- Assertions --------------------------------------------------------------
  .assert_lavaan_fit(fit, require_converged = TRUE, require_latent = TRUE, forbid_multilevel = TRUE)
  if (is.null(info)) info <- model_info(fit)

  ov_all   <- info$observed_variables
  ov_cont  <- info$ov_continuous
  ov_ord   <- info$ov_ordinal
  lats     <- info$latent_variables
  n_groups <- info$n_groups %||% 1L

  # -- Utility: safe builder for mgcv::gam calls -------------------------------
  build_gam_call <- function(fml, dat, fam, user_args) {
    # Lock fields we control
    user_args$formula <- NULL
    user_args$data    <- NULL
    user_args$family  <- NULL
    args <- c(list(formula = fml, data = dat, family = fam), user_args)
    do.call(mgcv::gam, args)
  }

  # -- 1) Model-based curves via augment() ------------------------------------
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

  # -- Helpers -----------------------------------------------------------------
  # Numeric rescale to [0,1]
  rescale_01 <- function(x) {
    r <- range(x, na.rm = TRUE, finite = TRUE)
    if (!is.finite(r[1]) || !is.finite(r[2]) || r[2] == r[1]) return(rep(0.5, length(x)))
    (x - r[1]) / (r[2] - r[1])
  }
  # Smithson–Verkuilen adjustment: map [0,1] -> (0,1)
  sv01 <- function(y01) {
    n <- sum(!is.na(y01))
    if (n <= 1L) return(y01)
    ((y01 * (n - 1)) + 0.5) / n
  }
  # Safe inverse link for a mgcv family
  linkinv_fun <- function(fam) {
    if (!is.null(fam$linkinv)) fam$linkinv else stop("Family lacks linkinv.")
  }

  # -- Robust extraction of loadings (factor->item), incl. group-specific -----
  PE_raw <- lavaan::parameterEstimates(fit, standardized = FALSE)
  PE_sel <- PE_raw[PE_raw$op == "=~" & PE_raw$rhs %in% ov_all, , drop = FALSE]
  keep_cols <- intersect(c("group","lhs","rhs","est"), names(PE_sel))
  PE <- PE_sel[, keep_cols, drop = FALSE]
  if (!"group" %in% names(PE)) PE$group <- 1L

  tol <- 1e-10

  # Overall (group-agnostic) item -> latents
  loads_by_item_overall <- tapply(seq_len(nrow(PE)), PE$rhs, function(idx) {
    lhs <- PE$lhs[idx]; est <- PE$est[idx]
    unique(lhs[is.finite(est) & abs(est) > tol])
  })
  loads_by_item_overall <- loads_by_item_overall[ov_all]

  # Group-specific: list keyed by group id -> (named list item -> latents)
  PE_by_group <- split(PE, PE$group)
  loads_by_item_by_group <- lapply(PE_by_group, function(PEg) {
    if (!nrow(PEg)) return(setNames(vector("list", length(ov_all)), ov_all))
    li <- tapply(seq_len(nrow(PEg)), PEg$rhs, function(idx) {
      lhs <- PEg$lhs[idx]; est <- PEg$est[idx]
      unique(lhs[is.finite(est) & abs(est) > tol])
    })
    li[ov_all[!ov_all %in% names(li)]] <- list(NULL)  # ensure all items present
    li[ov_all]
  })

  # --- Optimization 3): Skip per-group lookup when structure doesn't vary ----
  # Compare group-specific sets with overall sets; if all identical, skip group lookup
  set_equal <- function(a, b) {
    if (length(a) != length(b)) return(FALSE)
    all(sort(a) == sort(b))
  }
  vary_by_group <- any(vapply(ov_all, function(j) {
    overall <- loads_by_item_overall[[j]]
    if (is.null(overall)) return(TRUE)  # be conservative
    any(!vapply(loads_by_item_by_group, function(li) {
      set_equal(li[[j]] %||% character(0L), overall %||% character(0L))
    }, logical(1L)))
  }, logical(1L)))

  # -- mgcv setup --------------------------------------------------------------
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

  # -- Parallel plan -----------------------------------------------------------
  reset_plan <- .set_future_plan(plan = plan, workers = workers, cluster = cluster)
  on.exit(reset_plan(), add = TRUE)

  # -- 2) Empirical GAM curves (parallel over items within each group) --------
  for (g in seq_len(n_groups)) {
    rows_g <- if (".gid" %in% names(out)) which(out$.gid == g) else seq_len(nrow(out))
    if (!length(rows_g)) next

    eta_cols <- intersect(lats, names(out))
    if (!length(eta_cols)) stop("Factor scores missing in data; cannot build empirical curves.")

    # Freeze per-group data to keep futures light
    out_g <- out[rows_g, c(ov_all, eta_cols), drop = FALSE]

    # --- Optimization 4): cache formulas per unique set of predictors --------
    fml_cache <- new.env(parent = emptyenv())  # key: paste(rel_lats, collapse="|") -> formula

    item_results <- furrr::future_map(
      ov_all,
      function(j) {
        # -- choose relevant latents (group-specific -> overall -> all FS)
        rel_lats <- if (vary_by_group) {
          li_g <- loads_by_item_by_group[[as.character(g)]] %||% loads_by_item_by_group[[g]]
          (if (!is.null(li_g)) li_g[[j]] else NULL) %||% loads_by_item_overall[[j]]
        } else {
          loads_by_item_overall[[j]]
        }
        rel_lats <- intersect(rel_lats %||% eta_cols, eta_cols)
        if (!length(rel_lats) || !j %in% names(out_g)) return(NULL)

        df_g <- out_g[, c(j, rel_lats), drop = FALSE]
        fam_j   <- fam_for_item(j)
        invlink <- inv_for_item(j)
        rj      <- rng_by_item[[j]]

        # Build formula once per unique set of predictors
        key <- paste(rel_lats, collapse = "|")
        fml <- get0(key, envir = fml_cache, inherits = FALSE)
        if (is.null(fml)) {
          sm_terms <- paste(sprintf("s(%s)", rel_lats), collapse = " + ")
          fml <- stats::as.formula(paste("y_resp ~", sm_terms))
          assign(key, fml, envir = fml_cache)
        }

        # Fit data (complete cases only)
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

        # Prepare response ONCE (avoid repeated SV/rescale work)
        if (j %in% ov_cont) {
          y_fit <- df_ok[[j]]
          back_to_original <- function(mu_link) invlink(mu_link)
        } else {
          y_fit <- sv01(rescale_01(df_ok[[j]]))
          to_original_scale <- function(p01) { d <- rj[2] - rj[1]; p01 * d + rj[1] }
          back_to_original  <- function(mu_link) to_original_scale(invlink(mu_link))
        }
        X_fit <- df_ok[, rel_lats, drop = FALSE]

        if (nrow(X_fit) >= length(rel_lats) + 2L) {
          # Pick user args per item type
          user_args <- if (j %in% ov_cont) gam_args_cont else gam_args_ord

          fit_gam <- build_gam_call(
            fml,
            dat = cbind(y_resp = y_fit, X_fit),
            fam = fam_j,
            user_args = user_args
          )

          pr <- predict(
            fit_gam,
            newdata = df_g[, rel_lats, drop = FALSE],  # only predictors needed
            type = "link",
            se.fit = TRUE
          )

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
        packages = c("mgcv","stats")
      ),
      .progress = progress
    )

    # Write back results sequentially
    for (res in item_results) {
      if (is.null(res)) next
      out[rows_g, res$nm_est] <- res$mu
      out[rows_g, res$nm_lwr] <- res$lwr
      out[rows_g, res$nm_upr] <- res$upr
    }
  }

  # -- Type cleanup: keep .rid/.gid integer, everything numeric as double ------
  is_num <- vapply(out, is.numeric, logical(1L))
  keep_int <- names(out) %in% c(".rid", ".gid")
  to_double <- is_num & !keep_int
  for (nm in names(out)[to_double]) out[[nm]] <- as.double(out[[nm]])

  tibble::as_tibble(out)
}
