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
prepare_original_old <- function(fit,
                             data     = NULL,
                             info     = NULL,
                             level    = 0.95,
                             fam_cont = stats::gaussian(),
                             fam_ord  = mgcv::betar(link = "logit"),
                             ...) {

  # --- Assertions -------------------------------------------------------------
  .assert_lavaan_fit(fit, require_converged = TRUE, require_latent = TRUE, forbid_multilevel = TRUE)
  if (is.null(info)) info <- model_info(fit)

  ov_all   <- info$observed_variables
  ov_cont  <- info$ov_continuous
  ov_ord   <- info$ov_ordinal
  lats     <- info$latent_variables
  n_groups <- info$n_groups %||% 1L

  # --- 1) MODEL-BASED curves via augment() with requested prefixes -----------
  out <- augment(
    fit            = fit,
    data           = data,
    info           = info,
    # keep only essentials
    yhat           = TRUE,
    ci             = TRUE,
    resid          = FALSE,
    se_yhat        = FALSE,
    ystar          = FALSE,
    pr             = FALSE,
    # FS SEs for continuous-only models
    se_fs          = TRUE,
    # column layout unchanged; we only change prefixes:
    prefix_yhat    = "m_est_yhat_",
    prefix_ci      = c("m_lwr_yhat_", "m_upr_yhat_"),
    # keep the rest as defaults
    col_layout     = "by_type"
  )

  # --- Helpers for empirical curves ------------------------------------------
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

  # --- Robust extraction of nonzero loadings across groups ---------------------
  PE_raw <- lavaan::parameterEstimates(fit, standardized = FALSE)

  # Keep only loadings (op == "=~") for items present in the model
  PE_sel <- PE_raw[PE_raw$op == "=~" & PE_raw$rhs %in% ov_all, , drop = FALSE]

  # Safe column selection (some lavaan builds omit `group` for single-group)
  keep_cols <- intersect(c("group","lhs","rhs","est"), names(PE_sel))
  PE <- PE_sel[, keep_cols, drop = FALSE]

  # If `group` is absent (single-group), add it
  if (!"group" %in% names(PE)) PE$group <- 1L

  # Build per-item list of relevant factors (nonzero loadings)
  tol <- 1e-10
  loads_by_item <- tapply(seq_len(nrow(PE)), PE$rhs, function(idx) {
    lhs <- PE$lhs[idx]
    est <- PE$est[idx]
    unique(lhs[is.finite(est) & abs(est) > tol])
  })
  # Preserve model order
  loads_by_item <- loads_by_item[ov_all]

  # --- 2) EMPIRICAL curves via mgcv::gam() per item & group -------------------
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is required for empirical curves.")
  }
  zcrit <- stats::qnorm(1 - (1 - level) / 2)

  # Choose which family to use per item (continuous vs ordinal)
  fam_for_item <- function(j) if (j %in% ov_cont) fam_cont else fam_ord
  inv_for_item <- function(j) linkinv_fun(fam_for_item(j))

  # For ordinal items, we'll return to original scale [min,max] (usually 1..K)
  rng_by_item <- lapply(ov_all, function(j) {
    if (!j %in% names(out)) return(c(0,1))
    x <- out[[j]]
    if (is.factor(x)) c(1, nlevels(x)) else range(as.numeric(x), na.rm = TRUE, finite = TRUE)
  })
  names(rng_by_item) <- ov_all

  # Per-group loop (fit separate GAMs to allow group-specific smooths)
  for (g in seq_len(n_groups)) {
    rows_g <- if (".gid" %in% names(out)) which(out$.gid == g) else seq_len(nrow(out))
    if (!length(rows_g)) next

    # factor scores available?
    eta_cols <- intersect(lats, names(out))
    if (!length(eta_cols)) stop("Factor scores missing in data; cannot build empirical curves.")

    # Build per-item empirical curves
    for (j in ov_all) {
      if (!j %in% names(out)) next  # skip if not present in lavPredict data

      # Relevant factors for item j (fall back to all lats if unknown)
      rel_lats <- loads_by_item[[j]]
      if (is.null(rel_lats) || !length(rel_lats)) rel_lats <- eta_cols
      rel_lats <- rel_lats[rel_lats %in% eta_cols]
      if (!length(rel_lats)) next

      # Data for group g
      df_g <- out[rows_g, c(j, rel_lats), drop = FALSE]

      # Set up response & family-dependent preparation
      fam_j   <- fam_for_item(j)
      invlink <- inv_for_item(j)

      # Response for GAM:
      # - continuous: raw yj
      # - ordinal: rescaled to (0,1) with SV adjustment
      if (j %in% ov_cont) {
        y_resp <- df_g[[j]]
        # For gaussian identity, no transform needed; for others, we still fit on link scale and inverse later
        back_to_original <- function(mu_link) invlink(mu_link)  # likely identity
        # No extra range mapping needed
        to_original_scale <- identity
      } else {
        # ordinal
        rj   <- rng_by_item[[j]]
        y01  <- rescale_01(df_g[[j]])
        y01  <- sv01(y01)  # (0,1) for betar
        y_resp <- y01

        # Inverse link gives predicted proportion in (0,1); map back to [min,max]
        to_original_scale <- function(p01) {
          d <- rj[2] - rj[1]
          p01 * d + rj[1]
        }
        back_to_original <- function(mu_link) to_original_scale(invlink(mu_link))
      }

      # Build formula: additive smooths for each relevant factor
      sm_terms <- paste0("s(", rel_lats, ")", collapse = " + ")
      fml <- stats::as.formula(paste0("y_resp ~ ", sm_terms))

      # Fit GAM (NA handling: drop rows with NA in y or predictors)
      ok <- stats::complete.cases(df_g[, c(j, rel_lats), drop = FALSE])
      df_fit <- df_g[ok, , drop = FALSE]
      y_fit  <- if (j %in% ov_cont) df_fit[[j]] else sv01(rescale_01(df_fit[[j]]))
      X_fit  <- df_fit[, rel_lats, drop = FALSE]

      if (nrow(X_fit) >= length(rel_lats) + 2L) {
        fit_gam <- mgcv::gam(fml,
                             data   = cbind(y_resp = y_fit, X_fit),
                             family = fam_j,
                             ...)
        pr <- predict(fit_gam,
                      newdata = df_g[, rel_lats, drop = FALSE],
                      type = "link",
                      se.fit = TRUE)

        # CIs on link scale -> inverse link -> (optionally) back to original scale
        mu  <- back_to_original(pr$fit)
        lwr <- back_to_original(pr$fit - zcrit * pr$se.fit)
        upr <- back_to_original(pr$fit + zcrit * pr$se.fit)

      } else {
        # Not enough data to fit GAM -> return NAs
        mu  <- rep(NA_real_, nrow(df_g))
        lwr <- upr <- mu
      }

      # Column names and assignment
      nm_est <- paste0("e_est_yhat_", j)
      nm_lwr <- paste0("e_lwr_yhat_", j)
      nm_upr <- paste0("e_upr_yhat_", j)

      out[rows_g, nm_est] <- as.numeric(mu)
      out[rows_g, nm_lwr] <- as.numeric(lwr)
      out[rows_g, nm_upr] <- as.numeric(upr)
    } # items
  }   # groups

  # Ensure numeric columns are double (keep .rid/.gid integer if present)
  is_num <- vapply(out, is.numeric, logical(1L))
  keep_int <- names(out) %in% c(".rid", ".gid")
  to_double <- is_num & !keep_int
  for (nm in names(out)[to_double]) out[[nm]] <- as.double(out[[nm]])

  tibble::as_tibble(out)
}
