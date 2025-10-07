#' Augment ordinal SEM data with latent predictions, probabilities, and residuals
#'
#' @description
#' Adds per-indicator latent predictions for ordinal outcomes (expected scores; `yhat`),
#' optionally their delta-method confidence intervals, per-category probabilities,
#' latent linear predictors (`y*`), and residuals defined as *observed category (numeric) minus expected score*.
#' Works for single- and multi-group `lavaan` models with **only ordinal indicators**.
#'
#' @details
#' For each ordinal indicator \eqn{j}, the function:
#' \itemize{
#'   \item computes the latent linear predictor \eqn{y^*_j = \nu_j + \lambda_{j\cdot}\eta} (from factor scores and measurement parameters),
#'   \item obtains category probabilities \eqn{P(Y_j = k)} using probit links implied by thresholds,
#'   \item returns the expected score \eqn{E[Y_j]} on the numeric scoring of categories
#'         (uses actual numeric labels if available; otherwise scores \eqn{1, \ldots, m}),
#'   \item optionally returns delta-method CIs for \eqn{E[Y_j]} by propagating the sampling covariance
#'         of thresholds, intercept/mean, loadings (and residual variance under \code{theta} parametrization),
#'   \item and optionally returns residuals \eqn{R_j = \text{obs}_j - \widehat{E[Y_j]}} mapped to the same numeric scale.
#' }
#'
#' Multi-group models are supported; outputs are row-bound across groups and include
#' group identifiers. This function expects that the fitted model is intended for
#' ordered-categorical indicators and relies on \code{lavaan::lavInspect(..., "est")}
#' and factor-score predictions from a parallel-safe helper (\code{lavPredict_parallel}).
#'
#' @section Parameterization:
#' The routine respects the model's parameterization reported by \code{model_info(fit)$parameterization}.
#' Under the \code{"theta"} parameterization, residual SDs enter the probit z-scores and the delta-method
#' gradient includes the residual variance term; under \code{"delta"} they do not.
#'
#' @param fit A fitted \code{lavaan} object for ordinal indicators.
#' @param ystar Logical. If \code{TRUE}, append latent linear predictors \eqn{y^*} for each ordinal indicator (prefixed by \code{prefix_ystar}).
#' @param yhat Logical. If \code{TRUE} (default), append expected scores \eqn{E[Y]} for each ordinal indicator (prefixed by \code{prefix_yhat}).
#' @param pr Logical. If \code{TRUE}, append per-category probabilities for each ordinal indicator
#'   (columns named \code{<prefix_pr><level><sep><item>} where \code{<level>} are observed or implied levels).
#' @param ci Logical. If \code{TRUE}, compute delta-method confidence intervals for \eqn{E[Y]} and append
#'   lower/upper columns using \code{prefix_ci}.
#' @param resid Logical. If \code{TRUE}, append residuals \eqn{R = \text{observed} - \widehat{E[Y]}}
#'   for each ordinal indicator (prefixed by \code{prefix_resid}). Residuals are computed on the numeric
#'   scoring used for \eqn{E[Y]} (actual numeric labels if available, otherwise \code{1..m}).
#' @param level Confidence level for delta-method intervals; default \code{0.95}.
#' @param prefix_ystar,prefix_yhat,prefix_pr Character scalars with prefixes for the corresponding outputs.
#' @param prefix_ci Character vector of length 2 with prefixes for lower/upper CI columns for \eqn{E[Y]}
#'   (default \code{c(".yhat_lwr_", ".yhat_upr_")}).
#' @param prefix_resid Character scalar; prefix for residual columns (default \code{".resid"}).
#' @param sep Character separator used when forming probability column names between level and item (default \code{"__"}).
#'
#' @return A tibble with:
#' \itemize{
#'   \item leading factor-score columns (from \code{lavPredict_parallel}),
#'   \item \code{.gid} (group index) and \code{group} (group label),
#'   \item optional blocks according to toggles: \code{ystar}, \code{yhat} (\& optional CIs), \code{pr}, and \code{resid}.
#' }
#' Column naming:
#' \itemize{
#'   \item \code{<prefix_ystar><item>} for \eqn{y^*},
#'   \item \code{<prefix_yhat><item>} for \eqn{E[Y]},
#'   \item \code{<prefix_ci[1]><item>} and \code{<prefix_ci[2]><item>} for CI bounds,
#'   \item \code{<prefix_pr><level><sep><item>} for probabilities,
#'   \item \code{<prefix_resid><item>} for residuals.
#' }
#'
#' @section Assumptions and notes:
#' \itemize{
#'   \item Only ordinal indicators are supported (\code{model_info(fit)$ov_ordinal}).
#'   \item Expected scores and residuals use numeric scoring of categories. If observed categories are non-numeric
#'         (e.g., strings), they are mapped to \code{1..m} in the order of levels.
#'   \item Delta-method CIs require a valid parameter covariance matrix (\code{lavInspect(fit, "vcov")}).
#'   \item The function relies on helpers present in the package: \code{.assert_lavaan_fit}, \code{model_info},
#'         and \code{lavPredict_parallel}.
#' }
#'
#' @examples
#' \dontrun{
#' library(lavaan)
#' # Suppose `fit` is an ordinal CFA model
#' aug <- augment_ordinal_5(
#'   fit,
#'   yhat        = TRUE,
#'   ci          = TRUE,
#'   resid       = TRUE,
#'   pr          = FALSE,
#'   prefix_yhat = ".yhat_",
#'   prefix_resid= ".resid"
#' )
#' dplyr::glimpse(aug)
#' }
#'
#' @seealso \code{\link{lavPredict_parallel}}, \code{\link{model_info}}, \code{\link[lavaan]{lavInspect}}
#' @export


augment_ordinal_old <- function(fit,
                            ystar           = TRUE,
                            yhat            = TRUE,
                            pr              = TRUE,
                            ci              = TRUE,
                            level           = 0.95,
                            resid           = TRUE,
                            prefix_ystar    = ".ystar_",
                            prefix_yhat     = ".yhat_",
                            prefix_pr       = ".pr_",
                            prefix_ci       = c(".yhat_lwr_", ".yhat_upr_"),
                            prefix_resid    = ".resid_",
                            sep             = "__"

) {

  # -- Strict assertions -------------------------------------------------------
  # requires: converged fit, latent variables, and measurement (lambda + nu)
  .assert_lavaan_fit(
    fit,
    require_converged     = TRUE,
    require_latent        = TRUE,
    require_meanstructure = NA,
    require_measurement   = "lambda+nu",
    forbid_multilevel     = FALSE
  )

  # -- Model info --------------------------------------------------------------
  info         <- model_info(fit)
  ov_all       <- info$observed_variables
  ov_ord       <- info$ov_ordinal
  n_groups     <- info$n_groups %||% 1L
  group_labels <- info$group_labels %||% as.character(seq_len(n_groups))
  param        <- info$parameterization %||% "delta"

  fs_and_ov <- lavPredict_parallel(fit, return_type = "list")
  fs_list   <- if (is.data.frame(fs_and_ov)) list(fs_and_ov) else fs_and_ov

  est      <- lavaan::lavInspect(fit, "est")
  est_list <- if (is.list(est) && "lambda" %in% names(est)) list(est) else est

  th_raw   <- tryCatch(lavaan::lavInspect(fit, "th"), error = function(e) NULL)
  th_list  <- if (is.null(th_raw)) replicate(n_groups, setNames(numeric(0), character(0)), simplify = FALSE)
  else if (is.list(th_raw)) th_raw else list(th_raw)
  if (length(th_list) == 1L && n_groups > 1L) th_list <- rep(th_list, n_groups)
  names(th_list) <- group_labels

  # -- Cache for delta method (used only if yhat && ci == TRUE) ----------------
  coef_names <- names(lavaan::coef(fit))
  V_full     <- lavaan::lavInspect(fit, "vcov")
  zcrit      <- stats::qnorm(1 - (1 - level)/2)

  # -- Helpers -----------------------------------------------------------------
  .align_Lambda_nu <- function(lam_g, nu_g, ov_subset, ov_all, latent_order) {
    if (!is.null(colnames(lam_g))) lam_g <- lam_g[, latent_order, drop = FALSE]
    if (!is.null(rownames(lam_g))) {
      lam_g <- lam_g[ov_subset, , drop = FALSE]
    } else {
      row_idx <- match(ov_subset, ov_all); lam_g <- lam_g[row_idx, , drop = FALSE]
    }
    if (is.null(nu_g) || length(nu_g) == 0L) {
      nu_g <- rep(0, length(ov_subset)); names(nu_g) <- ov_subset
    } else if (!is.null(names(nu_g))) {
      nu_g <- nu_g[ov_subset]
    } else {
      row_idx <- match(ov_subset, ov_all); nu_g <- nu_g[row_idx]
    }
    list(lambda = lam_g, nu = nu_g)
  }

  .get_tau_for <- function(th_vec, j) {
    nm <- names(th_vec); if (is.null(nm)) return(numeric(0))
    sel <- startsWith(nm, paste0(j, "|"))
    tau <- th_vec[sel]
    if (!length(tau)) return(numeric(0))
    idx <- suppressWarnings(as.integer(gsub("^.*\\|t(\\d+)$", "\\1", names(tau))))
    idx[is.na(idx)] <- seq_along(tau)
    tau[order(idx)]
  }

  .levels_for_item <- function(df, th_vec, j) {
    if (!is.null(df) && j %in% names(df)) {
      yo <- df[[j]]
      if (is.factor(yo)) {
        return(levels(yo))
      } else {
        lev <- sort(unique(yo))
        return(as.character(lev))
      }
    } else {
      tau_j <- .get_tau_for(th_vec, j)
      if (!length(tau_j)) return(character(0))
      m <- length(tau_j) + 1L
      return(as.character(seq_len(m)))
    }
  }

  .needs_double_from_special <- function(x) {
    inherits(x, "lvn.mtrx") || inherits(x, "lavaan.matrix")
  }

  # -- Accumulators (created only if needed) -----------------------------------
  ystar_list   <- if (isTRUE(ystar)) vector("list", n_groups) else NULL
  probs_list   <- if (isTRUE(pr))    vector("list", n_groups) else NULL
  yexp_ci_list <- if (isTRUE(yhat) || isTRUE(resid)) vector("list", n_groups) else NULL

  # latent scores per group
  eta_list <- lapply(fs_list, function(df) as.matrix(df[, info$latent_variables, drop = FALSE]))

  for (g in seq_len(n_groups)) {
    est_g <- est_list[[g]]
    lam_g <- est_g$lambda
    nu_g  <- est_g$nu %||% numeric(0)

    al <- .align_Lambda_nu(lam_g, nu_g,
                           ov_subset    = ov_ord,
                           ov_all       = ov_all,
                           latent_order = info$latent_variables)

    eta_g    <- eta_list[[g]]                                 # N x q
    ystar_g  <- sweep(eta_g %*% t(al$lambda), 2, al$nu, `+`)  # N x p_ord
    colnames(ystar_g) <- ov_ord

    # residual SDs in theta parametrization
    sig_eps <- NULL
    if (identical(param, "theta") && (pr || yhat || resid)) {
      theta_g <- est_g$theta
      if (!is.null(theta_g)) {
        if (!is.null(rownames(theta_g))) {
          sig_eps <- sqrt(diag(theta_g)[ov_ord])
        } else {
          idx <- match(ov_ord, ov_all); sig_eps <- sqrt(diag(theta_g)[idx])
        }
        sig_eps[!is.finite(sig_eps) | sig_eps <= 0] <- 1
      } else sig_eps <- rep(1, length(ov_ord))
    }

    # Store y* only if requested
    if (isTRUE(ystar)) {
      ystar_df <- tibble::as_tibble(ystar_g)
      names(ystar_df) <- paste0(prefix_ystar, names(ystar_df))
      ystar_list[[g]] <- ystar_df
    }

    # If neither pr, yhat, nor resid is requested, we can skip the per-item loop entirely
    if (!(pr || yhat || resid)) {
      if (!is.null(yexp_ci_list)) yexp_ci_list[[g]] <- tibble::tibble()
      next
    }

    # interleaved yhat (+/- CI) and/or residuals accumulator for this group
    yexp_ci_cols <- list()
    probs_g      <- if (isTRUE(pr)) NULL else FALSE  # FALSE -> sentinel = "do not keep"

    # factor names aligned to lambda columns
    lam_cols <- if (!is.null(colnames(al$lambda))) colnames(al$lambda) else info$latent_variables

    for (jj in seq_along(ov_ord)) {
      j       <- ov_ord[jj]
      ystar_j <- ystar_g[, jj]

      # levels & thresholds; required if pr, yhat, or resid
      lev   <- .levels_for_item(fs_list[[g]], th_list[[g]], j)
      tau_j <- as.numeric(.get_tau_for(th_list[[g]], j))
      if (length(tau_j) && length(lev) != (length(tau_j) + 1L)) lev <- as.character(seq_len(length(tau_j) + 1L))
      if (!length(lev) || !length(tau_j)) next

      m       <- length(lev)
      tau_pad <- c(-Inf, tau_j, Inf)

      # z-matrix: N x (m+2)
      if (identical(param, "theta")) {
        se_j <- if (is.null(sig_eps)) 1 else sig_eps[jj]
        zmat <- outer(ystar_j, tau_pad, function(y, t) (t - y) / se_j)
      } else {
        zmat <- outer(ystar_j, tau_pad, function(y, t) (t - y))
      }

      # Probabilities (needed for yhat and residuals; optionally kept if pr = TRUE)
      Phi     <- stats::pnorm(zmat)
      probs_j <- Phi[, -1, drop = FALSE] - Phi[, -ncol(Phi), drop = FALSE]  # N x m
      if (isTRUE(pr)) {
        colnames(probs_j) <- paste0(prefix_pr, make.names(lev, allow_ = TRUE), sep, j)
        probs_g <- if (is.null(probs_g)) as.data.frame(probs_j, check.names = FALSE) else
          cbind(probs_g, as.data.frame(probs_j, check.names = FALSE))
      }

      # Expected scores (yhat) and/or residuals
      if (isTRUE(yhat) || isTRUE(resid)) {
        # Map category labels to numeric scores (use actual numeric levels if possible)
        scores_j <- { sc <- suppressWarnings(as.numeric(lev)); if (anyNA(sc)) seq_len(m) else sc }
        yexp_j   <- as.numeric(probs_j %*% scores_j)  # N-vector

        # Start column container for this indicator
        col_core <- tibble::tibble(.rows = length(yexp_j))

        # (Optional) confidence intervals for yhat
        if (isTRUE(yhat) && ci == TRUE) {
          # --- Vectorized gradients (unchanged logic) ------------------------
          phi_mat    <- stats::dnorm(zmat)                  # N x (m+2)
          phi_inner  <- phi_mat[, -1, drop = FALSE]         # z_k
          phi_prev   <- phi_mat[, -ncol(phi_mat), drop = FALSE]  # z_{k-1}

          # Base term B = sum_k s_k * (phi(z_k) - phi(z_{k-1}))  -> N-vector
          B <- as.numeric((phi_inner[, 1:m, drop = FALSE] - phi_prev[, 1:m, drop = FALSE]) %*% matrix(scores_j, ncol = 1))

          # dmu/d nu_j
          if (identical(param, "theta")) {
            se_j <- if (is.null(sig_eps)) 1 else sig_eps[jj]
            d_nu <- (-1 / se_j) * B
          } else {
            d_nu <- (-1) * B
          }

          # dmu/d lambda_jr for each factor r
          eta_mat    <- eta_g[, lam_cols, drop = FALSE]
          if (identical(param, "theta")) {
            dz_lam_mat <- -eta_mat / (if (is.null(sig_eps)) 1 else sig_eps[jj])
          } else {
            dz_lam_mat <- -eta_mat
          }
          d_lam_mat <- dz_lam_mat * B

          # dmu/d tau_{jr}, r=1..m-1
          s_diff    <- scores_j[1:(m - 1)] - scores_j[2:m]
          phi_r_mat <- phi_mat[, 1 + seq_len(m - 1), drop = FALSE]
          if (identical(param, "theta")) {
            se_j <- if (is.null(sig_eps)) 1 else sig_eps[jj]
            d_tau_mat <- sweep(phi_r_mat, 2, (s_diff / se_j), `*`)
          } else {
            d_tau_mat <- sweep(phi_r_mat, 2, s_diff, `*`)
          }

          # dmu/d theta_jj (theta only)
          if (identical(param, "theta")) {
            se_j_local <- if (is.null(sig_eps)) 1 else sig_eps[jj]
            dz_sigma   <- -(zmat / se_j_local)
            A <- phi_inner * dz_sigma[, -1, drop = FALSE]
            C <- phi_prev  * dz_sigma[, -ncol(dz_sigma), drop = FALSE]
            d_sigma <- as.numeric((A[, 1:m, drop = FALSE] - C[, 1:m, drop = FALSE]) %*% matrix(scores_j, ncol = 1))
            d_var   <- d_sigma * (1 / (2 * se_j_local))
          }

          # Build gradient G and SEs
          th_names  <- paste0(j, "|t", seq_len(m - 1))
          nu_name   <- { nm <- paste0("nu[", j, "]"); if (nm %in% coef_names) nm else {
            alt <- paste0(j, "~1"); if (alt %in% coef_names) alt else NA_character_
          }}
          lam_names <- paste0(j, "=~", lam_cols)
          var_name  <- if (identical(param, "theta")) paste0(j, "~~", j) else NA_character_

          keep <- c(th_names, nu_name, lam_names, var_name)
          keep <- keep[!is.na(keep) & keep %in% coef_names]
          idx  <- match(keep, coef_names)
          Vsub <- V_full[idx, idx, drop = FALSE]

          G_parts <- list()
          if (length(th_names)) {
            sel <- th_names %in% keep
            if (any(sel)) G_parts[["tau"]] <- d_tau_mat[, sel, drop = FALSE]
          }
          if (!is.na(nu_name) && nu_name %in% keep) {
            G_parts[["nu"]] <- matrix(d_nu, ncol = 1); colnames(G_parts[["nu"]]) <- nu_name
          }
          lam_keep <- lam_names[lam_names %in% keep]
          if (length(lam_keep)) {
            sel <- lam_cols %in% sub("^.*=~", "", lam_keep)
            G_parts[["lam"]] <- d_lam_mat[, sel, drop = FALSE]; colnames(G_parts[["lam"]]) <- lam_keep
          }
          if (identical(param, "theta") && !is.na(var_name) && var_name %in% keep) {
            G_parts[["var"]] <- matrix(d_var, ncol = 1); colnames(G_parts[["var"]]) <- var_name
          }

          if (length(G_parts)) {
            G  <- do.call(cbind, G_parts)
            GV <- G %*% Vsub
            se <- sqrt(rowSums(GV * G))
          } else {
            se <- rep(NA_real_, length(yexp_j))
          }

          yexp_lwr_j <- yexp_j - zcrit * se
          yexp_upr_j <- yexp_j + zcrit * se

          # Add yhat + CI columns
          col_core[[paste0(prefix_yhat, j)]]    <- yexp_j
          col_core[[paste0(prefix_ci[[1]], j)]] <- yexp_lwr_j
          col_core[[paste0(prefix_ci[[2]], j)]] <- yexp_upr_j

        } else if (isTRUE(yhat)) {
          # No CI -> only expected score column
          col_core[[paste0(prefix_yhat, j)]] <- yexp_j
        }

        # (Optional) residuals: observed minus expected (mapped to same numeric scale)
        if (isTRUE(resid)) {
          yobs <- fs_list[[g]][[j]]
          if (is.factor(yobs)) {
            idx_obs <- match(as.character(yobs), lev)
            yobs_num <- idx_obs
          } else {
            idx_obs <- match(as.character(yobs), lev)
            if (anyNA(idx_obs)) {
              idx2 <- match(suppressWarnings(as.numeric(yobs)), suppressWarnings(as.numeric(lev)))
              idx_obs[is.na(idx_obs)] <- idx2[is.na(idx_obs)]
            }
            yobs_num <- idx_obs
          }
          resid_j <- as.numeric(yobs_num) - yexp_j
          col_core[[paste0(prefix_resid, j)]] <- resid_j
        }

        yexp_ci_cols[[length(yexp_ci_cols) + 1L]] <- col_core
      } # end if (yhat || resid)
    } # end items loop

    if (isTRUE(pr)) {
      probs_list[[g]] <- if (is.null(probs_g)) tibble::tibble() else tibble::as_tibble(probs_g)
    }
    if (isTRUE(yhat) || isTRUE(resid)) {
      yexp_ci_list[[g]] <- if (length(yexp_ci_cols)) dplyr::bind_cols(yexp_ci_cols) else tibble::tibble()
    }
  } # end groups

  # -- Final bind: FS | (ystar?) | (yhat/CI or resid?) | (pr?) -----------------
  fs_tbl_parts <- vector("list", n_groups)
  for (g in seq_len(n_groups)) {
    fs_g <- tibble::as_tibble(fs_list[[g]])
    fs_g <- tibble::add_column(fs_g, .gid = g, .before = 1)
    fs_g <- tibble::add_column(fs_g, group = group_labels[g], .after = ".gid")
    fs_tbl_parts[[g]] <- fs_g
  }
  fs_tbl <- dplyr::bind_rows(fs_tbl_parts)

  parts <- list(fs_tbl)
  if (isTRUE(ystar))                 parts <- c(parts, list(dplyr::bind_rows(ystar_list)))
  if (isTRUE(yhat) || isTRUE(resid)) parts <- c(parts, list(dplyr::bind_rows(yexp_ci_list)))
  if (isTRUE(pr))                    parts <- c(parts, list(dplyr::bind_rows(probs_list)))

  out <- do.call(dplyr::bind_cols, parts)

  # Coerce special classes -> double (keep tibbles clean)
  out <- dplyr::mutate(out, dplyr::across(.cols = where(.needs_double_from_special), .fns = ~ as.double(.x)))

  tibble::as_tibble(out)
}
