#' Augment ordinal SEM data with latent predictions, probabilities, and residuals
#'
#' Adds per-indicator latent predictions for ordinal outcomes (expected scores; `yhat`),
#' optionally their delta-method confidence intervals, per-category probabilities,
#' latent linear predictors (`y*`), and residuals defined as observed category (numeric)
#' minus expected score. Works for single- and multi-group lavaan models with **only ordinal indicators**.
#'
#' @param fit A fitted `lavaan` object for ordinal indicators.
#' @param ystar Logical; append latent linear predictors y* (default TRUE).
#' @param yhat Logical; append expected scores E[Y] (default TRUE).
#' @param pr Logical; append per-category probabilities (default TRUE).
#' @param ci Logical; add delta-method CIs for E[Y] (default TRUE).
#' @param level Confidence level for CIs (default 0.95).
#' @param resid Logical; append residuals (obs - E[Y]) in numeric scoring (default TRUE).
#' @param se_yhat Logical; append delta-method SEs of E[Y] (default TRUE).
#' @param prefix_ystar,prefix_yhat,prefix_pr Character prefixes for y*, E[Y], and probabilities.
#' @param prefix_ci Length-2 character vector with prefixes for lower/upper CI columns for E[Y]
#'   (default c(".yhat_lwr_", ".yhat_upr_")).
#' @param prefix_resid Prefix for residual columns (default ".resid_").
#' @param prefix_se_yhat Prefix for SE(E[Y]) columns (default ".se_yhat_").
#' @param sep Separator in probability column names between level and item (default "__").
#'
#' @return A tibble with factor scores, `.gid` and group label, and selected outputs.
#' @export
augment_ordinal <- function(fit,
                            ystar           = TRUE,
                            yhat            = TRUE,
                            pr              = TRUE,
                            ci              = TRUE,
                            level           = 0.95,
                            resid           = TRUE,
                            se_yhat         = TRUE,
                            prefix_ystar    = ".ystar_",
                            prefix_yhat     = ".yhat_",
                            prefix_pr       = ".pr_",
                            prefix_ci       = c(".yhat_lwr_", ".yhat_upr_"),
                            prefix_resid    = ".resid_",
                            prefix_se_yhat  = ".se_yhat_",
                            sep             = "__") {

  # -- Assertions --------------------------------------------------------------
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

  # -- Delta-method prerequisites (used if ci or se) ---------------------------
  coef_names <- names(lavaan::coef(fit))
  V_full     <- if (ci || se_yhat) lavaan::lavInspect(fit, "vcov") else NULL
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
      if (is.factor(yo)) return(levels(yo))
      lev <- sort(unique(yo)); return(as.character(lev))
    } else {
      tau_j <- .get_tau_for(th_vec, j)
      if (!length(tau_j)) return(character(0))
      m <- length(tau_j) + 1L
      as.character(seq_len(m))
    }
  }
  .needs_double_from_special <- function(x) inherits(x, "lvn.mtrx") || inherits(x, "lavaan.matrix")

  # -- Accumulators ------------------------------------------------------------
  ystar_list   <- if (ystar) vector("list", n_groups) else NULL
  probs_list   <- if (pr)    vector("list", n_groups) else NULL
  yexp_ci_list <- if (yhat || resid) vector("list", n_groups) else NULL

  # latent scores per group
  eta_list <- lapply(fs_list, function(df) as.matrix(df[, info$latent_variables, drop = FALSE]))

  for (g in seq_len(n_groups)) {
    est_g <- est_list[[g]]
    lam_g <- est_g$lambda
    nu_g  <- est_g$nu %||% numeric(0)

    al <- .align_Lambda_nu(lam_g, nu_g, ov_subset = ov_ord, ov_all = ov_all, latent_order = info$latent_variables)

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

    if (ystar) {
      ystar_df <- tibble::as_tibble(ystar_g)
      names(ystar_df) <- paste0(prefix_ystar, names(ystar_df))
      ystar_list[[g]] <- ystar_df
    }

    if (!(pr || yhat || resid)) {
      if (!is.null(yexp_ci_list)) yexp_ci_list[[g]] <- tibble::tibble()
      next
    }

    yexp_ci_cols <- list()
    probs_g      <- if (pr) NULL else FALSE  # sentinel

    lam_cols <- if (!is.null(colnames(al$lambda))) colnames(al$lambda) else info$latent_variables

    for (jj in seq_along(ov_ord)) {
      j       <- ov_ord[jj]
      ystar_j <- ystar_g[, jj]

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

      Phi     <- stats::pnorm(zmat)
      probs_j <- Phi[, -1, drop = FALSE] - Phi[, -ncol(Phi), drop = FALSE]  # N x m
      if (pr) {
        colnames(probs_j) <- paste0(prefix_pr, make.names(lev, allow_ = TRUE), sep, j)
        probs_g <- if (is.null(probs_g)) as.data.frame(probs_j, check.names = FALSE) else
          cbind(probs_g, as.data.frame(probs_j, check.names = FALSE))
      }

      if (yhat || resid) {
        # Expected scores E[Y]
        scores_j <- { sc <- suppressWarnings(as.numeric(lev)); if (anyNA(sc)) seq_len(m) else sc }
        yexp_j   <- as.numeric(probs_j %*% scores_j)  # N-vector

        # Container for this item
        col_core <- tibble::tibble(.rows = length(yexp_j))

        # Do we need delta-method pieces? (for CI and/or SE of yhat)
        need_dm <- (ci && yhat) || (se_yhat && yhat)

        if (need_dm) {
          # Derivative building blocks
          phi_mat    <- stats::dnorm(zmat)                      # N x (m+2)
          phi_inner  <- phi_mat[, -1, drop = FALSE]
          phi_prev   <- phi_mat[, -ncol(phi_mat), drop = FALSE]

          # B term used in d/d(nu,lambda)
          B <- as.numeric((phi_inner[, 1:m, drop = FALSE] - phi_prev[, 1:m, drop = FALSE]) %*% matrix(scores_j, ncol = 1))

          # dmu/d nu_j
          if (identical(param, "theta")) {
            se_j <- if (is.null(sig_eps)) 1 else sig_eps[jj]
            d_nu <- (-1 / se_j) * B
          } else {
            d_nu <- (-1) * B
          }

          # dmu/d lambda_jr
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

          # Build gradient G and SE(yhat)
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
            se_vec <- sqrt(rowSums(GV * G))
          } else {
            se_vec <- rep(NA_real_, length(yexp_j))
          }

          # yhat column (if requested)
          if (yhat) {
            col_core[[paste0(prefix_yhat, j)]] <- yexp_j
            if (se_yhat) {
              col_core[[paste0(prefix_se_yhat, j)]] <- se_vec
            }
            if (ci) {
              col_core[[paste0(prefix_ci[[1]], j)]] <- yexp_j - zcrit * se_vec
              col_core[[paste0(prefix_ci[[2]], j)]] <- yexp_j + zcrit * se_vec
            }
          }

        } else if (yhat) {
          # No delta-method -> only E[Y]
          col_core[[paste0(prefix_yhat, j)]] <- yexp_j
        }

        # Residuals
        if (resid) {
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

    if (pr) {
      probs_list[[g]] <- if (is.null(probs_g)) tibble::tibble() else tibble::as_tibble(probs_g)
    }
    if (yhat || resid) {
      yexp_ci_list[[g]] <- if (length(yexp_ci_cols)) dplyr::bind_cols(yexp_ci_cols) else tibble::tibble()
    }
  } # end groups

  # -- Final bind: FS | (ystar?) | (yhat/SE/CI/resid?) | (pr?) -----------------
  fs_tbl_parts <- vector("list", n_groups)
  for (g in seq_len(n_groups)) {
    fs_g <- tibble::as_tibble(fs_list[[g]])
    fs_g <- tibble::add_column(fs_g, .gid = g, .before = 1)
    # create correct group label column named `.group`
    fs_g <- tibble::add_column(fs_g, .group = group_labels[g], .after = ".gid")
    fs_tbl_parts[[g]] <- fs_g
  }
  fs_tbl <- dplyr::bind_rows(fs_tbl_parts)

  # Relocate factor-score columns immediately after `.group`
  eta_cols <- info$latent_variables
  eta_cols <- eta_cols[eta_cols %in% names(fs_tbl)]  # be robust if some scores are absent
  if (length(eta_cols)) {
    fs_tbl <- dplyr::relocate(fs_tbl, dplyr::all_of(eta_cols), .after = ".group")
  }

  parts <- list(fs_tbl)
  if (ystar)         parts <- c(parts, list(dplyr::bind_rows(ystar_list)))
  if (yhat || resid) parts <- c(parts, list(dplyr::bind_rows(yexp_ci_list)))
  if (pr)            parts <- c(parts, list(dplyr::bind_rows(probs_list)))

  out <- do.call(dplyr::bind_cols, parts)

  # Coerce special classes -> double
  out <- dplyr::mutate(out, dplyr::across(.cols = where(.needs_double_from_special), .fns = ~ as.double(.x)))

  tibble::as_tibble(out)
}
