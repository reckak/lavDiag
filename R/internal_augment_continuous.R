#' Add predicted indicator values, residuals, and (optionally) CIs for continuous indicators
#'
#' Computes \eqn{\hat{y} = \nu + \Lambda \eta} for observed **continuous**
#' indicators and returns a tibble with:
#'   - `.rid`, `.gid`, `.group`, factor scores (+ optional FS SEs),
#'   - predicted values (`prefix_yhat`),
#'   - residuals (`prefix_resid`),
#'   - optional delta-method SEs for \eqn{\hat y} (`prefix_se_yhat`) and CIs.
#'
#' Note: Factor-score SEs are available only for continuous-only models. If any
#' ordinal indicator is present, `se_fs = TRUE` is ignored (with a warning).
#'
#' @param fit A fitted `lavaan` object with `meanstructure = TRUE`.
#' @param data Optional: data frame (single-group) or list of data frames (multi-group)
#'   from `lavPredict_parallel(fit, return_type = "list"|"data")`. If `NULL`, computed internally.
#' @param info Optional: result of `model_info(fit)`. If `NULL`, computed internally.
#' @param yhat Logical; include predicted values \eqn{\hat y}. Default `TRUE`.
#' @param resid Logical; include residuals (observed - \eqn{\hat y}). Default `TRUE`.
#' @param ci Logical; add Wald CIs for \eqn{\hat y} via delta method. Default `TRUE`.
#' @param level Confidence level for CIs. Default `0.95`.
#' @param vcov_type Optional `lavaan::vcov()` type. If `NULL`, uses default.
#' @param se_fs Logical; request factor-score SEs (continuous-only). Default `TRUE`.
#' @param se_yhat Logical; add delta-method SEs for \eqn{\hat y}. Default `TRUE`.
#' @param prefix_yhat,prefix_resid Character prefixes for predictions and residuals.
#' @param prefix_ci Length-2 vector with prefixes for lower/upper CI columns.
#' @param prefix_se_fs Prefix for factor-score SE columns (e.g., ".se_").
#' @param prefix_se_yhat Prefix for SE(\eqn{\hat y}) columns (e.g., ".se_yhat_").
#'
#' @return A tibble with `.rid`, `.gid`, `.group`, factor scores (and optional FS SEs),
#'   observed variables, and selected augmentation outputs for continuous indicators.
#' @noRd
.augment_continuous <- function(fit,
                                data            = NULL,
                                info            = NULL,
                                yhat            = TRUE,
                                ci              = TRUE,
                                level           = 0.95,
                                resid           = TRUE,
                                se_fs           = TRUE,
                                se_yhat         = TRUE,
                                prefix_yhat     = ".yhat_",
                                prefix_ci       = c(".yhat_lwr_", ".yhat_upr_"),
                                prefix_resid    = ".resid_",
                                prefix_se_fs   = ".se_",
                                prefix_se_yhat  = ".se_yhat_",
                                vcov_type       = NULL) {

  # -- Validate fit ------------------------------------------------------------
  .assert_lavaan_fit(
    fit,
    require_converged     = TRUE,
    require_meanstructure = TRUE,
    require_latent        = TRUE,
    forbid_multilevel     = TRUE
  )

  # -- Validate prefix_ci ------------------------------------------------------
  if (!is.character(prefix_ci) || length(prefix_ci) != 2L || any(!nzchar(prefix_ci))) {
    stop("`prefix_ci` must be two non-empty strings, e.g., c('.yhat_lwr_', '.yhat_upr_').", call. = FALSE)
  }
  prefix_yhat_lwr <- prefix_ci[1L]
  prefix_yhat_upr <- prefix_ci[2L]

  # --- Model metadata (optional injection) -----------------------------------
  if (is.null(info)) info <- model_info(fit)
  ov_all   <- info$observed_variables
  ov_cont  <- info$ov_continuous
  ov_ord   <- info$ov_ordinal
  if (length(ov_cont) == 0L) stop("No continuous observed indicators detected; nothing to predict.", call. = FALSE)

  n_groups     <- info$n_groups %||% 1L
  group_labels <- info$group_labels %||% as.character(seq_len(n_groups))

  # --- Factor scores + observed variables (one pass) --------------------------
  # FS SEs only if continuous-only; otherwise ignore and warn
  request_fs_se <- isTRUE(se_fs) && length(ov_ord) == 0L
  if (isTRUE(se_fs) && length(ov_ord) > 0L) {
    warning("Factor-score SEs are only available for continuous-only models; disabling `se_fs` because at least one ordinal indicator is present.")
  }

  if (is.null(data)) {
    data <- lavPredict_parallel(
      fit,
      return_type  = "list",
      se           = request_fs_se,
      prefix_se_fs = prefix_se_fs
    )
  }
  fs_list <- if (is.data.frame(data)) list(data) else data
  if (length(fs_list) != n_groups) {
    stop("lavPredict returned ", length(fs_list), " group table(s), expected ", n_groups, ".", call. = FALSE)
  }

  # --- Row ids used in estimation (.rid from case.idx) ------------------------
  case_idx <- tryCatch(lavaan::lavInspect(fit, "case.idx"), error = function(e) NULL)

  # --- Measurement parameters -------------------------------------------------
  est <- lavaan::lavInspect(fit, "est")
  est_list <- if (is.list(est) && "lambda" %in% names(est)) list(est) else est

  # --- Extract eta-hat per group ---------------------------------------------
  eta_list <- lapply(fs_list, function(df) as.matrix(df[, info$latent_variables, drop = FALSE]))

  # Decide whether we need yhat internally
  yhat_needed <- isTRUE(yhat) || isTRUE(resid) || isTRUE(ci) || isTRUE(se_yhat)

  # --- Compute predictions for continuous indicators (when needed) ------------
  yhat_list <- NULL
  if (yhat_needed) {
    yhat_list <- vector("list", n_groups)
    names(yhat_list) <- group_labels

    for (g in seq_len(n_groups)) {
      lam_g <- est_list[[g]]$lambda
      nu_g  <- est_list[[g]]$nu

      # Align Lambda: columns = latent order, rows = observed continuous order
      if (!is.null(colnames(lam_g))) lam_g <- lam_g[, info$latent_variables, drop = FALSE]
      if (!is.null(rownames(lam_g))) {
        lam_g <- lam_g[ov_cont, , drop = FALSE]
      } else {
        row_idx <- match(ov_cont, ov_all)
        lam_g   <- lam_g[row_idx, , drop = FALSE]
      }

      # Align nu to continuous indicators
      if (!is.null(names(nu_g))) {
        nu_g <- nu_g[ov_cont]
      } else {
        row_idx <- match(ov_cont, ov_all)
        nu_g    <- nu_g[row_idx]
      }

      eta_g  <- eta_list[[g]]                  # N x n_latent
      yhat_g <- eta_g %*% t(lam_g)             # N x n_cont
      yhat_g <- sweep(yhat_g, 2, nu_g, `+`)    # add intercepts
      colnames(yhat_g) <- paste0(prefix_yhat, ov_cont)
      yhat_list[[g]] <- tibble::as_tibble(yhat_g)
    }
  }

  # --- Build FS table with .rid, .gid, .group --------------------------------
  fs_tbl_parts <- vector("list", n_groups)
  for (g in seq_len(n_groups)) {
    tmp <- tibble::as_tibble(fs_list[[g]])

    rid_g <- tryCatch({ if (is.list(case_idx)) case_idx[[g]] else case_idx }, error = function(e) NULL)
    if (is.null(rid_g) || length(rid_g) != nrow(tmp)) rid_g <- seq_len(nrow(tmp))

    tmp <- tibble::add_column(tmp, .rid = rid_g,          .before = 1)
    tmp <- tibble::add_column(tmp, .gid = g,              .before = 2)
    tmp <- tibble::add_column(tmp, .group = group_labels[g], .after = ".gid")
    fs_tbl_parts[[g]] <- tmp
  }
  out <- dplyr::bind_rows(fs_tbl_parts)

  # --- Place factor scores and FS SEs -----------------------------------------
  eta_names <- info$latent_variables
  fs_cols <- eta_names[eta_names %in% names(out)]
  anchor <- if (".group" %in% names(out)) ".group" else ".gid"
  if (length(fs_cols)) {
    out <- dplyr::relocate(out, dplyr::all_of(fs_cols), .after = dplyr::all_of(anchor))
  }

  # Move FS SE columns as a block right after FS (if we requested them)
  if (isTRUE(request_fs_se)) {
    se_cols <- paste0(prefix_se_fs, fs_cols)
    se_cols <- se_cols[se_cols %in% names(out)]
    if (length(fs_cols) && length(se_cols)) {
      out <- dplyr::relocate(out, dplyr::all_of(se_cols),
                             .after = dplyr::all_of(fs_cols[length(fs_cols)]))
    }
  }

  # --- Attach predictions if requested ----------------------------------------
  if (yhat_needed) {
    yhat_tbl <- dplyr::bind_rows(yhat_list)
    if (isTRUE(yhat)) {
      out <- dplyr::bind_cols(out, yhat_tbl)
    }
  }

  # --- Residuals --------------------------------------------------------------
  if (isTRUE(resid)) {
    y_obs  <- as.matrix(out[, ov_cont, drop = FALSE])
    y_pred <- if (isTRUE(yhat)) {
      as.matrix(out[, paste0(prefix_yhat, ov_cont), drop = FALSE])
    } else {
      as.matrix(dplyr::bind_rows(yhat_list)[, paste0(prefix_yhat, ov_cont), drop = FALSE])
    }
    resid_m <- y_obs - y_pred
    colnames(resid_m) <- paste0(prefix_resid, ov_cont)
    out <- dplyr::bind_cols(out, tibble::as_tibble(resid_m))
  }

  # --- Delta-method for yhat: CIs and/or SEs ----------------------------------
  if (isTRUE(ci) || isTRUE(se_yhat)) {
    V <- tryCatch({
      if (is.null(vcov_type)) lavaan::vcov(fit) else lavaan::vcov(fit, type = vcov_type)
    }, error = function(e) NULL)
    if (is.null(V)) stop("vcov(fit) is not available; cannot compute delta-method for yhat.", call. = FALSE)

    PT     <- lavaan::parTable(fit)
    latent <- info$latent_variables
    single_group <- (n_groups == 1L)
    if (isTRUE(ci)) z <- stats::qnorm((1 + level) / 2)

    get_free_lambda <- function(g, j, k) {
      if (single_group) {
        w <- which(PT$op == "=~" & PT$lhs == k & PT$rhs == j)
      } else {
        w <- which(PT$group == g & PT$op == "=~" & PT$lhs == k & PT$rhs == j)
      }
      if (!length(w)) return(0L)
      PT$free[w[1L]]
    }
    get_free_nu <- function(g, j) {
      if (single_group) {
        w <- which(PT$op == "~1" & PT$lhs == j)
      } else {
        w <- which(PT$group == g & PT$op == "~1" & PT$lhs == j)
      }
      if (!length(w)) return(0L)
      PT$free[w[1L]]
    }

    yhat_all <- if (isTRUE(yhat)) {
      as.matrix(out[, paste0(prefix_yhat, ov_cont), drop = FALSE])
    } else {
      as.matrix(dplyr::bind_rows(yhat_list)[, paste0(prefix_yhat, ov_cont), drop = FALSE])
    }

    # Container for SE_yhat
    se_yhat_df <- if (isTRUE(se_yhat)) as.data.frame(matrix(NA_real_, nrow(out), length(ov_cont),
                                                            dimnames = list(NULL, paste0(prefix_se_yhat, ov_cont))))
    else NULL

    # Pre-create CI columns if needed
    if (isTRUE(ci)) {
      for (j in ov_cont) {
        out[[paste0(prefix_yhat_lwr, j)]] <- NA_real_
        out[[paste0(prefix_yhat_upr, j)]] <- NA_real_
      }
    }

    for (g in seq_len(n_groups)) {
      rows_g <- which(out$.gid == g)
      if (!length(rows_g)) next

      eta_mat <- as.matrix(out[rows_g, latent, drop = FALSE])

      for (j in ov_cont) {
        lam_idx_raw <- vapply(latent, function(k) get_free_lambda(g, j, k), integer(1))
        nu_idx      <- get_free_nu(g, j)

        ids <- c(lam_idx_raw[lam_idx_raw > 0L], if (nu_idx > 0L) nu_idx else integer(0))

        if (!length(ids)) {
          se_row <- numeric(length(rows_g))
        } else {
          uniq <- unique(ids)

          if (any(lam_idx_raw > 0L)) {
            A_lam   <- eta_mat[, lam_idx_raw > 0L, drop = FALSE]   # N_g x (#free lambdas)
            col_ids <- lam_idx_raw[lam_idx_raw > 0L]
            A <- do.call(cbind, lapply(uniq, function(id) {
              if (id == nu_idx) {
                rep(0, nrow(eta_mat))
              } else {
                idx <- which(col_ids == id)
                if (length(idx)) rowSums(A_lam[, idx, drop = FALSE]) else rep(0, nrow(eta_mat))
              }
            }))
          } else {
            A <- matrix(0, nrow = length(rows_g), ncol = length(uniq))
          }

          if (nu_idx > 0L) {
            hit <- match(nu_idx, uniq, nomatch = 0L)
            if (hit > 0L) A[, hit] <- A[, hit] + 1
          }

          Vsub <- V[uniq, uniq, drop = FALSE]
          se2  <- rowSums((A %*% Vsub) * A)
          se_row <- sqrt(pmax(se2, 0))
        }

        if (isTRUE(se_yhat)) {
          se_yhat_df[rows_g, paste0(prefix_se_yhat, j)] <- se_row
        }
        if (isTRUE(ci)) {
          mu <- as.numeric(yhat_all[rows_g, paste0(prefix_yhat, j), drop = TRUE])
          out[rows_g, paste0(prefix_yhat_lwr, j)] <- mu - z * se_row
          out[rows_g, paste0(prefix_yhat_upr, j)] <- mu + z * se_row
        }
      }
    }

    if (isTRUE(se_yhat)) {
      out <- dplyr::bind_cols(out, tibble::as_tibble(se_yhat_df))
    }
  }

  # --- Place CI columns immediately after their .yhat_* (if yhat is visible) --
  if (isTRUE(ci) && isTRUE(yhat)) {
    for (j in ov_cont) {
      y_col   <- paste0(prefix_yhat, j)
      lwr_col <- paste0(prefix_yhat_lwr, j)
      upr_col <- paste0(prefix_yhat_upr, j)
      if (all(c(y_col, lwr_col, upr_col) %in% names(out))) {
        out <- dplyr::relocate(out, dplyr::all_of(c(lwr_col, upr_col)),
                               .after = dplyr::all_of(y_col))
      }
    }
  }

  # Cast numerics
  out <- dplyr::mutate(out, dplyr::across(dplyr::where(is.numeric), as.double))

  tibble::as_tibble(out)
}
