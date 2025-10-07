#' Add predicted indicator values, residuals, and (optionally) CIs for continuous indicators
#'
#' Computes \eqn{\hat{y} = \nu + \Lambda \eta} for observed **continuous**
#' indicators and returns a tibble with (a) original factor-score output
#' (including observed variables), (b) predicted values prefixed by
#' \code{prefix_yhat}, (c) residuals (observed - predicted) prefixed by
#' \code{prefix_resid}, and (optionally) (d) Wald-type confidence intervals
#' for predicted values via a delta method that treats factor scores as fixed.
#' If \code{SE = TRUE}, and the model is continuous-only, factor-score
#' standard errors are attached by delegating to \code{lavPredict_parallel(se = TRUE)}.
#'
#' @param fit A fitted \code{lavaan} object with \code{meanstructure = TRUE}.
#' @param prefix_yhat Character scalar, prefix for predicted columns. Default \code{".yhat_"}.
#' @param prefix_resid Character scalar, prefix for residual columns. Default \code{".resid_"}.
#' @param ci Logical; if \code{TRUE}, adds Wald CIs for \eqn{\hat{y}} via delta method. Default \code{FALSE}.
#' @param level Confidence level, default \code{0.95}.
#' @param prefix_ci Length-2 character vector with lower/upper CI prefixes. Default \code{c(".yhat_lwr_", ".yhat_upr_")}.
#' @param vcov_type Optional character passed to \code{lavaan::vcov(fit, type = vcov_type)}.
#'   If \code{NULL}, uses the default \code{vcov(fit)}.
#' @param SE Logical; if \code{TRUE}, request factor-score SEs (continuous models only).
#'   Passed through to \code{lavPredict_parallel(se = SE)}. Default \code{FALSE}.
#' @param prefix_se Character; prefix for factor-score SE columns (e.g., \code{".se_"}).
#'   Passed to \code{lavPredict_parallel(se_prefix = prefix_se)}. Default \code{".se_"}.
#'
#' @return A \code{tibble}. For multi-group models, includes a \code{group}
#'   column with group labels. Predicted/residual columns are added **only**
#'   for continuous indicators; ordinal indicators are skipped. If \code{ci = TRUE},
#'   lower/upper CI columns are added for each predicted indicator. If \code{SE = TRUE}
#'   (continuous models), columns named \code{paste0(prefix_se, <factor>)} are present
#'   in the factor-score block.
#'
#' @export
augment_continuous_old <- function(fit,
                               prefix_yhat  = ".yhat_",
                               prefix_resid = ".resid_",
                               ci           = FALSE,
                               level        = 0.95,
                               prefix_ci    = c(".yhat_lwr_", ".yhat_upr_"),
                               vcov_type    = NULL,
                               se           = FALSE,
                               prefix_se    = ".se_") {

  # -- Defensive check ---------------------------------------------------------
  # requires: converged fit, meanstructure, latent vars
  .assert_lavaan_fit(
    fit,
    require_converged     = TRUE,
    require_meanstructure = TRUE,
    require_latent        = TRUE
  )

  # -- Validate prefix_ci ------------------------------------------------------
  if (!is.character(prefix_ci) || length(prefix_ci) != 2L || any(!nzchar(prefix_ci))) {
    stop("`prefix_ci` must be two non-empty strings, e.g., c('.yhat_lwr_', '.yhat_upr_').", call. = FALSE)
  }
  prefix_yhat_lwr <- prefix_ci[1L]
  prefix_yhat_upr <- prefix_ci[2L]

  # --- Model metadata ---------------------------------------------------------
  info    <- model_info(fit)
  ov_all  <- info$observed_variables
  ov_ord  <- lavaan::lavNames(fit, type = "ov.ord")
  ov_cont <- setdiff(ov_all, ov_ord)
  if (length(ov_cont) == 0L) stop("No continuous observed indicators detected; nothing to predict.", call. = FALSE)
  if (length(ov_ord) > 0L) {
    warning("Skipping ordinal indicators: ", paste(ov_ord, collapse = ", "),
            ". Predicted values and residuals are computed only for continuous indicators.")
  }

  n_groups     <- info$n_groups
  group_labels <- if (!is.null(info$group_labels)) info$group_labels else as.character(seq_len(n_groups))

  # --- Factor scores + observed variables (delegate SE handling) --------------
  # Note: lavPredict_parallel(se=TRUE) adds .se_<factor> columns for continuous models.
  fs_and_ov <- lavPredict_parallel(
    fit,
    return_type = "list",
    se          = isTRUE(se),
    prefix_se   = prefix_se
  )

  fs_list <- if (is.data.frame(fs_and_ov)) list(fs_and_ov) else fs_and_ov
  if (length(fs_list) != n_groups) {
    stop("lavPredict returned ", length(fs_list), " group table(s), expected ", n_groups, ".", call. = FALSE)
  }

  # --- Measurement parameters from lavInspect(..., "est") ---------------------
  est <- lavaan::lavInspect(fit, "est")
  est_list <- if (is.list(est) && "lambda" %in% names(est)) list(est) else est

  # --- Extract eta-hat per group ---------------------------------------------
  eta_list <- lapply(fs_list, function(df) as.matrix(df[, info$latent_variables, drop = FALSE]))

  # --- Compute predictions for continuous indicators --------------------------
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

  # --- Bind and compute residuals (explicit numeric .gid 1..G) ----------------
  fs_tbl_parts <- vector("list", n_groups)
  for (g in seq_len(n_groups)) {
    tmp <- tibble::as_tibble(fs_list[[g]])
    tmp <- tibble::add_column(tmp, .gid  = g,                 .before = 1)
    tmp <- tibble::add_column(tmp, group = group_labels[g],   .after  = 1)  # <-- right after .gid
    fs_tbl_parts[[g]] <- tmp
  }
  fs_tbl <- dplyr::bind_rows(fs_tbl_parts)
  fs_tbl$group <- group_labels[fs_tbl$.gid]

  yhat_tbl <- dplyr::bind_rows(yhat_list)
  out <- dplyr::bind_cols(fs_tbl, yhat_tbl)

  # Vectorized residuals
  y_obs   <- as.matrix(out[, ov_cont, drop = FALSE])
  y_pred  <- as.matrix(out[, paste0(prefix_yhat, ov_cont), drop = FALSE])
  resid_m <- y_obs - y_pred
  colnames(resid_m) <- paste0(prefix_resid, ov_cont)
  out <- dplyr::bind_cols(out, tibble::as_tibble(resid_m))

  # --- Optional: delta-method CIs for yhat -----------------------------------
  if (isTRUE(ci)) {
    z <- stats::qnorm((1 + level) / 2)

    V <- tryCatch({
      if (is.null(vcov_type)) lavaan::vcov(fit) else lavaan::vcov(fit, type = vcov_type)
    }, error = function(e) NULL)
    if (is.null(V)) stop("vcov(fit) is not available; cannot compute delta-method CIs.", call. = FALSE)

    PT     <- lavaan::parTable(fit)
    latent <- info$latent_variables
    single_group <- (n_groups == 1L)

    # Helpers to fetch 'free' indices; tolerate PT$group being NA for single-group fits
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

    # Pre-create CI columns
    for (j in ov_cont) {
      out[[paste0(prefix_yhat_lwr, j)]] <- NA_real_
      out[[paste0(prefix_yhat_upr, j)]] <- NA_real_
    }

    # Loop by explicit numeric .gid (robust across labelings)
    for (g in seq_len(n_groups)) {
      rows_g <- which(out$.gid == g)
      if (!length(rows_g)) next

      eta_mat <- as.matrix(out[rows_g, latent, drop = FALSE])

      for (j in ov_cont) {
        lam_idx_raw <- vapply(latent, function(k) get_free_lambda(g, j, k), integer(1))
        nu_idx      <- get_free_nu(g, j)

        ids <- c(lam_idx_raw[lam_idx_raw > 0L], if (nu_idx > 0L) nu_idx else integer(0))
        if (!length(ids)) {
          se <- numeric(length(rows_g))
        } else {
          uniq <- unique(ids)

          if (any(lam_idx_raw > 0L)) {
            A_lam   <- eta_mat[, lam_idx_raw > 0L, drop = FALSE]
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
          se2  <- rowSums((A %*% Vsub) * A)  # g' V g per row
          se   <- sqrt(pmax(se2, 0))
        }

        mu <- as.numeric(out[rows_g, paste0(prefix_yhat, j), drop = TRUE])
        out[rows_g, paste0(prefix_yhat_lwr, j)] <- mu - z * se
        out[rows_g, paste0(prefix_yhat_upr, j)] <- mu + z * se
      }
    }
  }

  # Cast any <lvn.mtrx> etc. to plain doubles
  out <- dplyr::mutate(
    out,
    dplyr::across(
      dplyr::where(is.numeric),
      as.double
    )
  )

  tibble::as_tibble(out)
}
