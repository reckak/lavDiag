#' Add predicted indicator values, residuals, and (optionally) CIs for continuous indicators
#'
#' Computes \eqn{\hat{y} = \nu + \Lambda \eta} for observed **continuous**
#' indicators and returns a tibble with (a) original factor-score output
#' (including observed variables), (b) predicted values prefixed by
#' \code{prefix_yhat}, (c) residuals (observed - predicted) prefixed by
#' \code{prefix_resid}, and (optionally) (d) Wald-type confidence intervals
#' for predicted values via a delta method that treats factor scores as fixed.
#'
#' @param fit A fitted \code{lavaan} object with \code{meanstructure = TRUE}.
#' @param prefix_yhat Character scalar, prefix for predicted columns. Default \code{".yhat_"}.
#' @param prefix_resid Character scalar, prefix for residual columns. Default \code{".resid_"}.
#' @param ci One of \code{c("none","delta")}. If \code{"delta"}, adds Wald CIs
#'   using \code{vcov(fit)} and a delta method w.r.t. \eqn{\nu,\Lambda}. Default \code{"none"}.
#' @param level Confidence level, default \code{0.95}.
#' @param prefix_ci Character vector of length 2 with prefixes for lower/upper CI columns. Default: \code{c(".yhat_lwr_", ".yhat_upr_")}.
#' @return A \code{tibble}. For multi-group models, includes a \code{group}
#'   column with group labels. Predicted/residual columns are added **only**
#'   for continuous indicators; ordinal indicators are skipped. If \code{ci != "none"},
#'   lower/upper CI columns are added for each predicted indicator.
#'
#' @details
#' CIs are computed as \eqn{\hat y_{ij} \pm z_{\alpha/2} \cdot \mathrm{SE}(\hat y_{ij})},
#' where the standard error uses a delta method with gradient
#' \eqn{\partial \hat y_{ij} / \partial \nu_j = 1} and
#' \eqn{\partial \hat y_{ij} / \partial \lambda_{jk} = \hat\eta_{ik}} for loadings
#' on the latent variables of indicator \eqn{j}. The covariance of free parameters
#' is taken from \code{vcov(fit)}; fixed parameters contribute zero.
#'
#' @examples
#' # fit <- lavaan::cfa('F =~ y1 + y2 + y3', data = dat, meanstructure = TRUE)
#' # out <- augment(fit, ci = "delta", level = 0.95)
#' # dplyr::select(out, dplyr::starts_with(".yhat_"))
#'
#' @export
augment2 <- function(fit,
                     prefix_yhat  = ".yhat_",
                     prefix_resid = ".resid_",
                     ci           = c("delta","none"),
                     level        = 0.95,
                     prefix_ci    = c(".yhat_lwr_", ".yhat_upr_")) {

  ci <- match.arg(ci)

  # -- Validate lavaan fit (use the existing helper from functions.txt) -------
  .assert_lavaan_fit(
    fit,
    require_converged     = TRUE,
    require_meanstructure = TRUE,
    require_latent        = TRUE
  )

  # -- Validate prefix_ci ------------------------------------------------------
  if (!is.character(prefix_ci) || length(prefix_ci) != 2L) {
    stop("`prefix_ci` must be a character vector of length 2, e.g., c('.yhat_lwr_', '.yhat_upr_').",
         call. = FALSE)
  }
  if (any(!nzchar(prefix_ci))) {
    stop("Both elements of `prefix_ci` must be non-empty strings.", call. = FALSE)
  }
  prefix_yhat_lwr <- prefix_ci[1L]
  prefix_yhat_upr <- prefix_ci[2L]

  # --- Model metadata ---------------------------------------------------------
  info    <- model_info(fit)
  ov_all  <- info$observed_variables
  ov_ord  <- lavaan::lavNames(fit, type = "ov.ord")
  ov_cont <- setdiff(ov_all, ov_ord)
  if (length(ov_cont) == 0L) {
    stop("No continuous observed indicators detected; nothing to predict.", call. = FALSE)
  }
  if (length(ov_ord) > 0L) {
    warning("Skipping ordinal indicators: ", paste(ov_ord, collapse = ", "),
            ". Predicted values and residuals are computed only for continuous indicators.")
  }

  n_groups     <- info$n_groups
  group_labels <- if (!is.null(info$group_labels)) info$group_labels else as.character(seq_len(n_groups))

  # --- Factor scores + observed variables per group ---------------------------
  fs_and_ov <- lavaan::lavPredict(
    fit, transform = FALSE, append.data = TRUE,
    assemble = FALSE, drop.list.single.group = FALSE
  )
  fs_list <- if (is.data.frame(fs_and_ov)) list(fs_and_ov) else fs_and_ov
  if (length(fs_list) != n_groups) {
    stop("lavPredict returned ", length(fs_list), " group table(s), expected ", n_groups, ".", call. = FALSE)
  }

  # --- Measurement parameters ONLY from lavInspect(..., "est") ----------------
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

  # --- Bind and compute residuals --------------------------------------------
  fs_tbl <- lapply(fs_list, tibble::as_tibble) |> dplyr::bind_rows(.id = "group")
  idx <- suppressWarnings(as.integer(fs_tbl$group))
  if (all(!is.na(idx)) && length(group_labels) >= max(idx)) fs_tbl$group <- group_labels[idx]

  yhat_tbl <- dplyr::bind_rows(yhat_list)
  out <- dplyr::bind_cols(fs_tbl, yhat_tbl)

  for (nm in ov_cont) {
    obs   <- nm
    yhat  <- paste0(prefix_yhat, nm)
    resid <- paste0(prefix_resid, nm)
    if (!obs %in% names(out))  stop("Observed variable '", obs,  "' not found in lavPredict output.", call. = FALSE)
    if (!yhat %in% names(out)) stop("Predicted column '", yhat, "' is missing; please report a bug.", call. = FALSE)
    out[[resid]] <- as.numeric(out[[obs]]) - as.numeric(out[[yhat]])
  }

  # --- Optional: delta-method CIs for yhat -----------------------------------
  if (ci == "delta") {
    # z-quantile
    z <- stats::qnorm( (1 + level) / 2 )

    # Free-parameter covariance
    V <- tryCatch(lavaan::vcov(fit), error = function(e) NULL)
    if (is.null(V)) stop("vcov(fit) is not available; cannot compute delta-method CIs.", call. = FALSE)

    PT <- lavaan::parTable(fit)

    # Helpers to fetch 'free' indices
    get_free_lambda <- function(g, j, k) {
      w <- which(PT$group == g & PT$op == "=~" & PT$lhs == k & PT$rhs == j)
      if (length(w) == 0L) return(0L)
      PT$free[w[1L]]
    }
    get_free_nu <- function(g, j) {
      w <- which(PT$group == g & PT$op == "~1" & PT$lhs == j)
      if (length(w) == 0L) return(0L)
      PT$free[w[1L]]
    }

    latent <- info$latent_variables
    map_list <- vector("list", n_groups)
    names(map_list) <- group_labels

    for (g in seq_len(n_groups)) {
      idx_j <- lapply(ov_cont, function(j) {
        lam_idx <- vapply(latent, function(k) get_free_lambda(g, j, k), integer(1))
        nu_idx  <- get_free_nu(g, j)
        list(lam_idx = lam_idx, nu_idx = nu_idx)
      })
      names(idx_j) <- ov_cont
      map_list[[g]] <- idx_j
    }

    # Pre-create lwr/upr columns
    for (j in ov_cont) {
      lwr_nm <- paste0(prefix_yhat_lwr, j)
      upr_nm <- paste0(prefix_yhat_upr, j)
      out[[lwr_nm]] <- NA_real_
      out[[upr_nm]] <- NA_real_
    }

    grp_col <- if ("group" %in% names(out)) out$group else rep(group_labels[1L], nrow(out))
    for (g_lab in group_labels) {
      g <- match(g_lab, group_labels)
      rows_g <- which(grp_col == g_lab)
      if (!length(rows_g)) next

      eta_mat <- as.matrix(out[rows_g, latent, drop = FALSE])

      for (j in ov_cont) {
        lam_idx <- map_list[[g]][[j]]$lam_idx
        nu_idx  <- map_list[[g]][[j]]$nu_idx
        idx_use <- c(lam_idx[lam_idx > 0L], if (nu_idx > 0L) nu_idx else integer(0))
        if (length(idx_use) == 0L) {
          se <- rep(0, length(rows_g))
        } else {
          free_k <- which(lam_idx > 0L)
          A <- if (length(free_k)) as.matrix(eta_mat[, free_k, drop = FALSE]) else matrix(nrow = length(rows_g), ncol = 0)
          if (nu_idx > 0L) A <- cbind(A, 1.0)
          Vsub <- V[idx_use, idx_use, drop = FALSE]
          se2 <- rowSums((A %*% Vsub) * A)
          se  <- sqrt(pmax(se2, 0))
        }

        yhat_nm <- paste0(prefix_yhat, j)
        lwr_nm  <- paste0(prefix_yhat_lwr, j)
        upr_nm  <- paste0(prefix_yhat_upr, j)

        mu  <- as.numeric(out[rows_g, yhat_nm, drop = TRUE])
        out[rows_g, lwr_nm] <- mu - z * se
        out[rows_g, upr_nm] <- mu + z * se
      }
    }
  }

  # Cast any <lvn.mtrx> etc. to plain doubles
  out <- dplyr::mutate(out, dplyr::across(where(is.numeric), as.double))

  tibble::as_tibble(out)
}
