#' Augment ordinal SEM output with factor scores, original data, and latent predictions (y*)
#'
#' Works for single- and multi-group lavaan models with **only ordinal indicators**.
#' For each ordinal indicator j, computes latent continuous predictor
#'   y*_j = nu_j + lambda_j• %*% eta
#' and returns it as `.ystar_<item>`.
#'
#' @param fit A fitted lavaan object. Intended for ordinal (ordered-categorical) models.
#' @param prefix_ystar Prefix for latent predicted values (default ".ystar_").
#'
#' @return A tibble containing:
#'   - factor scores for all latent variables,
#'   - original observed indicators (as returned by lavPredict(..., append.data=TRUE)),
#'   - `.ystar_<item>` columns for each ordinal indicator.
#'
#' @export
augment3 <- function(fit, prefix_ystar = ".ystar_") {
  # -- Basic validation: rely on project helper if available -------------------
  .assert_lavaan_fit(
    fit,
    require_converged     = TRUE,
    require_latent        = TRUE,
    # meanstructure is not strictly required for ordinal; intercepts are in thresholds,
    # but we still allow it either way
    require_meanstructure = NA
  )

  # -- Model info --------------------------------------------------------------
  info     <- model_info(fit)
  ov_all   <- info$observed_variables
  ov_ord   <- lavaan::lavNames(fit, type = "ov.ord")
  ov_cont  <- setdiff(ov_all, ov_ord)

  if (length(ov_ord) == 0L) {
    stop("augment3() expects an ordinal model: no ordinal indicators detected.", call. = FALSE)
  }
  if (length(ov_cont) > 0L) {
    stop("augment3() currently supports ordinal-only models. Detected continuous indicators: ",
         paste(ov_cont, collapse = ", "), call. = FALSE)
  }

  n_groups     <- info$n_groups
  group_labels <- if (!is.null(info$group_labels)) info$group_labels else as.character(seq_len(n_groups))

  # -- Factor scores + observed data per group --------------------------------
  fs_and_ov <- lavaan::lavPredict(
    fit,
    transform = FALSE,
    append.data = TRUE,
    assemble = FALSE,
    drop.list.single.group = FALSE
  )
  fs_list <- if (is.data.frame(fs_and_ov)) list(fs_and_ov) else fs_and_ov
  if (length(fs_list) != n_groups) {
    stop("lavPredict returned ", length(fs_list), " group table(s), expected ", n_groups, ".", call. = FALSE)
  }

  # -- Measurement parameters from lavInspect(..., "est") ----------------------
  est <- lavaan::lavInspect(fit, "est")
  est_list <- if (is.list(est) && "lambda" %in% names(est)) list(est) else est

  # -- Helper to align Lambda/nu to the ordinal indicators --------------------
  .align_Lambda_nu <- function(lam_g, nu_g, ov_subset, ov_all, latent_order) {
    # columns = latent order; rows = observed in ov_subset order
    if (!is.null(colnames(lam_g))) lam_g <- lam_g[, latent_order, drop = FALSE]
    if (!is.null(rownames(lam_g))) {
      lam_g <- lam_g[ov_subset, , drop = FALSE]
    } else {
      row_idx <- match(ov_subset, ov_all)
      lam_g   <- lam_g[row_idx, , drop = FALSE]
    }
    # nu may be absent/zero-length in ordinal probit; treat missing as 0
    if (is.null(nu_g) || length(nu_g) == 0L) {
      nu_g <- rep(0, length(ov_subset))
      names(nu_g) <- ov_subset
    } else if (!is.null(names(nu_g))) {
      nu_g <- nu_g[ov_subset]
    } else {
      row_idx <- match(ov_subset, ov_all)
      nu_g    <- nu_g[row_idx]
    }
    list(lambda = lam_g, nu = nu_g)
  }

  # -- Compute y* for ordinal indicators --------------------------------------
  ystar_list <- vector("list", n_groups)
  names(ystar_list) <- group_labels

  # factor scores matrices per group
  eta_list <- lapply(fs_list, function(df) as.matrix(df[, info$latent_variables, drop = FALSE]))

  for (g in seq_len(n_groups)) {
    lam_g <- est_list[[g]]$lambda
    nu_g  <- est_list[[g]]$nu %||% numeric(0)

    al <- .align_Lambda_nu(lam_g, nu_g, ov_subset = ov_ord,
                           ov_all = ov_all, latent_order = info$latent_variables)

    eta_g   <- eta_list[[g]]                          # N x n_latent
    ystar_g <- sweep(eta_g %*% t(al$lambda), 2, al$nu, `+`)  # N x n_ord
    colnames(ystar_g) <- paste0(prefix_ystar, ov_ord)
    ystar_list[[g]] <- tibble::as_tibble(ystar_g)
  }

  # -- Bind outputs with explicit numeric .gid to avoid coercion issues --------
  fs_tbl_parts <- vector("list", n_groups)
  for (g in seq_len(n_groups)) {
    fs_tbl_parts[[g]] <- tibble::add_column(
      tibble::as_tibble(fs_list[[g]]),
      .gid = g,
      .before = 1
    )
  }
  fs_tbl <- dplyr::bind_rows(fs_tbl_parts)
  fs_tbl$group <- group_labels[fs_tbl$.gid]

  ystar_tbl <- dplyr::bind_rows(ystar_list)

  out <- dplyr::bind_cols(fs_tbl, ystar_tbl)

  # -- Cast custom numeric-like classes (e.g., <lvn.mtrx>) to plain doubles ----
  out <- dplyr::mutate(
    out,
    dplyr::across(
      dplyr::where(is.numeric),
      as.double
    )
  )

  tibble::as_tibble(out)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
