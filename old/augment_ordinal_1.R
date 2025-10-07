#' Augment ordinal SEM output with factor scores, original data, latent predictions (y*),
#' and per-category probabilities (robust to missing appended columns)
#'
#' Works for single- and multi-group lavaan models with **only ordinal indicators**.
#' For each ordinal indicator j, computes:
#'   - latent predictor:         y*_j = nu_j + lambda_j• %*% eta   (scaled to unit-variance if theta)
#'   - category probabilities:   .pr_<level>__<item>  (per observed or threshold-inferred level)
#'
#' @param fit A fitted lavaan object intended for ordinal (ordered-categorical) models.
#' @param prefix_ystar Prefix for latent predicted values (default ".ystar_").
#'
#' @return A tibble containing factor scores, original indicators (if present in appended data),
#'         `.ystar_<item>`, and `.pr_<level>__<item>` columns plus `.gid` (numeric group id)
#'         and `group` (group label).
#' @export
augment_ordinal <- function(fit, prefix_ystar = ".ystar_") {
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
  info    <- model_info(fit)  # your internal helper (guaranteed to work)

  # Observed variables (all / ordinal / continuous)
  ov_all  <- info$observed_variables
  ov_ord  <- info$ov_ordinal
  ov_cont <- info$ov_continuous

  # Basic validations for this function's scope
  if (length(ov_ord) == 0L) {
    stop("augment4() expects an ordinal model: no ordinal indicators detected.", call. = FALSE)
  }
  if (length(ov_cont) > 0L) {
    stop("augment4() currently supports ordinal-only models. Detected continuous indicators: ",
         paste(ov_cont, collapse = ", "), call. = FALSE)
  }

  # Grouping + parameterization (delta/theta)
  n_groups     <- info$n_groups %||% 1L
  group_labels <- info$group_labels %||% as.character(seq_len(n_groups))
  param        <- info$parameterization %||% "delta"

  # -- Factor scores + observed data per group --------------------------------
  fs_and_ov <- lavPredict_parallel(
    fit,
    return_type = "list"
  )

  fs_list <- if (is.data.frame(fs_and_ov)) list(fs_and_ov) else fs_and_ov
  if (length(fs_list) != n_groups) {
    stop("lavPredict returned ", length(fs_list), " group table(s), expected ", n_groups, ".", call. = FALSE)
  }

  # Factor scores per group (matrix N_g x n_latent)
  eta_list <- lapply(fs_list, function(df) as.matrix(df[, info$latent_variables, drop = FALSE]))

  # -- Measurement parameters --------------------------------------------------
  est <- lavaan::lavInspect(fit, "est")
  # est can be either a single list (single-group) or a list-of-lists (multi-group)
  est_list <- if (is.list(est) && "lambda" %in% names(est)) list(est) else est
  if (length(est_list) != n_groups) {
    stop("Expected ", n_groups, " 'est' blocks; got ", length(est_list), ".", call. = FALSE)
  }

  # Thresholds (vector or list across groups)
  th_raw <- tryCatch(lavaan::lavInspect(fit, "th"), error = function(e) NULL)
  th_list <- if (is.null(th_raw)) {
    replicate(n_groups, setNames(numeric(0), character(0)), simplify = FALSE)
  } else if (is.list(th_raw)) {
    th_raw
  } else {
    list(th_raw)
  }
  if (length(th_list) != n_groups) {
    # recycle single vector across groups if needed
    if (length(th_list) == 1L && n_groups > 1L) {
      th_list <- rep(th_list, n_groups)
    } else {
      stop("Unexpected thresholds structure: got ", length(th_list), " elements for ", n_groups, " groups.", call. = FALSE)
    }
  }
  names(th_list) <- group_labels

  # -- Helpers -----------------------------------------------------------------
  # Align Lambda/nu to given indicators
  .align_Lambda_nu <- function(lam_g, nu_g, ov_subset, ov_all, latent_order) {
    # columns = latent variables: enforce order
    if (!is.null(colnames(lam_g))) lam_g <- lam_g[, latent_order, drop = FALSE]
    # rows = observed variables
    if (!is.null(rownames(lam_g))) {
      lam_g <- lam_g[ov_subset, , drop = FALSE]
    } else {
      row_idx <- match(ov_subset, ov_all); lam_g <- lam_g[row_idx, , drop = FALSE]
    }
    # intercepts/threshold centers (nu)
    if (is.null(nu_g) || length(nu_g) == 0L) {
      nu_g <- rep(0, length(ov_subset)); names(nu_g) <- ov_subset
    } else if (!is.null(names(nu_g))) {
      nu_g <- nu_g[ov_subset]
    } else {
      row_idx <- match(ov_subset, ov_all); nu_g <- nu_g[row_idx]
    }
    list(lambda = lam_g, nu = nu_g)
  }

  # Extract thresholds for variable j, ordered by numeric suffix t1 < t2 < ...
  .get_tau_for <- function(th_vec, j) {
    nm <- names(th_vec); if (is.null(nm)) return(numeric(0))
    sel <- startsWith(nm, paste0(j, "|"))
    tau <- th_vec[sel]
    if (!length(tau)) return(numeric(0))
    # order by numeric index after '|t'
    idx <- suppressWarnings(as.integer(gsub("^.*\\|t(\\d+)$", "\\1", names(tau))))
    idx[is.na(idx)] <- seq_along(tau)
    tau[order(idx)]
  }

  # Try to get levels from appended data; otherwise infer 1..m from thresholds
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

  # Predicate for gentle numeric coercion (do not convert integers)
  .needs_double_from_special <- function(x) {
    inherits(x, "lvn.mtrx") || inherits(x, "lavaan.matrix")
  }

  # -- Compute y* and probabilities -------------------------------------------
  ystar_list <- vector("list", n_groups)
  probs_list <- vector("list", n_groups)
  names(ystar_list) <- names(probs_list) <- group_labels

  for (g in seq_len(n_groups)) {
    est_g <- est_list[[g]]
    lam_g <- est_g$lambda
    nu_g  <- est_g$nu %||% numeric(0)

    al <- .align_Lambda_nu(
      lam_g, nu_g,
      ov_subset   = ov_ord,
      ov_all      = ov_all,
      latent_order = info$latent_variables
    )

    eta_g    <- eta_list[[g]]                                # N x n_latent
    ystar_g  <- sweep(eta_g %*% t(al$lambda), 2, al$nu, `+`) # N x n_ord
    colnames(ystar_g) <- ov_ord

    # --- NOTE: for theta parametrization, DO NOT scale y* here.
    # We will scale the (tau - y*) difference by residual SD item-by-item below.
    sig_eps <- NULL
    if (identical(param, "theta")) {
      theta_g <- est_g$theta
      if (!is.null(theta_g)) {
        # diag(theta) -> residual variances; align to ov_ord
        if (!is.null(rownames(theta_g))) {
          sig_eps <- sqrt(diag(theta_g)[ov_ord])
        } else {
          idx <- match(ov_ord, ov_all)
          sig_eps <- sqrt(diag(theta_g)[idx])
        }
        sig_eps[!is.finite(sig_eps) | sig_eps <= 0] <- 1
      } else {
        sig_eps <- rep(1, length(ov_ord))
      }
    }

    N_g <- nrow(ystar_g)
    probs_g <- NULL

    for (jj in seq_along(ov_ord)) {
      j <- ov_ord[jj]
      ystar_j <- ystar_g[, jj]

      # levels from data if present; otherwise infer from thresholds
      lev <- .levels_for_item(fs_list[[g]], th_list[[g]], j)
      tau_j <- as.numeric(.get_tau_for(th_list[[g]], j))

      # If levels present but count doesn't match thresholds, fall back to thresholds
      if (length(tau_j) && length(lev) != (length(tau_j) + 1L)) {
        warning("Levels in appended data for '", j, "' (", paste(lev, collapse = ", "),
                ") do not match thresholds (", length(tau_j) + 1L, " categories). ",
                "Using threshold-implied levels 1..m instead (group '", group_labels[g], "').")
        lev <- as.character(seq_len(length(tau_j) + 1L))
      }

      # If still no levels or thresholds, skip probabilities for this item
      if (!length(lev)) {
        warning("Levels for '", j, "' not found; skipping probability columns for this item (group '",
                group_labels[g], "').")
        next
      }

      # If no thresholds available, we cannot compute probabilities
      if (!length(tau_j)) {
        warning("Thresholds for '", j, "' are not available; skipping probability columns (group '",
                group_labels[g], "').")
        next
      }

      m <- length(lev)
      if (length(tau_j) != (m - 1L)) {
        warning("Threshold count for '", j, "' (", length(tau_j),
                ") does not match number of categories (", m,
                "); skipping probability columns for this item (group '",
                group_labels[g], "').")
        next
      }

      # padded thresholds
      tau_pad <- c(-Inf, tau_j, Inf)

      # --- Correct computation:
      # For 'theta' use z = (tau - y*) / sigma_eps_j; otherwise z = (tau - y*)
      if (identical(param, "theta")) {
        se_j  <- if (is.null(sig_eps)) 1 else sig_eps[jj]
        zmat  <- outer(ystar_j, tau_pad, function(y, t) (t - y) / se_j)
      } else {
        zmat  <- outer(ystar_j, tau_pad, function(y, t) (t - y))
      }

      diffs   <- stats::pnorm(zmat)
      probs_j <- diffs[, -1, drop = FALSE] - diffs[, -ncol(diffs), drop = FALSE]

      # Name columns: .pr_<level>__<item>
      colnames(probs_j) <- paste0(".pr_", make.names(lev, allow_ = TRUE), "__", j)

      probs_g <- if (is.null(probs_g)) {
        as.data.frame(probs_j, check.names = FALSE)
      } else {
        cbind(probs_g, as.data.frame(probs_j, check.names = FALSE))
      }
    }

    # y* to tibble with prefix
    ystar_df <- tibble::as_tibble(ystar_g)
    names(ystar_df) <- paste0(prefix_ystar, names(ystar_df))

    ystar_list[[g]] <- ystar_df
    probs_list[[g]] <- if (is.null(probs_g)) tibble::tibble() else tibble::as_tibble(probs_g)
  }

  # -- Bind outputs with explicit numeric .gid and group label -----------------
  fs_tbl_parts <- vector("list", n_groups)
  for (g in seq_len(n_groups)) {
    # put .gid and group label first; retain original fs+ov columns after
    fs_g  <- tibble::as_tibble(fs_list[[g]])
    fs_g  <- tibble::add_column(fs_g, .gid = g, .before = 1)
    fs_g  <- tibble::add_column(fs_g, group = group_labels[g], .after = ".gid")
    fs_tbl_parts[[g]] <- fs_g
  }
  fs_tbl <- dplyr::bind_rows(fs_tbl_parts)

  ystar_tbl <- dplyr::bind_rows(ystar_list)
  probs_tbl <- dplyr::bind_rows(probs_list)

  out <- dplyr::bind_cols(fs_tbl, ystar_tbl, probs_tbl)

  # -- Ensure plain doubles only for special classes (do not touch integers) ---
  out <- dplyr::mutate(
    out,
    dplyr::across(
      .cols = where(.needs_double_from_special),
      .fns  = ~ as.double(.x)
    )
  )

  tibble::as_tibble(out)
}

