#' Augment ordinal SEM with factor scores, original data, and selected predictions
#'
#' Works for single- and multi-group lavaan models with **only ordinal indicators**.
#' Select what to compute via `what`:
#'   - "star"  : latent predictor y*_j = nu_j + lambda_j• %*% eta  (scaled to unit-variance if theta)
#'   - "prob"  : per-category probabilities  .pr_<level>__<item>
#'   - "exp"   : expected score on the ordinal scale  .yexp_<item>  (coded 1..m)
#'   - "class" : most likely category  .class_<item>
#'
#' @param fit  Fitted lavaan object intended for ordinal (ordered-categorical) models.
#' @param what Character vector; subset of c("star","prob","exp","class"). Default "exp".
#' @param prefix_ystar Prefix for latent predictions y* (default ".ystar_").
#'
#' @return A tibble with factor scores, original indicators, and requested predictions.
#' @export
augment3 <- function(fit,
                     what = "exp",
                     prefix_ystar = ".ystar_") {
  # --- Validate inputs ---------------------------------------------------------
  choices <- c("star","prob","exp","class")
  if (!length(what)) what <- "exp"
  what <- match.arg(what, choices = choices, several.ok = TRUE)

  .assert_lavaan_fit(
    fit,
    require_converged     = TRUE,
    require_latent        = TRUE,
    require_meanstructure = NA
  )

  # --- Model info --------------------------------------------------------------
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
  param        <- lavaan::lavInspect(fit, "options")$parameterization %||% "delta"

  # --- Factor scores + original data ------------------------------------------
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

  # --- Measurement parameters --------------------------------------------------
  est <- lavaan::lavInspect(fit, "est")
  est_list <- if (is.list(est) && "lambda" %in% names(est)) list(est) else est

  # Thresholds per group (vector or list)
  th_raw <- tryCatch(lavaan::lavInspect(fit, "th"), error = function(e) NULL)
  th_list <- if (is.null(th_raw)) {
    replicate(n_groups, setNames(numeric(0), character(0)), simplify = FALSE)
  } else if (is.list(th_raw)) {
    th_raw
  } else {
    list(th_raw)
  }

  # --- Helpers ----------------------------------------------------------------
  # Align Lambda/nu to a given indicator subset
  .align_Lambda_nu <- function(lam_g, nu_g, ov_subset, ov_all, latent_order) {
    if (!is.null(colnames(lam_g))) lam_g <- lam_g[, latent_order, drop = FALSE]
    if (!is.null(rownames(lam_g))) {
      lam_g <- lam_g[ov_subset, , drop = FALSE]
    } else {
      row_idx <- match(ov_subset, ov_all); lam_g <- lam_g[row_idx, , drop = FALSE]
    }
    # Ordinal intercepts may be absent; treat as zeros
    if (is.null(nu_g) || length(nu_g) == 0L) {
      nu_g <- rep(0, length(ov_subset)); names(nu_g) <- ov_subset
    } else if (!is.null(names(nu_g))) {
      nu_g <- nu_g[ov_subset]
    } else {
      row_idx <- match(ov_subset, ov_all); nu_g <- nu_g[row_idx]
    }
    list(lambda = lam_g, nu = nu_g)
  }

  # Extract thresholds for variable j ordered by numeric suffix t1 < t2 < ...
  .get_tau_for <- function(th_vec, j) {
    nm <- names(th_vec); if (is.null(nm)) return(numeric(0))
    sel <- startsWith(nm, paste0(j, "|"))
    tau <- th_vec[sel]
    if (!length(tau)) return(numeric(0))
    idx <- suppressWarnings(as.integer(gsub("^.*\\|t(\\d+)$", "\\1", names(tau))))
    idx[is.na(idx)] <- seq_along(tau)
    as.numeric(tau[order(idx)])
  }

  # Levels from appended data if present; else infer 1..m from thresholds
  .levels_for_item <- function(df, th_vec, j) {
    if (!is.null(df) && j %in% names(df)) {
      y <- df[[j]]
      if (is.factor(y)) return(levels(y))
      lev <- sort(unique(y))
      return(as.character(lev))
    } else {
      tau_j <- .get_tau_for(th_vec, j)
      if (!length(tau_j)) return(character(0))
      m <- length(tau_j) + 1L
      return(as.character(seq_len(m)))
    }
  }

  # Factor scores matrices per group
  eta_list <- lapply(fs_list, function(df) as.matrix(df[, info$latent_variables, drop = FALSE]))

  # --- Compute requested quantities -------------------------------------------
  need_star <- "star" %in% what || "prob" %in% what || "exp" %in% what || "class" %in% what
  need_prob <- "prob" %in% what || "exp" %in% what || "class" %in% what
  need_exp  <- "exp"  %in% what
  need_cls  <- "class"%in% what

  ystar_list <- vector("list", n_groups)
  probs_list <- vector("list", n_groups)
  yexp_list  <- vector("list", n_groups)
  class_list <- vector("list", n_groups)
  names(ystar_list) <- names(probs_list) <- names(yexp_list) <- names(class_list) <- group_labels

  for (g in seq_len(n_groups)) {
    lam_g <- est_list[[g]]$lambda
    nu_g  <- est_list[[g]]$nu %||% numeric(0)
    al    <- .align_Lambda_nu(lam_g, nu_g, ov_subset = ov_ord,
                              ov_all = ov_all, latent_order = info$latent_variables)

    eta_g <- eta_list[[g]]

    # y* latent predictor
    ystar_g <- sweep(eta_g %*% t(al$lambda), 2, al$nu, `+`)  # N x n_ord
    colnames(ystar_g) <- ov_ord

    # theta parametrization: divide both y* and thresholds by residual sd
    tau_vec <- th_list[[g]]
    if (identical(param, "theta")) {
      theta_g <- est_list[[g]]$theta
      if (!is.null(theta_g)) {
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
      ystar_g <- sweep(ystar_g, 2, sig_eps, "/")
    } else {
      sig_eps <- rep(1, length(ov_ord))
    }

    # Store y* if needed
    if (need_star) {
      ystar_df <- tibble::as_tibble(ystar_g)
      names(ystar_df) <- paste0(prefix_ystar, names(ystar_df))
      ystar_list[[g]] <- ystar_df
    } else {
      ystar_list[[g]] <- tibble::tibble()
    }

    # Probabilities (+ E[Y], class) if needed
    if (need_prob) {
      N_g <- nrow(ystar_g)

      # Prepare holders for this group
      probs_g <- NULL
      yexp_g  <- NULL
      class_g <- NULL

      # Iterate ordinal indicators
      for (jj in seq_along(ov_ord)) {
        j <- ov_ord[jj]
        ystar_j <- ystar_g[, jj]

        lev <- .levels_for_item(fs_list[[g]], tau_vec, j)
        if (!length(lev)) {
          warning("Levels for '", j, "' not found; skipping probabilities for this item (group ", g, ").")
          next
        }
        m <- length(lev)

        tau_j <- .get_tau_for(tau_vec, j)
        if (identical(param, "theta") && length(tau_j)) {
          tau_j <- tau_j / sig_eps[jj]
        }
        if (length(tau_j) != (m - 1L)) {
          warning("Thresholds for '", j, "' mismatch (have ", length(tau_j),
                  ", need ", m - 1L, "); skipping probabilities for this item (group ", g, ").")
          next
        }
        tau_pad <- c(-Inf, as.numeric(tau_j), Inf)

        # P(Y=c) = Phi(tau_c - y*) - Phi(tau_{c-1} - y*), row-wise
        probs_j <- matrix(NA_real_, nrow = N_g, ncol = m)
        for (ii in seq_len(N_g)) {
          diffs <- stats::pnorm(tau_pad - ystar_j[ii])
          probs_j[ii, ] <- diffs[-1L] - diffs[-length(diffs)]
        }
        colnames(probs_j) <- paste0(".pr_", make.names(lev, allow_ = TRUE), "__", j)

        # Expected score on ordinal codes 1..m
        if (need_exp) {
          codes <- seq_len(m)
          yexp_j <- as.numeric(probs_j %*% codes)
        }

        # Most likely class (label)
        if (need_cls) {
          ml_idx <- max.col(probs_j, ties.method = "first")
          ml_lab <- lev[ml_idx]
        }

        # Bind per item (only what is requested)
        if ("prob" %in% what) {
          probs_g <- if (is.null(probs_g)) {
            as.data.frame(probs_j, check.names = FALSE)
          } else {
            cbind(probs_g, as.data.frame(probs_j, check.names = FALSE))
          }
        }
        if (need_exp) {
          nm <- paste0(".yexp_", j)
          if (is.null(yexp_g)) yexp_g <- setNames(data.frame(yexp_j), nm) else yexp_g[[nm]] <- yexp_j
        }
        if (need_cls) {
          nm <- paste0(".class_", j)
          fac <- factor(ml_lab, levels = lev)
          if (is.null(class_g)) class_g <- setNames(data.frame(fac), nm) else class_g[[nm]] <- fac
        }
      }

      probs_list[[g]] <- if (!is.null(probs_g)) tibble::as_tibble(probs_g) else tibble::tibble()
      yexp_list[[g]]  <- if (!is.null(yexp_g))  tibble::as_tibble(yexp_g)  else tibble::tibble()
      class_list[[g]] <- if (!is.null(class_g)) tibble::as_tibble(class_g) else tibble::tibble()
    } else {
      probs_list[[g]] <- tibble::tibble()
      yexp_list[[g]]  <- tibble::tibble()
      class_list[[g]] <- tibble::tibble()
    }
  }

  # --- Bind outputs with explicit numeric .gid --------------------------------
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
  probs_tbl <- dplyr::bind_rows(probs_list)
  yexp_tbl  <- dplyr::bind_rows(yexp_list)
  class_tbl <- dplyr::bind_rows(class_list)

  out <- dplyr::bind_cols(fs_tbl, ystar_tbl, probs_tbl, yexp_tbl, class_tbl)

  # --- Ensure plain doubles for numeric columns -------------------------------
  out <- dplyr::mutate(out, dplyr::across(dplyr::where(is.numeric), as.double))

  tibble::as_tibble(out)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
