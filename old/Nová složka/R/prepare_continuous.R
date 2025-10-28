#' Prepare data for item response curves for continuous indicators
#'
#' Builds per-(item, factor) "one-latent-at-a-time" predictions:
#'   m_est_{item}_{factor} = nu_item + lambda_{item,factor} * eta_factor,
#' holding all other latent variables at their mean (0 when transform = FALSE).
#' Also attaches Wald CIs via a delta method using only {nu_item, lambda_{item,factor}}.
#'
#' Output contains: .rid, .gid, .group, original observed variables, factor scores,
#' and the m_est / m_lwr / m_upr columns for all (item, factor) with nonzero loading.
#'
#' @param fit   Fitted lavaan object (single-level).
#' @param data  Optional newdata. If NULL, uses fitted data. Factor scores are
#'              obtained via lavPredict_parallel(append.data = TRUE).
#' @param info  Optional result of model_info(fit) to avoid recomputation.
#' @param level Confidence level for Wald CIs (default 0.95).
#' @return tibble
#' @export
prepare_continuous <- function(fit, data = NULL, info = NULL, level = 0.95) {
  # -- Guardrails --------------------------------------------------------------
  .assert_lavaan_fit(
    fit,
    require_converged     = TRUE,
    require_meanstructure = TRUE,
    require_latent        = TRUE,
    forbid_multilevel     = TRUE
  )

  # -- Model metadata ----------------------------------------------------------
  if (is.null(info)) info <- model_info(fit)
  ov_all        <- info$observed_variables
  ov_cont       <- info$ov_continuous
  ov_ord        <- info$ov_ordinal
  latent_names  <- info$latent_variables
  n_groups      <- info$n_groups %||% 1L
  group_labels  <- info$group_labels %||% as.character(seq_len(n_groups))

  if (length(ov_cont) == 0L) {
    stop("No continuous observed indicators detected; nothing to prepare.", call. = FALSE)
  }
  if (length(ov_ord) > 0L) {
    warning("Model contains ordinal indicators. They are ignored here (continuous-only IR curves).")
  }

  # -- Factor scores + observed variables (single pass) ------------------------
  # If caller supplied 'data', use it; otherwise use fitted data.
  # We always ask for append.data = TRUE so observed columns travel with scores.
  fs_and_ov <- if (is.null(data)) {
    lavPredict_parallel(fit, return_type = "data")
  } else {
    # Pass user 'data' straight through to lavaan via wrapper
    lavPredict_parallel(fit, return_type = "data", newdata = data)
  }

  # Ensure .gid / .group exist and factor scores are present
  if (!(".gid" %in% names(fs_and_ov))) {
    fs_and_ov$.gid <- 1L
  }
  if (!(".group" %in% names(fs_and_ov))) {
    fs_and_ov$.group <- if (n_groups == 1L) group_labels[1L] else as.character(fs_and_ov$.gid)
  }

  # Check that factor-score columns are present
  eta_cols <- intersect(latent_names, names(fs_and_ov))
  if (!length(eta_cols)) {
    stop("No latent factor-score columns detected in lavPredict output.", call. = FALSE)
  }

  # -- Measurement pieces per group -------------------------------------------
  est      <- lavaan::lavInspect(fit, "est")
  est_list <- if (is.list(est) && "lambda" %in% names(est)) list(est) else est

  # -- VCOV + delta-method helpers --------------------------------------------
  V <- tryCatch(lavaan::vcov(fit), error = function(e) NULL)
  if (is.null(V)) stop("vcov(fit) is not available; cannot compute delta-method CIs.", call. = FALSE)
  PT         <- lavaan::parTable(fit)
  single_grp <- (n_groups == 1L)
  zcrit      <- stats::qnorm(1 - (1 - level)/2)

  # Resolve free-id accessors (stable across MG)
  get_free_lambda <- function(g, j, k) {
    if (single_grp) { w <- which(PT$op == "=~" & PT$lhs == k & PT$rhs == j) }
    else            { w <- which(PT$group == g & PT$op == "=~" & PT$lhs == k & PT$rhs == j) }
    if (!length(w)) return(0L); PT$free[w[1L]]
  }
  get_free_nu <- function(g, j) {
    if (single_grp) { w <- which(PT$op == "~1" & PT$lhs == j) }
    else            { w <- which(PT$group == g & PT$op == "~1" & PT$lhs == j) }
    if (!length(w)) return(0L); PT$free[w[1L]]
  }

  # -- Pre-allocate containers for all (item, factor) pairs --------------------
  # We will generate columns only for loadings != 0 (finite and nonzero).
  out <- fs_and_ov

  # For each group, compute predictions and CIs
  for (g in seq_len(n_groups)) {
    rows_g <- which(out$.gid == g)
    if (!length(rows_g)) next

    est_g <- est_list[[g]]
    lam_g <- est_g$lambda
    nu_g  <- est_g$nu %||% numeric(0)

    # Align Lambda rows to continuous ov and columns to latent order
    if (!is.null(colnames(lam_g))) lam_g <- lam_g[, latent_names, drop = FALSE]
    if (!is.null(rownames(lam_g))) {
      lam_g <- lam_g[ov_cont, , drop = FALSE]
    } else {
      lam_g <- lam_g[match(ov_cont, ov_all), , drop = FALSE]
    }

    # Align nu to continuous ov
    if (is.null(nu_g) || length(nu_g) == 0L) {
      nu_g <- rep(0, length(ov_cont)); names(nu_g) <- ov_cont
    } else if (!is.null(names(nu_g))) {
      nu_g <- nu_g[ov_cont]
    } else {
      nu_g <- nu_g[match(ov_cont, ov_all)]
    }

    # Factor score matrix for this group (keep order = latent_names)
    eta_mat <- as.matrix(out[rows_g, latent_names, drop = FALSE])

    # Loop items and their nonzero loadings
    for (j_idx in seq_along(ov_cont)) {
      j <- ov_cont[j_idx]
      lam_row <- lam_g[j_idx, , drop = TRUE]
      nz_idx  <- which(is.finite(lam_row) & lam_row != 0)

      if (!length(nz_idx)) next  # nothing to plot for this item

      for (k_idx in nz_idx) {
        k_name <- colnames(lam_g)[k_idx] %||% latent_names[k_idx]

        # Point prediction: m_est_{j}_{k} = nu_j + lambda_{jk} * eta_k
        mu_vec <- as.numeric(nu_g[j_idx] + lam_row[k_idx] * eta_mat[, k_idx])

        # Delta-method SE using only free IDs for {nu_j, lambda_{jk}}
        id_nu  <- get_free_nu(g, j)
        id_lam <- get_free_lambda(g, j, k_name)
        ids    <- c(if (id_nu > 0L) id_nu else integer(0),
                    if (id_lam > 0L) id_lam else integer(0))

        if (length(ids)) {
          uniq <- unique(ids)
          Vsub <- V[uniq, uniq, drop = FALSE]

          # Gradient per row: [d/d(nu)=1 , d/d(lambda)=eta_k]
          A <- matrix(0, nrow = length(rows_g), ncol = length(uniq))
          if (id_nu > 0L) {
            hit <- match(id_nu, uniq, nomatch = 0L)
            if (hit > 0L) A[, hit] <- 1
          }
          if (id_lam > 0L) {
            hit <- match(id_lam, uniq, nomatch = 0L)
            if (hit > 0L) A[, hit] <- eta_mat[, k_idx]
          }

          # g' V g per row
          se2 <- rowSums((A %*% Vsub) * A)
          se  <- sqrt(pmax(se2, 0))
        } else {
          se <- rep(0, length(rows_g))  # both fixed => no sampling uncertainty
        }

        # Column names
        base <- paste0(j, "_", k_name)
        nm_est <- paste0("m_est_", base)
        nm_lwr <- paste0("m_lwr_", base)
        nm_upr <- paste0("m_upr_", base)

        # Write into output
        out[rows_g, nm_est] <- mu_vec
        out[rows_g, nm_lwr] <- mu_vec - zcrit * se
        out[rows_g, nm_upr] <- mu_vec + zcrit * se
      }
    }
  }

  # Cast numeric columns to double (keep .gid as integer)
  is_num <- vapply(out, is.numeric, logical(1L))
  keep_int <- names(out) %in% c(".gid")
  to_double <- is_num & !keep_int
  for (nm in names(out)[to_double]) out[[nm]] <- as.double(out[[nm]])

  tibble::as_tibble(out)
}
