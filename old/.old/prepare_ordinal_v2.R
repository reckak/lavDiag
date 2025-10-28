#' Prepare smooth latent grids + model-based item curves for **ordinal** indicators
#'
#' @description
#' Builds a synthetic, smooth grid of factor scores for each group and latent.
#' For every latent k and group g, it takes the empirical range of factor scores
#' (estimated via lavPredict), creates an evenly spaced sequence of length
#' `length.out`, and sets all *other* latents to zero (or to their empirical
#' group means). For each ordinal indicator–latent pair (j, k), it computes
#' one-latent-at-a-time expected scores E[Y_j] and delta-method CIs using
#' `vcov(fit)`. The output contains no original observed indicators — only the
#' synthetic latent grid and model-based curves.
#'
#' @param fit A converged `lavaan` object with an ordinal part (at least one `ov.ord`).
#' @param info Optional list returned by `model_info()`.
#' @param data Optional precomputed factor-score data (as returned by
#'   `lavPredict_parallel(fit, return_type = "data", se = TRUE)`), or a per-group
#'   list of such data.frames. Used only to extract FS ranges and (optionally)
#'   group-wise summaries of FS SEs; original OV are never used here.
#' @param level CI level (default 0.95).
#' @param vcov_type Optional `lavaan::vcov(type=)` override; `NULL` = default.
#' @param length.out Integer; grid length per *(group × latent)* block (default 100).
#' @param other_latents `"zero"` (default) or `"mean"`; how to fill non-target latents.
#' @param latent_var_as_factor If TRUE (default), `.latent_var` is an ordered factor.
#' @param se If TRUE (default), attach group-wise summaries of FS SEs to the grid.
#' @param se_summary `"median"` (default) or `"mean"`; summary for FS SEs.
#'
#' @return A tibble with:
#'   `.rid`, `.gid`, `.group`, `.latent_var`,
#'   latent columns (grid), optional FS SE columns (prefix `.se_<latent>`),
#'   and per-(indicator×latent) model curves:
#'   `m_est_<item>_<latent>`, `m_lwr_*`, `m_upr_*`.
#' @export
prepare_ordinal_old <- function(fit,
                            info                  = NULL,
                            data                  = NULL,
                            level                 = 0.95,
                            vcov_type             = NULL,
                            length.out            = 100L,
                            other_latents         = c("zero", "mean"),
                            latent_var_as_factor  = TRUE,
                            se                    = TRUE,
                            se_summary            = c("median","mean")) {

  # -- Assertions --------------------------------------------------------------
  .assert_lavaan_fit(
    fit,
    require_converged     = TRUE,
    require_latent        = TRUE,
    require_meanstructure = NA,         # meanstructure not required for ordinals
    require_measurement   = "lambda",   # loadings needed; nu for ordinals is 0
    forbid_multilevel     = TRUE
  )

  other_latents <- match.arg(other_latents, c("zero","mean"))
  se_summary    <- match.arg(se_summary, c("median","mean"))
  stopifnot(is.numeric(length.out), length(length.out) == 1L, is.finite(length.out))
  length.out <- as.integer(length.out)
  if (length.out < 2L) rlang::abort("`length.out` must be >= 2.")

  prefix_se <- ".se_"

  # -- Model info --------------------------------------------------------------
  if (is.null(info)) info <- model_info(fit)
  ov_all       <- info$observed_variables
  ov_ord       <- info$ov_ordinal
  latent_names <- info$latent_variables
  n_groups     <- info$n_groups %||% 1L
  gvar         <- info$group_var %||% ".group"
  g_labels     <- info$group_labels %||% as.character(seq_len(n_groups))
  param        <- info$parameterization %||% "delta"

  if (length(ov_ord) == 0L) rlang::abort("No ordinal indicators in the model.")
  if (length(latent_names) == 0L) rlang::abort("No latent variables detected.")

  # -- Factor scores only for ranges (and optional SE summaries) ---------------
  fs_df <- if (is.null(data)) {
    lavPredict_parallel(fit, return_type = "data",
                        se = isTRUE(se), prefix_se = prefix_se)
  } else data

  # Mini-normalization: allow `data` to be a per-group list
  if (is.list(fs_df) && !is.data.frame(fs_df)) {
    fs_list <- fs_df
    fs_df <- dplyr::bind_rows(lapply(seq_along(fs_list), function(i) {
      df_g <- tibble::as_tibble(fs_list[[i]])
      gid  <- if (length(names(fs_list))) match(names(fs_list)[i], g_labels) else i
      if (is.na(gid)) gid <- i
      if (!".gid" %in% names(df_g))   df_g$.gid   <- gid
      if (!".group" %in% names(df_g)) df_g$.group <- g_labels[gid]
      df_g
    }))
  }

  if (n_groups > 1L && !gvar %in% names(fs_df)) {
    if (".group" %in% names(fs_df)) gvar <- ".group" else
      rlang::abort(paste0("Expected group column '", gvar, "' not found in factor-score data."))
  }

  # Ensure .group/.gid are present
  if (n_groups > 1L) {
    fs_df$.group <- as.character(fs_df[[gvar]])
    fs_df$.gid   <- as.integer(factor(fs_df$.group, levels = g_labels))
  } else {
    fs_df$.group <- g_labels[1L]
    fs_df$.gid   <- 1L
  }

  # -- Range helper ------------------------------------------------------------
  rng_or_fallback <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(c(-3, 3))
    r <- range(x)
    if (!is.finite(r[1]) || !is.finite(r[2]) || r[1] == r[2]) c(-3, 3) else r
  }

  # -- Means for "other_latents = 'mean'" -------------------------------------
  means_by_g <- NULL
  if (identical(other_latents, "mean")) {
    means_by_g <- fs_df |>
      dplyr::group_by(.gid) |>
      dplyr::summarise(dplyr::across(dplyr::all_of(latent_names), ~mean(.x[is.finite(.x)]),
                                     .names = "{.col}"),
                       .groups = "drop")
  }

  # -- (Optional) group-wise summaries of FS SEs -------------------------------
  se_cols_exist <- isTRUE(se) && all(paste0(prefix_se, latent_names) %in% names(fs_df))
  se_summ <- NULL
  if (se_cols_exist) {
    fun <- if (se_summary == "median") stats::median else mean
    se_summ <- fs_df |>
      dplyr::group_by(.gid) |>
      dplyr::summarise(
        dplyr::across(dplyr::all_of(paste0(prefix_se, latent_names)),
                      ~fun(.x[is.finite(.x)]), .names = "{.col}"),
        .groups = "drop"
      )
  }

  # -- Build the synthetic grid ------------------------------------------------
  grid_list <- vector("list", n_groups * length(latent_names))
  idx <- 0L
  for (g in seq_len(n_groups)) {
    base_vec <- setNames(numeric(length(latent_names)), latent_names)
    if (!is.null(means_by_g)) {
      mu_row <- means_by_g[means_by_g$.gid == g, , drop = FALSE]
      if (nrow(mu_row) == 1L) for (ln in latent_names) base_vec[[ln]] <- as.numeric(mu_row[[ln]])
    }
    fs_g <- fs_df[fs_df$.gid == g, latent_names, drop = FALSE]
    for (k in seq_along(latent_names)) {
      kname <- latent_names[k]
      r_gk  <- rng_or_fallback(fs_g[[kname]])
      seq_k <- seq(from = r_gk[1], to = r_gk[2], length.out = length.out)
      mat   <- matrix(rep(base_vec, each = length.out), ncol = length(base_vec))
      colnames(mat) <- names(base_vec)
      mat[, kname]  <- seq_k

      tb <- tibble::tibble(
        .gid        = g,
        .group      = g_labels[g],
        .latent_var = if (isTRUE(latent_var_as_factor)) factor(kname, levels = latent_names) else kname,
        !!!as.data.frame(mat, stringsAsFactors = FALSE)
      )

      if (!is.null(se_summ)) {
        se_row <- se_summ[se_summ$.gid == g, , drop = FALSE]
        if (nrow(se_row) == 1L) {
          for (ln in latent_names) {
            nm <- paste0(prefix_se, ln)
            if (nm %in% names(se_row)) tb[[nm]] <- rep(as.numeric(se_row[[nm]]), length.out)
          }
        }
      }

      idx <- idx + 1L
      grid_list[[idx]] <- tb
    }
  }

  out <- dplyr::bind_rows(grid_list)
  out$.rid <- seq_len(nrow(out))
  out <- out |>
    dplyr::relocate(.rid, .before = 1) |>
    dplyr::relocate(.gid, .after = .rid) |>
    dplyr::relocate(.group, .after = .gid) |>
    dplyr::relocate(.latent_var, .after = .group) |>
    dplyr::relocate(dplyr::all_of(latent_names), .after = ".latent_var")

  if (!is.null(se_summ)) {
    se_cols <- paste0(prefix_se, latent_names); se_cols <- se_cols[se_cols %in% names(out)]
    if (length(se_cols)) out <- dplyr::relocate(out, dplyr::all_of(se_cols),
                                                .after = dplyr::all_of(latent_names[length(latent_names)]))
  }

  # -- Measurement pieces & thresholds for delta method ------------------------
  est   <- lavaan::lavInspect(fit, "est")
  est_l <- if (is.list(est) && "lambda" %in% names(est)) list(est) else est
  th0   <- tryCatch(lavaan::lavInspect(fit, "th"), error = function(e) NULL)
  th_l  <- if (is.null(th0)) replicate(n_groups, setNames(numeric(0), character(0)), simplify = FALSE)
  else if (is.list(th0)) th0 else list(th0)
  if (length(th_l) == 1L && n_groups > 1L) th_l <- rep(th_l, n_groups)

  V <- tryCatch({ if (is.null(vcov_type)) lavaan::vcov(fit) else lavaan::vcov(fit, type = vcov_type) },
                error = function(e) NULL)
  if (is.null(V)) rlang::abort("vcov(fit) is not available; cannot compute delta-method CIs for ordinals.")
  coef_names <- names(lavaan::coef(fit))
  pt <- lavaan::parTable(fit)
  zcrit <- stats::qnorm(1 - (1 - level)/2)

  # Helpers reused from augment_ordinal() logic -------------------------------
  .align_Lambda_nu <- function(lam_g, nu_g, ov_subset, ov_all, latent_order) {
    if (!is.null(colnames(lam_g))) lam_g <- lam_g[, latent_order, drop = FALSE]
    if (!is.null(rownames(lam_g))) lam_g <- lam_g[ov_subset, , drop = FALSE]
    else                           lam_g <- lam_g[match(ov_subset, ov_all), , drop = FALSE]
    if (is.null(nu_g) || length(nu_g) == 0L) { nu_g <- rep(0, length(ov_subset)); names(nu_g) <- ov_subset }
    else if (!is.null(names(nu_g))) nu_g <- nu_g[ov_subset]
    else                            nu_g <- nu_g[match(ov_subset, ov_all)]
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

  # -- Compute curves per group, item, latent ---------------------------------
  out2 <- out
  for (g in seq_len(n_groups)) {
    est_g <- est_l[[g]]
    lam_g <- est_g$lambda
    nu_g  <- est_g$nu %||% numeric(0)
    al    <- .align_Lambda_nu(lam_g, nu_g, ov_subset = ov_ord, ov_all = ov_all, latent_order = latent_names)

    # theta parametrization: residual SDs for y* (used to scale z)
    sig_eps <- NULL
    if (identical(param, "theta")) {
      theta_g <- est_g$theta
      if (!is.null(theta_g)) {
        if (!is.null(rownames(theta_g))) sig_eps <- sqrt(diag(theta_g)[ov_ord])
        else                             sig_eps <- sqrt(diag(theta_g)[match(ov_ord, ov_all)])
        sig_eps[!is.finite(sig_eps) | sig_eps <= 0] <- 1
      } else sig_eps <- rep(1, length(ov_ord))
    }

    # convenience: baseline eta for non-focal latents
    base_eta <- setNames(rep(0, length(latent_names)), latent_names)
    if (identical(other_latents, "mean") && !is.null(means_by_g)) {
      mu_row <- means_by_g[means_by_g$.gid == g, , drop = FALSE]
      if (nrow(mu_row) == 1L) for (ln in latent_names) base_eta[[ln]] <- as.numeric(mu_row[[ln]])
    }

    for (jj in seq_along(ov_ord)) {
      j     <- ov_ord[jj]
      tau_j <- as.numeric(.get_tau_for(th_l[[g]], j))
      if (!length(tau_j)) next
      m       <- length(tau_j) + 1L
      scores  <- seq_len(m)        # numeric scoring 1..m (safe fallback)
      tau_pad <- c(-Inf, tau_j, Inf)

      # pre-extract row of loadings for j
      lam_row <- al$lambda[jj, , drop = TRUE]
      nz_idx  <- which(is.finite(lam_row) & lam_row != 0)
      if (!length(nz_idx)) next

      # constant offset from non-focal latents held at baseline (for y*):
      mu_const_all <- sum(lam_row[nz_idx] * base_eta[nz_idx])

      for (k_idx in nz_idx) {
        k_name <- latent_names[k_idx]
        rows_k <- which(out2$.gid == g & out2$.latent_var == k_name)
        eta_k  <- out2[[k_name]][rows_k]

        # y* = nu_j + lambda_{jk} * eta_k + sum_{r!=k} lambda_{jr} * base_eta_r
        ystar <- (al$nu[jj] %||% 0) +
          lam_row[k_idx] * eta_k +
          (mu_const_all - lam_row[k_idx] * base_eta[k_name])

        # z = (tau - y*) / s  where s = 1 for "delta", or s = sigma_eps_j for "theta"
        if (identical(param, "theta")) {
          se_j <- if (is.null(sig_eps)) 1 else sig_eps[jj]
          zmat <- outer(ystar, tau_pad, function(y, t) (t - y) / se_j)
        } else {
          zmat <- outer(ystar, tau_pad, function(y, t) (t - y))
          se_j <- 1
        }

        Phi     <- stats::pnorm(zmat)
        probs_j <- Phi[, -1, drop = FALSE] - Phi[, -ncol(Phi), drop = FALSE]  # N x m
        mu_vec  <- as.numeric(probs_j %*% scores)

        # --- Delta-method gradient wrt (nu_j, lambda_{jk}, thresholds tau_{j,*}, and possibly other loadings at baseline=mean)
        phi_mat   <- stats::dnorm(zmat)                    # N x (m+2)
        phi_inner <- phi_mat[, -1, drop = FALSE]
        phi_prev  <- phi_mat[, -ncol(phi_mat), drop = FALSE]

        # dμ/dz = -B, with B = sum_c (score_c)*(φ(τ_c - z) - φ(τ_{c-1} - z))
        B <- as.numeric((phi_inner[, 1:m, drop = FALSE] - phi_prev[, 1:m, drop = FALSE]) %*% matrix(scores, ncol = 1))

        d_nu   <- (-1 / se_j) * B                      # ∂μ/∂ν_j  (ν_j=0 for ordinals; keep for completeness)
        d_lamk <- (-eta_k / se_j) * B                  # ∂μ/∂λ_{jk}
        # ∂μ/∂τ_{j,c} = (score_c - score_{c+1}) * φ( (τ_c - y*) / s ) * (1/s)
        d_tau  <- matrix(0, nrow = length(rows_k), ncol = m - 1L)
        if (m > 1L) {
          score_diff <- scores[1:(m-1)] - scores[2:m]
          d_tau <- sweep(phi_mat[, 2:(ncol(phi_mat)-1), drop = FALSE], 2, score_diff, `*`) * (1 / se_j)
        }

        # Include extra uncertainty from other loadings if baseline = "mean"
        extra_ids <- integer(0)
        extra_A   <- NULL
        if (identical(other_latents, "mean")) {
          other_r <- setdiff(nz_idx, k_idx)
          if (length(other_r)) {
            # columns correspond to λ_{jr}, derivative = (∂μ/∂z)*(∂z/∂λ_{jr}) = (-B)*(base_eta_r)/s
            extra_A   <- do.call(cbind, lapply(other_r, function(r) (-B / se_j) * rep(base_eta[r], length(rows_k))))
            names(extra_A) <- paste0("lam__", j, "__", latent_names[other_r])
          }
        }

        # -- Map derivatives to free-parameter indices in vcov -----------------
        # indices for nu_j, lambda_{jk}, thresholds τ_{j,*}
        get_id_nu <- function(g, jname) {
          ii <- which(pt$lhs == jname & pt$op == "~1" & pt$group == g & pt$free > 0)
          if (length(ii)) pt$free[ii][1] else 0L
        }
        get_id_lam <- function(g, jname, kname) {
          ii <- which(pt$lhs == kname & pt$op == "=~" & pt$rhs == jname & pt$group == g & pt$free > 0)
          if (length(ii)) pt$free[ii][1] else 0L
        }
        get_id_tau <- function(g, jname, m) {
          vapply(seq_len(m-1L), function(ti) {
            ii <- which(pt$lhs == jname & pt$op == "|" & pt$rhs == paste0("t", ti) & pt$group == g & pt$free > 0)
            if (length(ii)) pt$free[ii][1] else 0L
          }, integer(1))
        }

        id_nu   <- get_id_nu(g, j)
        id_lamk <- get_id_lam(g, j, k_name)
        id_tau  <- get_id_tau(g, j, m)

        # Build gradient matrix A (rows = grid rows_k; cols = selected free params, de-duplicated)
        A_core <- cbind(
          if (id_nu   > 0L) d_nu  else NULL,
          if (id_lamk > 0L) d_lamk else NULL,
          if (any(id_tau > 0L))    d_tau[, id_tau > 0L, drop = FALSE] else NULL,
          extra_A
        )
        cols_all <- c(if (id_nu > 0L) id_nu else integer(0),
                      if (id_lamk > 0L) id_lamk else integer(0),
                      id_tau[id_tau > 0L],
                      # map extra loadings to their indices
                      if (!is.null(extra_A)) {
                        sapply(setdiff(nz_idx, k_idx), function(r) {
                          get_id_lam(g, j, latent_names[r])
                        })
                      } else integer(0))

        if (length(cols_all)) {
          uniq <- unique(cols_all)
          pos  <- match(uniq, cols_all)
          A    <- A_core[, pos, drop = FALSE]
          Vsub <- V[uniq, uniq, drop = FALSE]
          se_mu <- sqrt(pmax(rowSums((A %*% Vsub) * A), 0))
        } else {
          se_mu <- rep(0, length(rows_k))
        }

        base  <- paste0(j, "_", k_name)
        nm_e  <- paste0("m_est_", base)
        nm_l  <- paste0("m_lwr_", base)
        nm_u  <- paste0("m_upr_", base)

        out2[rows_k, nm_e] <- mu_vec
        out2[rows_k, nm_l] <- mu_vec - zcrit * se_mu
        out2[rows_k, nm_u] <- mu_vec + zcrit * se_mu
      } # k_idx
    }   # jj
  }     # g

  # numeric hygiene
  is_num   <- vapply(out2, is.numeric, logical(1L))
  keep_int <- names(out2) %in% c(".gid", ".rid")
  to_dbl   <- is_num & !keep_int
  for (nm in names(out2)[to_dbl]) out2[[nm]] <- as.double(out2[[nm]])

  tibble::as_tibble(out2)
}
