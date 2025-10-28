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
                            se_summary            = c("median","mean"),
                            plan                  = c("auto","none","multisession","multicore","sequential","cluster"),
                            workers               = NULL,
                            cluster               = NULL) {

  # --- Setup & assertions ----------------------------------------------------
  .assert_lavaan_fit(
    fit,
    require_converged     = TRUE,
    require_latent        = TRUE,
    require_meanstructure = NA,
    require_measurement   = "lambda",
    forbid_multilevel     = TRUE
  )

  other_latents <- match.arg(other_latents, c("zero","mean"))
  se_summary    <- match.arg(se_summary, c("median","mean"))
  plan          <- match.arg(plan)
  stopifnot(is.numeric(length.out), length(length.out) == 1L, is.finite(length.out))
  length.out <- as.integer(length.out)
  if (length.out < 2L) rlang::abort("`length.out` must be >= 2.")

  prefix_se <- ".se_"

  # --- Model info ------------------------------------------------------------
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

  # --- Factor scores (allow list-per-group) ----------------------------------
  fs_df <- if (is.null(data)) {
    lavPredict_parallel(fit, return_type = "data", se = isTRUE(se), prefix_se = prefix_se)
  } else data

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

  if (n_groups > 1L) {
    fs_df$.group <- as.character(fs_df[[gvar]])
    fs_df$.gid   <- as.integer(factor(fs_df$.group, levels = g_labels))
  } else {
    fs_df$.group <- g_labels[1L]
    fs_df$.gid   <- 1L
  }

  # --- Range helper ----------------------------------------------------------
  rng_or_fallback <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(c(-3, 3))
    r <- range(x)
    if (!is.finite(r[1]) || !is.finite(r[2]) || r[1] == r[2]) c(-3, 3) else r
  }

  # Means for non-focal latents when other_latents == "mean"
  means_by_g <- NULL
  if (identical(other_latents, "mean")) {
    means_by_g <- fs_df |>
      dplyr::group_by(.gid) |>
      dplyr::summarise(dplyr::across(dplyr::all_of(latent_names), ~mean(.x[is.finite(.x)]), .names = "{.col}"), .groups = "drop")
  }

  # Optional group-wise summaries of FS SEs (relaxed: use any present SE cols)
  present_se <- intersect(paste0(prefix_se, latent_names), names(fs_df))
  se_summ <- NULL
  if (isTRUE(se) && length(present_se) > 0L) {
    fun <- if (se_summary == "median") stats::median else mean
    se_summ <- fs_df |>
      dplyr::group_by(.gid) |>
      dplyr::summarise(
        dplyr::across(dplyr::all_of(present_se), ~fun(.x[is.finite(.x)]), .names = "{.col}"),
        .groups = "drop"
      )
  }

  # --- Build the synthetic grid (vectorized) ---------------------------------
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
    se_cols <- paste0(prefix_se, latent_names); se_cols <- se_cols %in% names(out)
    if (any(se_cols)) out <- dplyr::relocate(out, dplyr::all_of(paste0(prefix_se, latent_names)[se_cols]), .after = dplyr::all_of(latent_names[length(latent_names)]))
  }

  # Pre-index rows per (g, k) once (used many times below) — optimized with integer codes
  if (is.factor(out$.latent_var)) {
    lv_levels <- levels(out$.latent_var)
    lv_int    <- as.integer(out$.latent_var)
    rows_idx <- lapply(seq_len(n_groups), function(g) {
      rg <- out$.gid == g
      setNames(lapply(latent_names, function(kname) {
        k_int <- match(kname, lv_levels)
        which(rg & lv_int == k_int)
      }), latent_names)
    })
  } else {
    # Fallback when `.latent_var` is not a factor
    rows_idx <- lapply(seq_len(n_groups), function(g) {
      rg <- out$.gid == g
      setNames(lapply(latent_names, function(kname) which(rg & out$.latent_var == kname)), latent_names)
    })
  }

  # --- Measurement pieces for delta method -----------------------------------
  est   <- lavaan::lavInspect(fit, "est")
  est_l <- if (is.list(est) && "lambda" %in% names(est)) list(est) else est
  th0   <- tryCatch(lavaan::lavInspect(fit, "th"), error = function(e) NULL)
  th_l  <- if (is.null(th0)) replicate(n_groups, setNames(numeric(0), character(0)), simplify = FALSE) else if (is.list(th0)) th0 else list(th0)
  if (length(th_l) == 1L && n_groups > 1L) th_l <- rep(th_l, n_groups)

  V <- tryCatch({ if (is.null(vcov_type)) lavaan::vcov(fit) else lavaan::vcov(fit, type = vcov_type) }, error = function(e) NULL)
  if (is.null(V)) rlang::abort("vcov(fit) is not available; cannot compute delta-method CIs for ordinals.")
  pt <- lavaan::parTable(fit)
  zcrit <- stats::qnorm(1 - (1 - level)/2)

  # Numerical safety: ensure free-parameter IDs do not exceed vcov dimension
  npar_V   <- ncol(V)
  max_free <- suppressWarnings(max(pt$free[pt$free > 0], 0L))
  if (is.finite(max_free) && max_free > npar_V) {
    rlang::warn("Some free parameter indices exceed vcov dimension; those IDs will be ignored in SE computation.")
  }

  # Align lambda/nu per group once; precompute residual SDs if theta
  align_by_g <- vector("list", n_groups)
  sig_eps_by_g <- vector("list", n_groups)
  for (g in seq_len(n_groups)) {
    est_g <- est_l[[g]]
    lam_g <- est_g$lambda
    nu_g  <- est_g$nu %||% numeric(0)

    if (!is.null(colnames(lam_g))) lam_g <- lam_g[, latent_names, drop = FALSE]
    if (!is.null(rownames(lam_g))) lam_g <- lam_g[ov_ord, , drop = FALSE] else lam_g <- lam_g[match(ov_ord, ov_all), , drop = FALSE]
    if (is.null(nu_g) || length(nu_g) == 0L) { nu_g <- rep(0, length(ov_ord)); names(nu_g) <- ov_ord }
    else if (!is.null(names(nu_g))) nu_g <- nu_g[ov_ord] else nu_g <- nu_g[match(ov_ord, ov_all)]

    align_by_g[[g]] <- list(lambda = lam_g, nu = nu_g)

    if (identical(param, "theta")) {
      theta_g <- est_g$theta
      if (!is.null(theta_g)) {
        if (!is.null(rownames(theta_g))) sig_eps <- sqrt(diag(theta_g)[ov_ord]) else sig_eps <- sqrt(diag(theta_g)[match(ov_ord, ov_all)])
        sig_eps[!is.finite(sig_eps) | sig_eps <= 0] <- 1
      } else sig_eps <- rep(1, length(ov_ord))
      sig_eps_by_g[[g]] <- sig_eps
    } else {
      sig_eps_by_g[[g]] <- rep(1, length(ov_ord))
    }
  }

  # Pre-map parameter IDs for all (g, j, k, tau) with clamping to `npar_V`
  id_map <- vector("list", n_groups)
  for (g in seq_len(n_groups)) {
    id_nu <- setNames(integer(length(ov_ord)), ov_ord)
    id_lam <- array(0L, dim = c(length(ov_ord), length(latent_names)), dimnames = list(ov_ord, latent_names))
    id_tau <- lapply(ov_ord, function(j) integer(0))

    for (jj in seq_along(ov_ord)) {
      j <- ov_ord[jj]
      ii <- which(pt$lhs == j & pt$op == "~1" & pt$group == g & pt$free > 0)
      id_val <- if (length(ii)) pt$free[ii][1] else 0L
      id_nu[j] <- if (id_val > 0L && id_val <= npar_V) id_val else 0L

      for (k in seq_along(latent_names)) {
        kname <- latent_names[k]
        ii <- which(pt$lhs == kname & pt$op == "=~" & pt$rhs == j & pt$group == g & pt$free > 0)
        id_val <- if (length(ii)) pt$free[ii][1] else 0L
        id_lam[j, kname] <- if (id_val > 0L && id_val <= npar_V) id_val else 0L
      }

      rows_tau <- which(pt$lhs == j & pt$op == "|" & pt$group == g & pt$free > 0)
      ids_tau  <- if (length(rows_tau)) pt$free[rows_tau] else integer(0)
      if (length(ids_tau)) ids_tau <- ids_tau[ids_tau > 0L & ids_tau <= npar_V]
      id_tau[[j]] <- ids_tau
    }

    id_map[[g]] <- list(id_nu = id_nu, id_lam = id_lam, id_tau = id_tau)
  }

  # Baseline means for non-focal latents per group
  base_eta_by_g <- lapply(seq_len(n_groups), function(g) {
    base <- setNames(rep(0, length(latent_names)), latent_names)
    if (identical(other_latents, "mean") && !is.null(means_by_g)) {
      mu_row <- means_by_g[means_by_g$.gid == g, , drop = FALSE]
      if (nrow(mu_row) == 1L) for (ln in latent_names) base[[ln]] <- as.numeric(mu_row[[ln]])
    }
    base
  })

  # Helper to fetch thresholds per item quickly
  get_tau_for <- function(th_vec, j) {
    nm <- names(th_vec); if (is.null(nm)) return(numeric(0))
    sel <- startsWith(nm, paste0(j, "|"))
    tau <- th_vec[sel]
    if (!length(tau)) return(numeric(0))
    idx <- suppressWarnings(as.integer(gsub("^.*\\|t(\\d+)$", "\\1", names(tau))))
    idx[is.na(idx)] <- seq_along(tau)
    as.numeric(tau[order(idx)])
  }

  # --- Parallelization over (group, item) tasks ------------------------------
  reset_plan <- NULL
  if (plan != "none") {
    reset_plan <- .set_future_plan(plan = plan, workers = workers, cluster = cluster)
    on.exit({ if (!is.null(reset_plan)) reset_plan() }, add = TRUE)
  }

  tasks <- expand.grid(g = seq_len(n_groups), jj = seq_along(ov_ord), KEEP.OUT.ATTRS = FALSE)

  compute_one <- function(g, jj) {
    j        <- ov_ord[jj]
    lam_row  <- align_by_g[[g]]$lambda[jj, , drop = TRUE]
    nz_idx   <- which(is.finite(lam_row) & lam_row != 0)
    if (!length(nz_idx)) return(NULL)

    tau_j <- get_tau_for(th_l[[g]], j)
    if (!length(tau_j)) return(NULL)
    m       <- length(tau_j) + 1L
    scores  <- seq_len(m)
    tau_pad <- c(-Inf, tau_j, Inf)

    nu_j    <- (align_by_g[[g]]$nu[jj] %||% 0)
    sig_eps <- sig_eps_by_g[[g]][jj]

    base_eta <- base_eta_by_g[[g]]
    mu_const_all <- sum(lam_row[nz_idx] * base_eta[nz_idx])

    out_parts <- vector("list", length(nz_idx))

    for (i in seq_along(nz_idx)) {
      k_idx  <- nz_idx[i]
      k_name <- latent_names[k_idx]
      rows_k <- rows_idx[[g]][[k_name]]
      eta_k  <- out[[k_name]][rows_k]

      # y* on the grid for the focal latent k
      ystar <- nu_j + lam_row[k_idx] * eta_k + (mu_const_all - lam_row[k_idx] * base_eta[k_name])

      # z = (tau - y*) / s
      if (identical(param, "theta")) {
        zmat <- outer(ystar, tau_pad, function(y, t) (t - y) / sig_eps)
        inv_s <- 1 / sig_eps
      } else {
        zmat <- outer(ystar, tau_pad, function(y, t) (t - y))
        inv_s <- 1
      }

      Phi     <- stats::pnorm(zmat)
      probs_j <- Phi[, -1, drop = FALSE] - Phi[, -ncol(Phi), drop = FALSE]
      mu_vec  <- as.numeric(probs_j %*% scores)

      # Delta-method pieces
      phi_mat   <- stats::dnorm(zmat)
      phi_inner <- phi_mat[, -1, drop = FALSE]
      phi_prev  <- phi_mat[, -ncol(phi_mat), drop = FALSE]
      B <- as.numeric((phi_inner[, 1:m, drop = FALSE] - phi_prev[, 1:m, drop = FALSE]) %*% matrix(scores, ncol = 1))

      d_nu   <- (-inv_s)   * B
      d_lamk <- (-eta_k*inv_s) * B
      d_tau  <- matrix(0, nrow = length(rows_k), ncol = m - 1L)
      if (m > 1L) {
        score_diff <- scores[1:(m-1)] - scores[2:m]
        d_tau <- sweep(phi_mat[, 2:(ncol(phi_mat)-1), drop = FALSE], 2, score_diff, `*`) * inv_s
      }

      # Gradient assembly
      A_mats <- list(); A_ids <- integer(0)
      id_nu   <- id_map[[g]]$id_nu[j]
      id_lamk <- id_map[[g]]$id_lam[j, k_name]
      ids_tau_all <- id_map[[g]]$id_tau[[j]]
      ids_tau <- if (length(ids_tau_all)) ids_tau_all[seq_len(min(length(ids_tau_all), m-1))] else integer(0)

      if (id_nu > 0L)   { A_mats[[length(A_mats)+1L]] <- matrix(d_nu,   ncol = 1L); A_ids <- c(A_ids, id_nu) }
      if (id_lamk > 0L) { A_mats[[length(A_mats)+1L]] <- matrix(d_lamk, ncol = 1L); A_ids <- c(A_ids, id_lamk) }
      if (length(ids_tau)) { A_mats[[length(A_mats)+1L]] <- d_tau; A_ids <- c(A_ids, ids_tau) }

      if (identical(other_latents, "mean")) {
        other_r <- setdiff(nz_idx, k_idx)
        if (length(other_r)) {
          add_cols <- list(); add_ids <- integer(0)
          for (r in other_r) {
            id_lr <- id_map[[g]]$id_lam[j, latent_names[r]]
            if (id_lr > 0L) {
              add_cols[[length(add_cols)+1L]] <- (-B*inv_s) * base_eta[ latent_names[r] ]
              add_ids  <- c(add_ids, id_lr)
            }
          }
          if (length(add_ids)) {
            A_mats[[length(A_mats)+1L]] <- do.call(cbind, add_cols)
            A_ids  <- c(A_ids, add_ids)
          }
        }
      }

      se_mu <- if (length(A_ids)) {
        uniq <- unique(A_ids)
        pos  <- match(uniq, A_ids)
        A    <- do.call(cbind, A_mats)[, pos, drop = FALSE]
        Vsub <- V[uniq, uniq, drop = FALSE]
        sqrt(pmax(rowSums((A %*% Vsub) * A), 0))
      } else {
        rep(0, length(rows_k))
      }

      base  <- paste0(j, "_", k_name)
      nm_e  <- paste0("m_est_", base)
      nm_l  <- paste0("m_lwr_", base)
      nm_u  <- paste0("m_upr_", base)

      out_parts[[i]] <- list(rows = rows_k, nm_e = nm_e, nm_l = nm_l, nm_u = nm_u,
                             est = mu_vec, lwr = mu_vec - zcrit * se_mu, upr = mu_vec + zcrit * se_mu)
    }

    out_parts
  }

  # Choose runner; require future.apply if a parallel plan is requested
  if (plan == "none") {
    parts <- lapply(seq_len(nrow(tasks)), function(ii) compute_one(tasks$g[ii], tasks$jj[ii]))
  } else {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      rlang::abort("Parallel plan requested but 'future.apply' is not installed. Install it or set plan='none'.")
    }
    parts <- future.apply::future_lapply(seq_len(nrow(tasks)), function(ii) compute_one(tasks$g[ii], tasks$jj[ii]))
  }

  # --- Scoped writes: pre-scan names and precreate columns -------------------
  out2 <- out
  nm_e_all <- character(0); nm_l_all <- character(0); nm_u_all <- character(0)
  for (chunk in parts) {
    if (is.null(chunk)) next
    for (piece in chunk) {
      if (is.null(piece)) next
      nm_e_all <- c(nm_e_all, piece$nm_e)
      nm_l_all <- c(nm_l_all, piece$nm_l)
      nm_u_all <- c(nm_u_all, piece$nm_u)
    }
  }
  nm_e_all <- unique(nm_e_all); nm_l_all <- unique(nm_l_all); nm_u_all <- unique(nm_u_all)
  for (nm in nm_e_all) if (is.null(out2[[nm]])) out2[[nm]] <- rep(NA_real_, nrow(out2))
  for (nm in nm_l_all) if (is.null(out2[[nm]])) out2[[nm]] <- rep(NA_real_, nrow(out2))
  for (nm in nm_u_all) if (is.null(out2[[nm]])) out2[[nm]] <- rep(NA_real_, nrow(out2))

  use_vctrs <- requireNamespace("vctrs", quietly = TRUE)
  for (chunk in parts) {
    if (is.null(chunk)) next
    for (piece in chunk) {
      if (is.null(piece)) next
      if (use_vctrs) {
        out2[[piece$nm_e]] <- vctrs::vec_assign(out2[[piece$nm_e]], piece$rows, piece$est)
        out2[[piece$nm_l]] <- vctrs::vec_assign(out2[[piece$nm_l]], piece$rows, piece$lwr)
        out2[[piece$nm_u]] <- vctrs::vec_assign(out2[[piece$nm_u]], piece$rows, piece$upr)
      } else {
        out2[[piece$nm_e]][piece$rows] <- piece$est
        out2[[piece$nm_l]][piece$rows] <- piece$lwr
        out2[[piece$nm_u]][piece$rows] <- piece$upr
      }
    }
  }

  # numeric hygiene
  is_num   <- vapply(out2, is.numeric, logical(1L))
  keep_int <- names(out2) %in% c(".gid", ".rid")
  to_dbl   <- is_num & !keep_int
  for (nm in names(out2)[to_dbl]) out2[[nm]] <- as.double(out2[[nm]])

  tibble::as_tibble(out2)
}
