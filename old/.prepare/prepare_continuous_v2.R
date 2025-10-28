#' Prepare smooth latent grids + model-based item curves for **continuous** indicators
#' (Group-mapping fix: never trust incoming .gid; align strictly to lavaan's group order)
#'
#' @description
#' Builds a synthetic, smooth grid of factor scores for each group and latent.
#' For every latent \eqn{k} and group \eqn{g}, it takes the empirical range of
#' factor scores \eqn{\eta_k} (estimated from the original data, but not returned),
#' creates an evenly spaced sequence of length `length.out`, and sets all *other*
#' latents to zero (or to their empirical group means). For each indicator–latent
#' pair \eqn{(j,k)}, the function computes one-latent-at-a-time expectations
#' \deqn{\hat{y}_{jk} = \nu_j + \lambda_{jk}\,\eta_k}
#' and delta-method CIs using `vcov(fit)`. The output contains no original observed
#' indicators — only the synthetic latent grid and model-based curves.
#'
#' In addition, if `se = TRUE`, group-wise summaries of factor-score standard errors
#' are attached: for each *(group × latent)* block, the SE column is filled by a
#' summary (default **median**) of the original factor-score SEs obtained from
#' [lavPredict_parallel(se = TRUE)]. This provides a practical SE reference on the
#' grid, even though the grid rows are synthetic and not tied to real observations.
#'
#' @param fit A converged [`lavaan`] object with a mean structure and measurement
#'   parameters (lambda + nu).
#' @param info Optional list returned by [model_info()]. If `NULL`, it is computed.
#' @param data Optional precomputed factor-score data (a data.frame/tibble)
#'   as returned by `lavPredict_parallel(fit, return_type = "data", se = TRUE)`.
#'   If provided, the function **uses this object directly** to compute ranges and
#'   SE summaries and **does not** call `lavPredict_parallel()` again.
#' @param level Confidence level for delta-method CIs; default `0.95`.
#' @param vcov_type Optional character string forwarded to [lavaan::vcov()] as `type`
#'   (e.g., `"robust.sem"`). If `NULL` (default), `vcov(fit)` is used.
#' @param length.out Integer; grid length per *(group × latent)* block. Default `100`.
#' @param other_latents How to fill non-target latents when varying one latent:
#'   `"zero"` (default) or `"mean"` (group-wise empirical means of factor scores).
#' @param latent_var_as_factor Logical; if `TRUE` (default) `.latent_var` is a factor
#'   with levels ordered as in the model's latent names. If `FALSE`, it remains character.
#' @param se Logical; if `TRUE` (default), attach factor-score SE columns using
#'   a group-wise summary of original SEs.
#' @param se_summary One of `"median"` (default) or `"mean"`. Determines how original
#'   factor-score SEs are summarised within *(group × latent)* to fill the grid.
#'
#' @return A `tibble` with:
#' \itemize{
#'   \item `.rid` — row ID of the synthetic grid,
#'   \item `.gid`, `.group` — numeric ID and label of the lavaan group,
#'   \item `.latent_var` — name of the latent whose sequence varies in the row-block,
#'   \item latent columns (factor-score grid), **immediately followed** by factor-score SE columns
#'         (if `se = TRUE`, named `".se_<latent>"`),
#'   \item per-indicator/per-latent model-based expectations and CIs:
#'     \code{m_est_<indicator>_<latent>}, \code{m_lwr_<indicator>_<latent>},
#'     \code{m_upr_<indicator>_<latent>}.
#' }
#'
#' @section Notes:
#' * Only **continuous** indicators are used; ordinal indicators (if present) are ignored.
#' * Empirical ranges are computed from factor scores via [lavPredict_parallel()], but the
#'   original observed data are not returned.
#' * If a latent's empirical range in a group is degenerate (`min == max`) or non-finite,
#'   a fallback range `[-3, 3]` is used.
#' * CIs use the delta method for the free parameters \eqn{\{\nu_j, \lambda_{jk}\}} based on
#'   the relevant submatrix of `vcov(fit)`.
#'
#' @export
prepare_continuous_old <- function(fit,
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
    require_meanstructure = TRUE,
    require_latent        = TRUE,
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
  ov_cont      <- info$ov_continuous
  latent_names <- info$latent_variables
  n_groups     <- info$n_groups %||% 1L
  gvar         <- info$group_var %||% ".group"
  g_labels     <- info$group_labels %||% as.character(seq_len(n_groups))

  if (length(ov_cont) == 0L) rlang::abort("No continuous indicators in the model.")
  if (length(latent_names) == 0L) rlang::abort("No latent variables detected.")

  # -- Factor scores for ranges (and optional SE summaries) --------------------
  fs_df <- if (is.null(data)) {
    lavPredict_parallel(fit, return_type = "data", se = isTRUE(se), prefix_se = prefix_se)
  } else {
    data
  }

  # Mini-normalization: allow `data` to be a per-group list
  if (is.list(fs_df) && !is.data.frame(fs_df)) {
    fs_list <- fs_df
    fs_df <- dplyr::bind_rows(lapply(seq_along(fs_list), function(i) {
      df_g <- tibble::as_tibble(fs_list[[i]])
      # Prefer matching by *label* to lavaan's g_labels; fallback to i
      this_lab <- if (length(names(fs_list))) names(fs_list)[i] else NA_character_
      gid <- match(this_lab, g_labels)
      if (is.na(gid)) gid <- i
      df_g$.group <- if (!is.null(this_lab) && !is.na(this_lab)) this_lab else g_labels[gid]
      df_g$.gid   <- gid
      df_g
    }))
  }

  # -- Canonical group mapping (FIX): never trust incoming .gid ----------------
  # Always rebuild .gid from the observed *label* using lavaan's group order.
  # 1) Identify a label column to use (prefer gvar; else .group if present).
  if (!(gvar %in% names(fs_df))) {
    if (".group" %in% names(fs_df)) gvar <- ".group" else gvar <- NA_character_
  }
  if (is.na(gvar)) {
    # Single-group fallback with a synthetic label
    fs_df$.group <- g_labels[1L]
  } else {
    fs_df$.group <- as.character(fs_df[[gvar]])
  }
  # 2) Rebuild .gid strictly via match to lavaan's group labels
  fs_df$.gid <- as.integer(match(fs_df$.group, g_labels))

  # 3) Guard: if any NA (label not in lavaan), coerce to first group with warning
  if (anyNA(fs_df$.gid)) {
    bad <- unique(fs_df$.group[is.na(fs_df$.gid)])
    warning(sprintf(
      "Unknown group label(s) in factor-score data: %s. Coercing to the first lavaan group.",
      paste(bad, collapse = ", ")
    ))
    fs_df$.gid[is.na(fs_df$.gid)] <- 1L
    fs_df$.group[is.na(match(fs_df$.group, g_labels))] <- g_labels[1L]
  }

  # 4) For single-group models, enforce (.gid = 1, .group = g_labels[1])
  if (n_groups == 1L) {
    fs_df$.gid   <- 1L
    fs_df$.group <- g_labels[1L]
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
                                     .names = "{.col}"), .groups = "drop")
    # Sanitize non-finite group means: replace with 0 and warn
    for (i in seq_len(nrow(means_by_g))) {
      gid_i <- means_by_g$.gid[i]
      for (ln in latent_names) if (!is.finite(means_by_g[[ln]][i])) {
        warning(sprintf("Replacing non-finite mean for latent '%s' in group %s with 0.", ln, g_labels[gid_i]))
        means_by_g[[ln]][i] <- 0
      }
    }
  }

  # -- (Optional) group-wise summaries of FS SEs -------------------------------
  se_cols_exist <- isTRUE(se) && all(paste0(prefix_se, latent_names) %in% names(fs_df))
  se_summ <- NULL
  if (se_cols_exist) {
    fun <- if (se_summary == "median") stats::median else mean
    se_summ <- fs_df |>
      dplyr::group_by(.gid) |>
      dplyr::summarise(dplyr::across(dplyr::all_of(paste0(prefix_se, latent_names)), ~fun(.x[is.finite(.x)]), .names = "{.col}"), .groups = "drop")
  }

  # -- Build the synthetic grid ------------------------------------------------
  grid_list <- vector("list", n_groups * length(latent_names))
  idx <- 0L
  for (g in seq_len(n_groups)) {
    # Base fill for non-target latents (zeros or group means)
    base_vec <- stats::setNames(numeric(length(latent_names)), latent_names)
    if (!is.null(means_by_g)) {
      mu_row <- means_by_g[means_by_g$.gid == g, , drop = FALSE]
      if (nrow(mu_row) == 1L) for (ln in latent_names) base_vec[[ln]] <- as.numeric(mu_row[[ln]])
    }
    # Ranges per latent in this group
    fs_g <- fs_df[fs_df$.gid == g, latent_names, drop = FALSE]
    for (k in seq_along(latent_names)) {
      kname <- latent_names[k]
      r_gk  <- rng_or_fallback(fs_g[[kname]])
      seq_k <- seq(from = r_gk[1], to = r_gk[2], length.out = length.out)

      # Materialize rows; vary only k-th latent
      mat <- matrix(rep(base_vec, each = length.out), ncol = length(base_vec))
      colnames(mat) <- names(base_vec)
      mat[, kname]  <- seq_k

      tb <- tibble::tibble(
        .gid        = g,
        .group      = g_labels[g],   # <- use lavaan's canonical label directly (no data-driven remap)
        .latent_var = if (isTRUE(latent_var_as_factor)) factor(kname, levels = latent_names) else kname,
        !!!as.data.frame(mat, stringsAsFactors = FALSE)
      )

      if (!is.null(se_summ)) {
        se_row <- se_summ[se_summ$.gid == g, , drop = FALSE]
        if (nrow(se_row) == 1L) for (ln in latent_names) {
          nm <- paste0(prefix_se, ln)
          if (nm %in% names(se_row)) tb[[nm]] <- rep(as.numeric(se_row[[nm]]), length.out)
        }
      }

      idx <- idx + 1L
      grid_list[[idx]] <- tb
    }
  }

  out <- dplyr::bind_rows(grid_list)
  out$.rid <- seq_len(nrow(out))

  # Column order ---------------------------------------------------------------
  out <- out |>
    dplyr::relocate(.rid, .before = 1) |>
    dplyr::relocate(.gid, .after = .rid) |>
    dplyr::relocate(.group, .after = .gid) |>
    dplyr::relocate(.latent_var, .after = .group) |>
    dplyr::relocate(dplyr::all_of(latent_names), .after = ".latent_var")

  if (!is.null(se_summ)) {
    se_cols <- paste0(prefix_se, latent_names)
    se_cols <- se_cols[se_cols %in% names(out)]
    if (length(se_cols)) out <- dplyr::relocate(out, dplyr::all_of(se_cols), .after = dplyr::all_of(latent_names[length(latent_names)]))
  }

  # -- Measurement pieces & vcov for delta method ------------------------------
  .build_lambda_nu_from_par <- function(par_tbl, info) {
    G  <- info$n_groups
    ov <- info$observed_variables
    lv <- info$latent_variables
    Lambda <- vector("list", G)
    Nu     <- vector("list", G)
    for (g in seq_len(G)) {
      Lg  <- matrix(0, nrow = length(ov), ncol = length(lv), dimnames = list(ov, lv))
      Nug <- stats::setNames(numeric(length(ov)), ov)
      wl  <- which(par_tbl$group == g & par_tbl$op == "=~" & par_tbl$lhs %in% lv & par_tbl$rhs %in% ov)
      if (length(wl)) for (i in wl) Lg[par_tbl$rhs[i], par_tbl$lhs[i]] <- par_tbl$est[i]
      wn  <- which(par_tbl$group == g & par_tbl$op == "~1" & par_tbl$lhs %in% ov)
      if (length(wn)) for (i in wn) Nug[par_tbl$lhs[i]] <- par_tbl$est[i]
      Lambda[[g]] <- Lg
      Nu[[g]]     <- Nug
    }
    list(lambda = Lambda, nu = Nu)
  }

  V <- tryCatch({ if (is.null(vcov_type)) lavaan::vcov(fit) else lavaan::vcov(fit, type = vcov_type) }, error = function(e) NULL)
  if (is.null(V)) rlang::abort("vcov(fit) is not available; cannot compute delta-method CIs.")
  PT    <- lavaan::parTable(fit)
  meas  <- .build_lambda_nu_from_par(PT, info)
  zcrit <- stats::qnorm(1 - (1 - level)/2)
  single_grp <- (n_groups == 1L)

  get_free_lambda <- function(g, j, k) {
    if (single_grp) w <- which(PT$op == "=~" & PT$lhs == k & PT$rhs == j) else w <- which(PT$group == g & PT$op == "=~" & PT$lhs == k & PT$rhs == j)
    if (!length(w)) return(0L) else PT$free[w[1L]]
  }
  get_free_nu <- function(g, j) {
    if (single_grp) w <- which(PT$op == "~1" & PT$lhs == j) else w <- which(PT$group == g & PT$op == "~1" & PT$lhs == j)
    if (!length(w)) return(0L) else PT$free[w[1L]]
  }

  out2 <- out
  for (g in seq_len(n_groups)) {
    rows_g <- which(out2$.gid == g)
    if (!length(rows_g)) next

    lam_g <- meas$lambda[[g]]
    nu_g  <- meas$nu[[g]]

    ov_all <- info$observed_variables
    if (!is.null(rownames(lam_g))) {
      idx_rows <- match(ov_all, rownames(lam_g), nomatch = 0L)
      lam_g    <- lam_g[idx_rows[idx_rows > 0L], , drop = FALSE]
      lam_g    <- lam_g[match(ov_cont, rownames(lam_g)), , drop = FALSE]
    } else {
      lam_g <- lam_g[match(ov_cont, ov_all), , drop = FALSE]
    }

    if (!is.null(colnames(lam_g))) lam_g <- lam_g[, latent_names, drop = FALSE]
    colnames(lam_g) <- latent_names

    if (is.null(nu_g) || length(nu_g) == 0L) {
      nu_g <- stats::setNames(rep(0, length(ov_cont)), ov_cont)
    } else if (!is.null(names(nu_g))) {
      nu_g <- nu_g[ov_cont]
    } else {
      nu_g <- nu_g[match(ov_cont, ov_all)]
    }

    for (j_idx in seq_along(ov_cont)) {
      j <- ov_cont[j_idx]
      lam_row <- lam_g[j_idx, , drop = TRUE]
      nz_idx  <- which(is.finite(lam_row) & lam_row != 0)
      if (!length(nz_idx)) next

      if (identical(other_latents, "mean")) {
        mu_row <- means_by_g[means_by_g$.gid == g, , drop = FALSE]
        if (nrow(mu_row) == 1L) {
          base_eta <- as.numeric(mu_row[, latent_names, drop = TRUE])
          names(base_eta) <- latent_names
        } else base_eta <- stats::setNames(rep(0, length(latent_names)), latent_names)
      } else base_eta <- stats::setNames(rep(0, length(latent_names)), latent_names)

      mu_const_all <- sum(lam_row[nz_idx] * base_eta[nz_idx])

      for (k_idx in nz_idx) {
        k_name <- latent_names[k_idx]
        rows_k <- which(out2$.gid == g & out2$.latent_var == k_name)
        eta_k  <- out2[[k_name]][rows_k]
        mu_vec <- as.numeric(nu_g[j_idx] + lam_row[k_idx] * eta_k + (mu_const_all - lam_row[k_idx] * base_eta[k_name]))

        id_nu  <- get_free_nu(g, j)
        id_lam <- get_free_lambda(g, j, k_name)

        extra_ids <- integer(0)
        extra_A   <- NULL
        if (identical(other_latents, "mean")) {
          other_r <- setdiff(nz_idx, k_idx)
          if (length(other_r)) {
            ids_tmp <- vapply(other_r, function(r) get_free_lambda(g, j, latent_names[r]), integer(1))
            keep    <- which(ids_tmp > 0L)
            if (length(keep)) {
              other_r  <- other_r[keep]
              extra_ids <- ids_tmp[keep]
              extra_A   <- do.call(cbind, lapply(other_r, function(r) rep(base_eta[r], length(rows_k))))
            }
          }
        }

        cols <- c(if (id_nu > 0L) id_nu else integer(0), if (id_lam > 0L) id_lam else integer(0), extra_ids)

        if (length(cols)) {
          uniq <- unique(cols)
          A_core <- cbind(if (id_nu  > 0L) rep(1, length(rows_k)) else NULL, if (id_lam > 0L) eta_k else NULL)
          if (!is.null(extra_A)) A_core <- cbind(A_core, extra_A)
          pos <- match(uniq, cols)
          A   <- A_core[, pos, drop = FALSE]
          Vsub <- V[uniq, uniq, drop = FALSE]
          se   <- sqrt(pmax(rowSums((A %*% Vsub) * A), 0))
        } else se <- rep(0, length(rows_k))

        base  <- paste0(j, "_", k_name)
        out2[rows_k, paste0("m_est_", base)] <- mu_vec
        out2[rows_k, paste0("m_lwr_", base)] <- mu_vec - zcrit * se
        out2[rows_k, paste0("m_upr_", base)] <- mu_vec + zcrit * se
      }
    }
  }

  is_num   <- vapply(out2, is.numeric, logical(1L))
  keep_int <- names(out2) %in% c(".gid", ".rid")
  to_dbl   <- is_num & !keep_int
  for (nm in names(out2)[to_dbl]) out2[[nm]] <- as.double(out2[[nm]])

  print(info$group_labels)

  print(unique(fs_df$gender))

  tibble::as_tibble(out2)
}
