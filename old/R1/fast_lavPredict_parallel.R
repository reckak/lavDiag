#' Fast & simple parallel wrapper for lavaan::lavPredict (ordinal-aware)
#'
#' - Uses "all-categories" dummy strategy to stabilize chunks for ordinal models.
#' - Optionally deduplicates input rows (distinct) only for large data.
#'
#' @param fit lavaan/blavaan object.
#' @param workers integer; number of workers (default: max(1, cores-1)).
#' @param plan one of c("auto","multisession","multicore","sequential").
#' @param chunk_size optional integer rows per chunk (default computed).
#' @param return_type One of "list" or "data".
#' @param progress logical; show furrr progress bar.
#' @param se logical; if TRUE and model is continuous-only, attach SE columns.
#' @param prefix_se_fs character; prefix for SE columns (default ".se_").
#' @param distinct one of c("auto","always","never"); use distinct() to shrink work.
#' @param distinct_threshold numeric; when distinct="auto", apply distinct() only if
#'   nrow(data) >= distinct_threshold (default 1e4).
#' @param ... extra args passed to lavaan::lavPredict().
#' @return tibble (single-group) or list/tibble (multi-group; see return_type).
#' Fast & simple parallel wrapper for lavaan::lavPredict (ordinal-aware)
#'
#' - Uses "all-categories" dummy strategy to stabilize chunks for ordinal models.
#' - Optionally deduplicates input rows (distinct) only for large data.
#'
#' @param fit lavaan/blavaan object.
#' @param workers integer; number of workers (default: max(1, cores-1)).
#' @param plan one of c("auto","multisession","multicore","sequential").
#' @param chunk_size optional integer rows per chunk (default computed).
#' @param return_type One of "list" or "data".
#' @param progress logical; show furrr progress bar.
#' @param se logical; if TRUE and model is continuous-only, attach SE columns.
#' @param prefix_se_fs character; prefix for SE columns (default ".se_").
#' @param distinct one of c("never","auto","always"); use distinct() to shrink work.
#' @param distinct_threshold numeric; when distinct="auto", apply distinct() only if
#'   nrow(data) >= distinct_threshold (default 5e4).
#' @param ... extra args passed to lavaan::lavPredict().
#' @return tibble (single-group) or list/tibble (multi-group; see return_type).
lavPredict_parallel <- function(fit,
                                workers            = NULL,
                                plan               = c("auto","multisession","multicore","sequential"),
                                chunk_size         = NULL,
                                return_type        = c("list","data"),
                                progress           = FALSE,
                                se                 = FALSE,
                                prefix_se_fs       = ".se_",
                                distinct           = c("never","auto","always"),
                                distinct_threshold = 5e4,
                                ...) {
  # -- Checks & meta -----------------------------------------------------------
  .assert_lavaan_fit(fit)
  plan        <- match.arg(plan)
  return_type <- match.arg(return_type)
  distinct    <- match.arg(distinct)

  info    <- model_info(fit)
  ov_ord  <- info$ov_ordinal
  is_mg   <- isTRUE(info$n_groups > 1L)
  gvar    <- info$group_var %||% ".group"

  if (is.null(workers)) workers <- max(1L, parallel::detectCores(logical = TRUE) - 1L)

  # Extract fitted data once (preserve original classes)
  if (!is_mg) {
    dat_original <- tibble::as_tibble(lavaan::lavInspect(fit, "data"))
  } else {
    dat_original <- lavaan::lavInspect(fit, "data") |>
      purrr::map(tibble::as_tibble) |>
      dplyr::bind_rows(.id = gvar)
  }

  # Keep truly present ordinal vars
  ov_ord  <- ov_ord[ov_ord %in% names(dat_original)]
  has_ord <- length(ov_ord) > 0L

  # -------- Fast serial path for continuous-only ------------------------------
  if (!has_ord || workers <= 1L) {
    args_base <- list(
      object      = fit,
      newdata     = dat_original,
      append.data = TRUE,
      assemble    = TRUE
    )
    if (is_mg) args_base$drop.list.single.group <- FALSE
    args_base <- c(args_base, rlang::dots_list(...))
    base <- do.call(lavaan::lavPredict, args_base)
    out  <- tibble::as_tibble(base)

    if (isTRUE(se)) {
      lv_names <- tryCatch(lavaan::lavNames(fit, type = "lv"), error = function(e) character())
      se_attr <- tryCatch(attr(base, "se"), error = function(e) NULL)
      if (is.null(se_attr)) {
        args_se <- list(object = fit, newdata = dat_original, append.data = FALSE, assemble = TRUE, se = TRUE)
        if (is_mg) args_se$drop.list.single.group <- FALSE
        args_se <- c(args_se, rlang::dots_list(...))
        base_se <- tryCatch(do.call(lavaan::lavPredict, args_se), error = function(e) NULL)
        se_attr <- if (!is.null(base_se)) attr(base_se, "se") else NULL
      }
      if (!is.null(se_attr)) {
        if (is.list(se_attr)) se_attr <- do.call(rbind, se_attr)
        se_mat <- as.matrix(se_attr)
        lv_pred_cols <- intersect(lv_names, colnames(out))
        if (length(lv_pred_cols)) {
          if (is.null(colnames(se_mat))) {
            if (ncol(se_mat) == length(lv_pred_cols)) colnames(se_mat) <- lv_pred_cols
            else colnames(se_mat) <- lv_pred_cols[seq_len(ncol(se_mat))]
          }
          n_out <- nrow(out); n_se <- nrow(se_mat)
          if (!identical(n_se, n_out)) {
            if (isTRUE(is_mg) && gvar %in% names(out)) {
              grp_fac <- factor(out[[gvar]], levels = unique(out[[gvar]]))
              cnt <- as.integer(table(grp_fac)); G <- length(cnt)
              if (n_se == G) {
                idx_rep <- rep(seq_len(G), times = cnt)
                se_mat  <- se_mat[idx_rep, , drop = FALSE]
              } else if (n_se == 1L) {
                se_mat <- se_mat[rep(1L, n_out), , drop = FALSE]
              } else {
                warning("SE matrix has unexpected number of rows (", n_se, "); skipping SE attachment.")
                se_mat <- NULL
              }
            } else if (n_se == 1L) {
              se_mat <- se_mat[rep(1L, n_out), , drop = FALSE]
            } else {
              warning("SE matrix has unexpected number of rows (", n_se, "); skipping SE attachment.")
              se_mat <- NULL
            }
          }
          if (!is.null(se_mat)) {
            keep <- intersect(colnames(se_mat), lv_pred_cols)
            if (length(keep)) {
              se_df <- as.data.frame(se_mat[, keep, drop = FALSE])
              names(se_df) <- paste0(prefix_se_fs, names(se_df))
              out <- dplyr::bind_cols(out, se_df)
            }
          }
        }
      } else if (has_ord) {
        warning("SE not available from lavaan for ordinal models (ignored).")
      }
    }

    if (!is_mg) return(out)
    if (return_type == "data") return(out)
    grp <- out[[gvar]]; out[[gvar]] <- NULL
    return(split(out, grp))
  }

  if (isTRUE(se)) {
    warning("`se=TRUE` is only supported for continuous-only models; ignoring for ordinal models.")
  }

  # -------- ORDINAL PARALLEL PATH (ALL-categories dummy) ----------------------

  key_cols <- c(if (is_mg) gvar else NULL, ov_ord)

  # decide distinct
  use_distinct <- switch(distinct,
                         "always" = TRUE,
                         "never"  = FALSE,
                         "auto"   = (nrow(dat_original) >= distinct_threshold))

  # Build ALL-categories dummy per ordinal var (per group if MG) — once
  build_all_dummy <- function(df_all, ord_vars, group_var = NULL) {
    if (!length(ord_vars) || nrow(df_all) == 0L) return(df_all[0, , drop = FALSE])

    if (is.null(group_var)) {
      proto <- df_all[1, , drop = FALSE]
      out_list <- vector("list", length(ord_vars))
      for (i in seq_along(ord_vars)) {
        v <- ord_vars[i]
        x_all <- df_all[[v]]
        cats  <- if (is.factor(x_all)) levels(x_all) else sort(unique(x_all))
        if (!length(cats)) { out_list[[i]] <- df_all[0, , drop = FALSE]; next }
        rows <- lapply(cats, function(cat) { r <- proto; r[[v]] <- .cast_like(cat, x_all); r })
        out_list[[i]] <- dplyr::bind_rows(rows)
      }
      dplyr::bind_rows(out_list) |> dplyr::distinct()
    } else {
      split_g <- split(df_all, df_all[[group_var]], drop = TRUE)
      out_g <- lapply(split_g, function(dg) {
        if (nrow(dg) == 0L) return(dg[0, , drop = FALSE])
        proto <- dg[1, , drop = FALSE]
        parts <- vector("list", length(ord_vars))
        gname <- dg[[group_var]][1]
        for (i in seq_along(ord_vars)) {
          v <- ord_vars[i]
          x_all <- df_all[df_all[[group_var]] == gname, v, drop = TRUE]
          cats  <- if (is.factor(x_all)) levels(x_all) else sort(unique(x_all))
          if (!length(cats)) { parts[[i]] <- dg[0, , drop = FALSE]; next }
          rows <- lapply(cats, function(cat) { r <- proto; r[[v]] <- .cast_like(cat, x_all); r })
          parts[[i]] <- dplyr::bind_rows(rows)
        }
        dplyr::bind_rows(parts) |> dplyr::distinct()
      })
      dplyr::bind_rows(out_g)
    }
  }

  # DRIVER data
  if (use_distinct) {
    dat_driver <- dplyr::distinct(dat_original[, key_cols, drop = FALSE])
  } else {
    dat_driver <- dat_original[, key_cols, drop = FALSE]
  }

  # Dummy dle "ALL"
  dummy <- build_all_dummy(dat_original[, key_cols, drop = FALSE],
                           ov_ord,
                           group_var = if (is_mg) gvar else NULL)
  dummy <- dummy[, key_cols, drop = FALSE]

  # Chunking: if not provided, make ≈ workers chunks (low overhead)
  n_rows <- nrow(dat_driver)
  if (n_rows == 0L) return(dat_original)
  if (is.null(chunk_size)) {
    n_chunks   <- max(1L, workers)
    chunk_size <- ceiling(n_rows / n_chunks)
  }
  idx <- split(seq_len(n_rows), ceiling(seq_along(seq_len(n_rows)) / chunk_size))
  chunked <- lapply(idx, function(ix) dat_driver[ix, , drop = FALSE])

  # Prepend dummy (once per chunk)
  if (nrow(dummy) > 0L) {
    chunked <- lapply(chunked, function(x) dplyr::bind_rows(dummy, x))
  }

  # Plan set/reset
  reset_plan <- .set_future_plan(plan = plan, workers = workers)
  on.exit(reset_plan(), add = TRUE)

  # Predict
  fopts <- list(append.data = TRUE, assemble = TRUE)
  if (is_mg) fopts$drop.list.single.group <- FALSE
  dots <- rlang::dots_list(...)

  pred_fun <- function(df_chunk, fit, fopts, dots) {
    do.call(lavaan::lavPredict, c(list(object = fit, newdata = df_chunk), fopts, dots))
  }

  out_list <- base::suppressWarnings(
    furrr::future_map(
      chunked,
      pred_fun,
      fit   = fit,
      fopts = fopts,
      dots  = dots,
      .progress = progress,
      .options  = furrr::furrr_options(
        seed     = TRUE,
        packages = c("lavaan","stats","dplyr","tibble"),
        globals  = FALSE
      )
    )
  )

  # Drop dummy rows
  dN <- nrow(dummy)
  for (i in seq_along(out_list)) {
    if (dN > 0L) out_list[[i]] <- out_list[[i]][-(seq_len(dN)), , drop = FALSE]
  }
  out_all <- do.call(rbind, out_list) |> tibble::as_tibble()
  pred_cols <- setdiff(names(out_all), key_cols)

  if (use_distinct) {
    # Distinct path: map via fast left_join() on keys
    out_joined <- dplyr::left_join(
      dat_original,
      out_all,
      by = key_cols
    )
  } else {
    # Fast path: rely on preserved row order (chunks then rbind)
    if (nrow(out_all) != nrow(dat_original)) {
      stop("lavPredict_parallel: row count mismatch without distinct().")
    }
    pred_mat   <- out_all[, pred_cols, drop = FALSE]
    out_joined <- dplyr::bind_cols(dat_original, pred_mat)
  }

  if (!is_mg) return(out_joined)
  if (return_type == "data") return(out_joined)
  grp <- out_joined[[gvar]]; out_joined[[gvar]] <- NULL
  split(out_joined, grp)
}

