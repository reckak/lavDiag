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

  info     <- model_info(fit)
  ov_all   <- info$observed_variables
  ov_ord   <- info$ov_ordinal
  ov_cont  <- info$ov_continuous
  is_mg    <- isTRUE(info$n_groups > 1L)
  gvar     <- info$group_var %||% ".group"
  has_ord  <- length(ov_ord)  > 0L
  has_cont <- length(ov_cont) > 0L

  if (is.null(workers)) workers <- max(1L, parallel::detectCores(logical = TRUE) - 1L)

  # Extract fitted data once (preserve original classes)
  if (!is_mg) {
    dat_original <- tibble::as_tibble(lavaan::lavInspect(fit, "data"))
  } else {
    dat_original <- lavaan::lavInspect(fit, "data") |>
      purrr::map(tibble::as_tibble) |>
      dplyr::bind_rows(.id = gvar)
  }

  ov_ord <- ov_ord[ov_ord %in% names(dat_original)]

  # -------- Fast serial path for continuous-only or single-worker ------------
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
      se_attr  <- tryCatch(attr(base, "se"), error = function(e) NULL)
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
    warning("`se=TRUE` is only supported for continuous-only models; ignoring for ordinal/mixed models.")
  }

  # -------- Helpers for dummy rows -------------------------------------------

  create_allcat_dummy <- function(df, ord_vars, group_var = NULL) {
    create_dummy_from_original(df, ord_vars, group_var = group_var)
  }

  # Minimal-cover variance dummy (2 prototype rows per group)
  create_continuous_variance_dummy_minimal <- function(df, cont_vars, group_var = NULL) {
    if (!length(cont_vars) || nrow(df) == 0L) return(df[0, , drop = FALSE])

    build_for_group <- function(dg) {
      if (nrow(dg) == 0L) return(dg[0, , drop = FALSE])
      proto1 <- dg[1, , drop = FALSE]
      proto2 <- dg[min(2L, nrow(dg)), , drop = FALSE]
      any_assigned <- FALSE

      for (i in seq_along(cont_vars)) {
        v <- cont_vars[i]; if (!v %in% names(dg)) next
        x <- dg[[v]]
        ux <- unique(x[is.finite(suppressWarnings(as.numeric(x))) & !is.na(x)])
        if (length(ux) < 2L) next
        if (i %% 2L == 1L) {
          proto1[[v]] <- .cast_like(ux[1], x)
          proto2[[v]] <- .cast_like(ux[2], x)
        } else {
          proto1[[v]] <- .cast_like(ux[2], x)
          proto2[[v]] <- .cast_like(ux[1], x)
        }
        any_assigned <- TRUE
      }

      if (!any_assigned) return(dg[0, , drop = FALSE])
      dplyr::bind_rows(proto1, proto2) |> dplyr::distinct()
    }

    if (is.null(group_var)) {
      build_for_group(df)
    } else {
      df |>
        dplyr::group_split(.data[[group_var]], .keep = TRUE) |>
        lapply(build_for_group) |>
        dplyr::bind_rows()
    }
  }

  # -------- Driver data & index ----------------------------------------------
  use_distinct <- switch(distinct,
                         "always" = TRUE,
                         "never"  = FALSE,
                         "auto"   = (nrow(dat_original) >= distinct_threshold))
  if (has_cont) use_distinct <- FALSE

  dat_driver <- if (use_distinct) dplyr::distinct(dat_original) else dat_original

  # Deterministic row index (NOTE: lavaan won't keep it; we'll re-attach it after prediction)
  dat_driver$.ix <- seq_len(nrow(dat_driver))

  # -------- Dummy rows --------------------------------------------------------
  dummy_ord  <- create_allcat_dummy(dat_original, ov_ord, group_var = if (is_mg) gvar else NULL)
  dummy_cont <- create_continuous_variance_dummy_minimal(dat_original, ov_cont, group_var = if (is_mg) gvar else NULL)

  dummy <- dplyr::bind_rows(dummy_ord, dummy_cont) |> dplyr::distinct()

  # Make dummy columns match driver (including .ix) and mark them as dummy via NA .ix
  missing_cols <- setdiff(names(dat_driver), names(dummy))
  if (length(missing_cols)) for (cc in missing_cols) dummy[[cc]] <- NA
  dummy <- dummy[, names(dat_driver), drop = FALSE]
  dummy$.ix <- NA_integer_

  # -------- Chunking ----------------------------------------------------------
  n_rows <- nrow(dat_driver)
  if (n_rows == 0L) return(dat_original)
  if (is.null(chunk_size)) {
    n_chunks   <- max(1L, workers)
    chunk_size <- ceiling(n_rows / n_chunks)
  }
  idx <- split(seq_len(n_rows), ceiling(seq_along(seq_len(n_rows)) / chunk_size))
  chunked <- lapply(idx, function(ix) dat_driver[ix, , drop = FALSE])

  # Prepend dummy to every chunk
  if (nrow(dummy) > 0L) {
    chunked <- lapply(chunked, function(x) dplyr::bind_rows(dummy, x))
  }

  # -------- Future plan -------------------------------------------------------
  reset_plan <- .set_future_plan(plan = plan, workers = workers)
  on.exit(reset_plan(), add = TRUE)

  # -------- Predict -----------------------------------------------------------
  fopts <- list(append.data = TRUE, assemble = TRUE)
  if (is_mg) fopts$drop.list.single.group <- FALSE
  dots <- rlang::dots_list(...)

  pred_fun <- function(df_chunk, fit, fopts, dots) {
    # Run lavaan prediction on the chunk
    pred <- do.call(lavaan::lavPredict, c(list(object = fit, newdata = df_chunk), fopts, dots))
    pred <- tibble::as_tibble(pred)
    # IMPORTANT: Re-attach the driver index lost by lavaan (same row order/length)
    # This allows us to drop dummy rows (NA .ix) and align back deterministically.
    pred$.ix <- df_chunk$.ix
    pred
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

  # -------- Drop dummy rows & bind -------------------------------------------
  drop_dummy <- function(df) df[!is.na(df$.ix), , drop = FALSE]
  out_list <- lapply(out_list, drop_dummy)

  out_all <- do.call(rbind, out_list) |> tibble::as_tibble()

  # -------- Deterministic merge-back -----------------------------------------
  if (use_distinct) {
    # Distinct path: safe left_join on all shared columns
    join_keys <- intersect(names(dat_driver), names(out_all))
    if (!length(join_keys)) stop("lavPredict_parallel: no common join keys.", call. = FALSE)
    out_joined <- dplyr::left_join(dat_original, out_all, by = join_keys)
  } else {
    # Order predictions by .ix and align to dat_driver
    ord_ix <- order(out_all$.ix)
    out_all <- out_all[ord_ix, , drop = FALSE]

    if (nrow(out_all) != nrow(dat_driver)) {
      stop("lavPredict_parallel: row count mismatch after dropping dummy rows.", call. = FALSE)
    }
    if (!identical(out_all$.ix, dat_driver$.ix)) {
      m <- match(dat_driver$.ix, out_all$.ix)
      if (anyNA(m)) stop("lavPredict_parallel: cannot align predictions to driver rows.", call. = FALSE)
      out_all <- out_all[m, , drop = FALSE]
    }

    pred_cols <- setdiff(names(out_all), names(dat_driver))
    pred_mat  <- out_all[, pred_cols, drop = FALSE]

    same_shape <- identical(nrow(dat_driver), nrow(dat_original)) &&
      identical(setdiff(names(dat_driver), ".ix"), names(dat_original))
    if (same_shape) {
      out_joined <- dplyr::bind_cols(dat_original, pred_mat)
    } else {
      base <- dat_driver; base$.ix <- NULL
      out_all2 <- dplyr::bind_cols(base, pred_mat)
      out_joined <- dplyr::left_join(dat_original, out_all2, by = intersect(names(dat_original), names(base)))
    }
  }

  # -------- Return ------------------------------------------------------------
  if (!is_mg) return(out_joined)
  if (return_type == "data") return(out_joined)
  grp <- out_joined[[gvar]]; out_joined[[gvar]] <- NULL
  split(out_joined, grp)
}
