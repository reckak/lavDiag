#' Fast & simple parallel wrapper for lavaan::lavPredict (ordinal-only)
#'
#' Parallelization is used only when the model has ordinal indicators (ov.ord)
#' and workers > 1. For continuous-only models it falls back to a single
#' serial call (lavPredict is fast enough there).
#'
#' Optionally, for continuous models only, attaches standard errors of factor scores
#' using lavaan's built-in support (no bootstrap).
#'
#' @param fit lavaan/blavaan object.
#' @param workers integer, number of workers (default: max(1, cores-1)).
#' @param plan one of c("auto","multisession","multicore","sequential").
#' @param chunk_size optional integer rows per chunk (default computed as ceiling(n_unique/workers)).
#' @param return_type One of "list" or "data".
#' @param progress logical; show furrr progress bar.
#' @param se logical; if TRUE and model is continuous-only, attach SE columns.
#' @param prefix_se_fs character; prefix for SE columns (default ".se_"), e.g., ".se_" -> ".se_<factor>".
#' @param ... extra args passed to lavaan::lavPredict().
#' @return tibble (single-group) or list/tibble (multi-group; see return_type).
lavPredict_parallel <- function(fit,
                                workers        = NULL,
                                plan           = c("auto","multisession","multicore","sequential"),
                                chunk_size     = NULL,
                                return_type    = c("list","data"),
                                progress       = FALSE,
                                se             = FALSE,
                                prefix_se_fs   = ".se_",
                                ...) {
  # -- Strict assertions & meta ------------------------------------------------
  .assert_lavaan_fit(fit)
  plan        <- match.arg(plan)
  return_type <- match.arg(return_type)

  info    <- model_info(fit)
  ov_ord  <- info$ov_ordinal
  is_mg   <- isTRUE(info$n_groups > 1L)
  gvar    <- info$group_var %||% ".group"

  # Extract fitted data once (preserve original classes)
  if (!is_mg) {
    dat_original <- lavaan::lavInspect(fit, "data") %>% tibble::as_tibble()
  } else {
    dat_original <- lavaan::lavInspect(fit, "data") %>%
      purrr::map(tibble::as_tibble) %>%
      dplyr::bind_rows(.id = gvar)
  }

  # Keep ov_ord that truly exist in data (no retyping!)
  ov_ord  <- ov_ord[ov_ord %in% names(dat_original)]
  has_ord <- length(ov_ord) > 0L

  # Decide parallelization: ordinal-only
  if (is.null(workers)) workers <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
  do_parallel <- has_ord && workers > 1L

  # Latent variable names to detect score columns
  lv_names <- tryCatch(lavaan::lavNames(fit, type = "lv"), error = function(e) character())

  # -- Fast serial path for continuous-only -----------------------------------
  if (!do_parallel) {
    # Build args without passing NULL for drop.list.single.group
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

    # Attach SE only if requested (continuous-only)
    if (isTRUE(se)) {
      # 1) Try attribute "se" directly
      se_attr <- tryCatch(attr(base, "se"), error = function(e) NULL)

      # 2) If not present, try a second call with se=TRUE (some lavaan versions need this)
      if (is.null(se_attr)) {
        args_se <- list(
          object      = fit,
          newdata     = dat_original,
          append.data = FALSE,
          assemble    = TRUE,
          se          = TRUE
        )
        if (is_mg) args_se$drop.list.single.group <- FALSE
        args_se <- c(args_se, rlang::dots_list(...))

        base_se <- tryCatch(do.call(lavaan::lavPredict, args_se), error = function(e) NULL)
        se_attr <- if (!is.null(base_se)) attr(base_se, "se") else NULL
      }

      # 3) Bind SE columns if obtained  ---- robust for MG + pairwise missing
      if (!is.null(se_attr)) {
        # Normalize to matrix
        if (is.list(se_attr)) se_attr <- do.call(rbind, se_attr)
        se_mat <- as.matrix(se_attr)

        # Score columns present in 'out'
        lv_pred_cols <- intersect(lv_names, colnames(out))
        if (length(lv_pred_cols)) {
          # Ensure column names (fallback to lv_pred_cols order if missing)
          if (is.null(colnames(se_mat))) {
            if (ncol(se_mat) == length(lv_pred_cols)) {
              colnames(se_mat) <- lv_pred_cols
            } else {
              colnames(se_mat) <- lv_pred_cols[seq_len(ncol(se_mat))]
            }
          }

          # --- shape alignment ---
          n_out <- nrow(out)
          n_se  <- nrow(se_mat)

          if (!identical(n_se, n_out)) {
            if (isTRUE(is_mg) && gvar %in% names(out)) {
              # counts per group in 'out'
              grp_fac <- factor(out[[gvar]], levels = unique(out[[gvar]]))
              cnt <- as.integer(table(grp_fac))
              G   <- length(cnt)

              if (n_se == G) {
                # Repeat each group's SE row 'cnt[g]' times
                idx_rep <- rep(seq_len(G), times = cnt)
                se_mat  <- se_mat[idx_rep, , drop = FALSE]
              } else if (n_se == 1L) {
                # Global 1×F SE -> repeat for all rows
                se_mat <- se_mat[rep(1L, n_out), , drop = FALSE]
              } else {
                warning("SE matrix has unexpected number of rows (", n_se,
                        "); skipping SE attachment.")
                se_mat <- NULL
              }
            } else if (n_se == 1L) {
              # Single-group global 1×F -> repeat
              se_mat <- se_mat[rep(1L, n_out), , drop = FALSE]
            } else {
              warning("SE matrix has unexpected number of rows (", n_se,
                      "); skipping SE attachment.")
              se_mat <- NULL
            }
          }

          if (!is.null(se_mat)) {
            # Keep only columns also present in the scores
            keep <- intersect(colnames(se_mat), lv_pred_cols)
            if (length(keep)) {
              se_df <- as.data.frame(se_mat[, keep, drop = FALSE])
              names(se_df) <- paste0(prefix_se_fs, names(se_df))
              # Safe bind (now N rows)
              out <- dplyr::bind_cols(out, se_df)
            }
          }
        }
      } else if (has_ord) {
        # Shouldn't happen in continuous path; just in case
        warning("SE not available from lavaan for ordinal models (ignored).")
      }   # <---- MISSING BRACE WAS HERE (closes: if (isTRUE(se)) )
    }

      if (!is_mg) return(out)
      if (return_type == "data") return(out)
      grp <- out[[gvar]]; out[[gvar]] <- NULL
      return(split(out, grp))
    }

  # -- Parallel path for ordinal models ---------------------------------------
  if (isTRUE(se)) {
    warning("`se=TRUE` is only supported for continuous-only models; ignoring for ordinal models.")
  }

  # 1) unique rows to avoid redundant predictions
  dat_unique <- dplyr::distinct(dat_original)

  # 2) dummy rows = real rows from original data that cover all categories
  #    (no type changes, no factors enforced)
  dummy <- create_dummy_from_original(dat_original, ov_ord, group_var = if (is_mg) gvar else NULL)

  # 3) chunking
  n_rows <- nrow(dat_unique)
  if (n_rows == 0L) return(dat_original)
  if (is.null(chunk_size)) chunk_size <- ceiling(n_rows / workers)
  idx_list <- split(seq_len(n_rows),
                    ceiling(seq_along(seq_len(n_rows)) / chunk_size))
  chunked <- lapply(idx_list, function(ix) dat_unique[ix, , drop = FALSE])

  # prepend dummy to each chunk
  chunked <- purrr::map(chunked, ~ dplyr::bind_rows(dummy, .x))

  # 4) plan set/reset
  reset_plan <- .set_future_plan(plan = plan, workers = workers)
  on.exit(reset_plan(), add = TRUE)

  # 5) predict in parallel
  fopts <- list(append.data = TRUE, assemble = TRUE)
  if (is_mg) fopts$drop.list.single.group <- FALSE

  # Capture dots once to avoid capturing the whole calling environment on workers
  dots <- rlang::dots_list(...)

  # Stateless worker
  pred_fun <- function(df_chunk, fit, fopts, dots) {
    do.call(lavaan::lavPredict,
            c(list(object = fit, newdata = df_chunk), fopts, dots))
  }

  out_list <- base::suppressWarnings(
    furrr::future_map(
      chunked,
      pred_fun,
      fit   = fit,          # pass explicitly => no globals lookup
      fopts = fopts,
      dots  = dots,
      .progress = progress,
      .options  = furrr::furrr_options(
        seed     = TRUE,
        packages = c("lavaan","stats","dplyr","tibble"),
        globals  = FALSE     # <- disable scanning/serializing outer env
      )
    )
  )

  # 6) drop dummy rows
  dN <- nrow(dummy)
  for (i in seq_along(out_list)) {
    out_list[[i]] <- out_list[[i]][-(seq_len(dN)), , drop = FALSE]
  }

  # 7) bind & merge back to original (intersect join keys to avoid missing-key errors)
  out_all   <- do.call(rbind, out_list) %>% tibble::as_tibble()
  join_keys <- intersect(colnames(dat_unique), colnames(out_all))
  if (!length(join_keys)) {
    stop("lavPredict_parallel: no common join keys between data and predictions.", call. = FALSE)
  }
  out_joined <- dplyr::left_join(dat_original, out_all, by = join_keys)

  if (!is_mg) return(out_joined)
  if (return_type == "data") return(out_joined)
  grp <- out_joined[[gvar]]; out_joined[[gvar]] <- NULL
  split(out_joined, grp)
}

# ----- Helpers ---------------------------------------------------------------

# Build a COMPLETE set of dummy rows so that every ordinal variable
# has ALL its categories represented in each chunk.
# - Single-group: one template row taken from data; for each ov.ord variable,
#   create one row per category (using levels() if factor, otherwise unique()).
# - Multi-group: do the same **per group** and ensure at least one row per group.
create_dummy_from_original <- function(dat_original, ov_ord, group_var = NULL) {
  # Nothing to do
  if (length(ov_ord) == 0L) return(dat_original[0, , drop = FALSE])

  ov_ord <- ov_ord[ov_ord %in% names(dat_original)]
  if (length(ov_ord) == 0L) return(dat_original[0, , drop = FALSE])

  # Helper: for a given data frame and one ordinal variable name,
  # produce rows covering ALL categories for that variable.
  build_rows_for_var <- function(df, var) {
    proto <- df[1, , drop = FALSE]  # template row with correct classes
    x     <- df[[var]]

    # Categories: prefer levels() if factor; otherwise distinct observed values
    cats <- if (is.factor(x)) levels(x) else sort(unique(x))

    # If no categories (all NA), return 0-row df with same cols
    if (length(cats) == 0L) return(df[0, , drop = FALSE])

    # Build one row per category, casting the value to the same type as column x
    rows <- lapply(cats, function(cat) {
      r <- proto
      r[[var]] <- .cast_like(cat, x)
      r
    })
    dplyr::bind_rows(rows)
  }

  if (is.null(group_var)) {
    # -------- SINGLE-GROUP --------
    if (nrow(dat_original) == 0L) return(dat_original[0, , drop = FALSE])

    parts <- lapply(ov_ord, function(v) build_rows_for_var(dat_original, v))
    dummy <- dplyr::bind_rows(parts) %>% dplyr::distinct()
    if (nrow(dummy) == 0L) dummy <- dat_original[0, , drop = FALSE]
    return(dummy)
  } else {
    # -------- MULTI-GROUP --------
    if (!group_var %in% names(dat_original)) {
      stop("Group column '", group_var, "' not found in data supplied to create_dummy_from_original().")
    }
    if (nrow(dat_original) == 0L) return(dat_original[0, , drop = FALSE])

    dummy_list <- dat_original %>%
      dplyr::group_split(.data[[group_var]], .keep = TRUE) %>%
      lapply(function(df_g) {
        parts_g <- lapply(ov_ord, function(v) build_rows_for_var(df_g, v))
        base_g  <- df_g %>% dplyr::slice_head(n = 1)
        dplyr::bind_rows(base_g, parts_g) %>% dplyr::distinct()
      })

    dummy <- dplyr::bind_rows(dummy_list)

    # Guarantee at least one row per group
    have_groups <- unique(dummy[[group_var]])
    missing_g   <- setdiff(unique(dat_original[[group_var]]), have_groups)
    if (length(missing_g)) {
      add <- dat_original %>%
        dplyr::group_by(.data[[group_var]]) %>%
        dplyr::filter(.data[[group_var]] %in% missing_g) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::ungroup()
      dummy <- dplyr::bind_rows(dummy, add)
    }

    dplyr::distinct(dummy)
  }
}

# Cast a single scalar value to have the "same type" as a prototype vector x
.cast_like <- function(value, x) {
  if (is.factor(x)) {
    return(factor(value, levels = levels(x), ordered = is.ordered(x)))
  }
  if (inherits(x, "Date"))    return(as.Date(value, origin = "1970-01-01"))
  if (inherits(x, "POSIXct")) return(as.POSIXct(value, tz = attr(x, "tzone") %||% "UTC", origin = "1970-01-01"))
  if (is.integer(x))          return(as.integer(value))
  if (is.numeric(x))          return(as.numeric(value))
  if (is.logical(x))          return(as.logical(value))
  if (is.character(x))        return(as.character(value))
  value
}
