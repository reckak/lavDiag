#' Fast & simple parallel wrapper for lavaan::lavPredict (ordinal-only)
#'
#' Parallelization is used only when the model has ordinal indicators (ov.ord)
#' and workers > 1. For continuous-only models it falls back to a single
#' serial call (lavPredict is fast enough there).
#'
#' @param fit lavaan/blavaan object.
#' @param workers integer, number of workers (default: max(1, cores-1)).
#' @param plan one of c("auto","multisession","multicore","sequential").
#' @param chunk_size optional integer rows per chunk (default computed as ceiling(n_unique/workers)).
#' @param return_type "auto" | "data" | "list". For multi-group, "auto" -> list.
#' @param progress logical; show furrr progress bar.
#' @param ... extra args passed to lavaan::lavPredict().
#' @return tibble (single-group) or list/tibble (multi-group; see return_type).
lavPredict_parallel <- function(fit,
                                workers     = NULL,
                                plan        = c("auto","multisession","multicore","sequential"),
                                chunk_size  = NULL,
                                return_type = c("list","data"),
                                progress    = FALSE,
                                ...) {
  # -- Checks & meta -----------------------------------------------------------
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
  ov_ord <- ov_ord[ov_ord %in% names(dat_original)]
  has_ord <- length(ov_ord) > 0L

  # Decide parallelization: ordinal-only
  if (is.null(workers)) workers <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
  do_parallel <- has_ord && workers > 1L

  # -- Fast serial path for continuous-only -----------------------------------
  if (!do_parallel) {
    if (!is_mg) {
      out <- lavaan::lavPredict(fit,
                                newdata = dat_original,
                                append.data = TRUE,
                                assemble = TRUE,
                                ...)
      return(tibble::as_tibble(out))
    } else {
      out <- lavaan::lavPredict(fit,
                                newdata = dat_original,
                                append.data = TRUE,
                                assemble = TRUE,
                                drop.list.single.group = FALSE,
                                ...)
      out <- tibble::as_tibble(out)
      if (return_type == "data") return(out)
      grp <- out[[gvar]]; out[[gvar]] <- NULL
      return(split(out, grp))
    }
  }

  # -- Parallel path for ordinal models ---------------------------------------
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

  pred_fun <- function(df_chunk) {
    do.call(lavaan::lavPredict, c(list(object = fit, newdata = df_chunk),
                                  fopts, list(...)))
  }
  out_list <- furrr::future_map(chunked, pred_fun, .progress = progress)

  # 6) drop dummy rows
  dN <- nrow(dummy)
  for (i in seq_along(out_list)) {
    out_list[[i]] <- out_list[[i]][-(seq_len(dN)), , drop = FALSE]
  }

  # 7) bind & merge back to original
  out_all <- do.call(rbind, out_list) %>% tibble::as_tibble()
  out_joined <- dplyr::left_join(dat_original, out_all, by = colnames(dat_unique))

  if (!is_mg) return(out_joined)
  if (return_type == "data") return(out_joined)
  grp <- out_joined[[gvar]]; out_joined[[gvar]] <- NULL
  split(out_joined, grp)
}

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
    out <- dplyr::bind_rows(rows)
    out
  }

  if (is.null(group_var)) {
    # -------- SINGLE-GROUP --------
    # One template row must exist; if df empty, return 0-row
    if (nrow(dat_original) == 0L) return(dat_original[0, , drop = FALSE])

    # For each ordinal var, make rows covering all categories and bind
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

    # Work group-wise: for each group, build coverage rows for every ov.ord
    dummy_list <- dat_original %>%
      dplyr::group_split(.data[[group_var]], .keep = TRUE) %>%
      lapply(function(df_g) {
        # If a group is present, ensure at least one base row
        parts_g <- lapply(ov_ord, function(v) build_rows_for_var(df_g, v))
        # If some var had 0 categories (all NA), we still need at least one row for the group
        base_g  <- df_g %>% dplyr::slice_head(n = 1)
        dplyr::bind_rows(base_g, parts_g) %>% dplyr::distinct()
      })

    dummy <- dplyr::bind_rows(dummy_list)

    # Final safety: guarantee at least one row per group
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

    dummy <- dplyr::distinct(dummy)
    return(dummy)
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


# tiny infix helper
`%||%` <- function(x, y) if (is.null(x)) y else x
