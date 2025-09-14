#' Parallelized lavPredict for Ordinal Models
#'
#' A parallelized version of lavaan::lavPredict optimized for models with ordinal
#' variables. Supports both single-group and multi-group models.
#'
#' @param object A fitted lavaan object
#' @param newdata Optional data.frame with new data for prediction. If NULL,
#'   uses the original model data.
#' @param chunk_size Integer. Number of rows per chunk for parallel processing.
#'   Default is 100.
#' @param cores Integer. Number of cores to use for parallel processing.
#'   Default is parallel::detectCores() - 1.
#' @param by_group Logical. For multi-group models, whether to chunk within
#'   each group separately (TRUE) or across all groups (FALSE). Default is TRUE.
#' @param ... Additional arguments passed to lavaan::lavPredict
#'
#' @return Predictions in the same format as lavaan::lavPredict
#'
#' @details
#' This function addresses the challenge of using lavPredict with ordinal variables
#' in parallel processing by ensuring each chunk contains all levels of ordinal
#' variables through dummy rows.
#'
#' For multi-group models, two strategies are available:
#' \itemize{
#'   \item by_group = TRUE: Process each group separately (recommended for groups
#'     with different ordinal levels)
#'   \item by_group = FALSE: Process all groups together (faster for similar groups)
#' }
#'
#' @examples
#' \dontrun{
#' # Single group model
#' predictions <- parallel_lavPredict(fit, chunk_size = 100)
#'
#' # Multi-group model
#' predictions <- parallel_lavPredict(fit_mg, by_group = TRUE)
#' }
#'
#' @import future
#' @import furrr
#' @import dplyr
#' @importFrom tidyr unnest
#' @importFrom lavaan lavInspect lavNames lavPredict
#' @importFrom purrr map transpose
#' @importFrom parallel detectCores
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
#' @export
parallel_lavPredict <- function(object,
                                newdata = NULL,
                                chunk_size = 100L,
                                cores = parallel::detectCores() - 1L,
                                by_group = TRUE,
                                ...) {

  # Input validation
  if (!inherits(object, "lavaan")) {
    stop("'object' must be a fitted lavaan model", call. = FALSE)
  }

  if (!is.null(newdata) && !is.data.frame(newdata)) {
    stop("'newdata' must be a data.frame or NULL", call. = FALSE)
  }

  chunk_size <- as.integer(chunk_size)
  if (chunk_size < 1) {
    stop("'chunk_size' must be a positive integer", call. = FALSE)
  }

  cores <- as.integer(cores)
  if (cores < 1) {
    stop("'cores' must be a positive integer", call. = FALSE)
  }

  # Check required packages
  check_packages(c("future", "furrr"))

  # Set up parallelization
  current_plan <- future::plan()
  if (!inherits(current_plan, "sequential")) {
    message("Future plan already set, using current configuration")
  } else {
    future::plan(future::multisession, workers = cores)
    on.exit(future::plan(future::sequential), add = TRUE)
  }

  # Get data
  if (is.null(newdata)) {
    dat <- lavaan::lavInspect(object, "data")
    if (is.list(dat)) {
      # Multi-group data is a list of matrices
      dat <- lapply(dat, tibble::as_tibble)
    } else {
      # Single group data is a matrix
      dat <- tibble::as_tibble(dat)
    }
  } else {
    dat <- tibble::as_tibble(newdata)
  }

  # Check if multigroup model
  ngroups <- lavaan::lavInspect(object, "ngroups")
  is_multigroup <- ngroups > 1L

  if (is_multigroup) {
    group_labels <- lavaan::lavInspect(object, "group.label")
    message(sprintf("Detected multi-group model with %d groups: %s",
                    ngroups, paste(group_labels, collapse = ", ")))
  }

  # Identify ordinal variables
  ov_ord <- lavaan::lavNames(object, "ov.ord")

  if (length(ov_ord) == 0L) {
    message("No ordinal variables found. Using standard lavPredict.")
    return(lavaan::lavPredict(object, newdata = newdata, ...))
  }

  # Process based on model type and strategy
  if (is_multigroup && by_group) {
    process_multigroup_by_group(object, dat, ov_ord, chunk_size, ...)
  } else {
    process_combined(object, dat, ov_ord, chunk_size, ...)
  }
}

#' Check Required Packages
#' @param packages Character vector of package names
#' @keywords internal
check_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace,
                              FUN.VALUE = logical(1L), quietly = TRUE)]

  if (length(missing) > 0L) {
    stop(sprintf("Required packages not available: %s",
                 paste(missing, collapse = ", ")), call. = FALSE)
  }
}

#' Process Multi-group Model by Group
#' @param object lavaan object
#' @param dat data
#' @param ov_ord ordinal variable names
#' @param chunk_size chunk size
#' @param ... additional arguments
#' @keywords internal
process_multigroup_by_group <- function(object, dat, ov_ord, chunk_size, ...) {

  group_labels <- lavaan::lavInspect(object, "group.label")
  group_var <- lavaan::lavInspect(object, "group")

  if (is.list(dat)) {
    # Data is already split by groups
    group_results <- vector("list", length(dat))
    names(group_results) <- group_labels

    for (g in seq_along(dat)) {
      group_data <- dat[[g]]
      message(sprintf("Processing group %d (%s) with %d rows",
                      g, group_labels[g], nrow(group_data)))

      group_results[[g]] <- process_single_group(
        object, group_data, ov_ord, chunk_size, g, ...)
    }

    return(group_results)

  } else if (!is.null(group_var)) {
    # Data has group variable
    group_results <- vector("list", length(group_labels))
    names(group_results) <- group_labels

    for (g in seq_along(group_labels)) {
      group_data <- dat %>%
        dplyr::filter(.data[[group_var]] == group_labels[g])

      message(sprintf("Processing group %d (%s) with %d rows",
                      g, group_labels[g], nrow(group_data)))

      if (nrow(group_data) > 0L) {
        group_results[[g]] <- process_single_group(
          object, group_data, ov_ord, chunk_size, g, ...)
      }
    }

    return(group_results)

  } else {
    stop("Cannot identify group structure for multi-group model", call. = FALSE)
  }
}

#' Process Single Group
#' @param object lavaan object
#' @param group_data data for this group
#' @param ov_ord ordinal variable names
#' @param chunk_size chunk size
#' @param group_id group identifier
#' @param ... additional arguments
#' @keywords internal
process_single_group <- function(object, group_data, ov_ord, chunk_size, group_id, ...) {

  n_rows <- nrow(group_data)
  if (n_rows == 0L) return(NULL)

  # Create chunks
  idx_list <- split(seq_len(n_rows), ceiling(seq_len(n_rows) / chunk_size))

  message(sprintf("  Splitting into %d chunks (size: %d)",
                  length(idx_list), chunk_size))

  # Create dummy rows for ordinal variables
  dummy_data <- create_dummy_rows(group_data, ov_ord)

  # Process chunks in parallel
  predictions <- process_chunks_parallel(object, group_data, dummy_data,
                                         idx_list, group_id, ...)

  # Combine results
  combine_predictions(predictions)
}

#' Process Combined (All Groups Together)
#' @param object lavaan object
#' @param dat data
#' @param ov_ord ordinal variable names
#' @param chunk_size chunk size
#' @param ... additional arguments
#' @keywords internal
process_combined <- function(object, dat, ov_ord, chunk_size, ...) {

  # Handle list data (multi-group)
  if (is.list(dat)) {
    dat <- do.call(rbind, dat)
  }

  # Ensure dat is a tibble
  if (is.matrix(dat)) {
    dat <- tibble::as_tibble(dat)
  }

  n_rows <- nrow(dat)
  idx_list <- split(seq_len(n_rows), ceiling(seq_len(n_rows) / chunk_size))

  message(sprintf("Processing %d rows in %d chunks (size: %d)",
                  n_rows, length(idx_list), chunk_size))

  # Create dummy rows
  dummy_data <- create_dummy_rows(dat, ov_ord)

  # Process chunks in parallel
  message("Running parallel predictions...")

  predictions <- process_chunks_parallel(object, dat, dummy_data, idx_list,
                                         group_id = NULL, ...)

  # Combine results
  result <- combine_predictions(predictions)

  message("Parallel prediction completed!")
  result
}

#' Create Dummy Rows for Ordinal Variables
#' @param data data.frame
#' @param ov_ord ordinal variable names
#' @return data.frame with dummy rows
#' @keywords internal
create_dummy_rows <- function(data, ov_ord) {

  # Get unique values for ordinal variables
  unique_valid <- function(x) {
    out <- unique(x)
    out[!is.na(out)]
  }

  dummy_ord <- data %>%
    dplyr::select(dplyr::all_of(ov_ord)) %>%
    furrr::future_map(unique_valid, .options = furrr::furrr_options(seed = TRUE)) %>%
    dplyr::bind_cols()

  if (nrow(dummy_ord) == 0L) {
    return(data[integer(0), , drop = FALSE])
  }

  # Create complete dummy rows
  dummy_full <- data[seq_len(nrow(dummy_ord)), , drop = FALSE] %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(ov_ord),
                                ~ dummy_ord[[dplyr::cur_column()]]))

  # Handle non-ordinal variables with mean values
  non_ord_vars <- setdiff(names(data), ov_ord)
  if (length(non_ord_vars) > 0L) {
    means <- data %>%
      dplyr::select(dplyr::all_of(non_ord_vars)) %>%
      dplyr::summarise(dplyr::across(dplyr::everything(),
                                     ~ mean(.x, na.rm = TRUE)))

    dummy_full <- dummy_full %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(non_ord_vars),
                                  ~ rep(means[[dplyr::cur_column()]],
                                        nrow(dummy_full))))
  }

  dummy_full
}

#' Process Chunks in Parallel
#' @param object lavaan object
#' @param data original data
#' @param dummy_data dummy rows
#' @param idx_list list of row indices for chunks
#' @param group_id group identifier (NULL for combined processing)
#' @param ... additional arguments
#' @keywords internal
process_chunks_parallel <- function(object, data, dummy_data, idx_list, group_id, ...) {

  n_dummy <- nrow(dummy_data)

  # Prepare chunks with dummy rows
  chunked_data <- idx_list %>%
    furrr::future_map(function(idx) {
      chunk_dat <- data[idx, , drop = FALSE]
      # Ensure both are tibbles/data.frames
      if (is.matrix(chunk_dat)) {
        chunk_dat <- tibble::as_tibble(chunk_dat)
      }
      if (is.matrix(dummy_data)) {
        dummy_data <- tibble::as_tibble(dummy_data)
      }
      dplyr::bind_rows(dummy_data, chunk_dat)
    }, .options = furrr::furrr_options(seed = TRUE))

  # Parallel prediction
  predictions <- chunked_data %>%
    furrr::future_map(function(chunk) {
      pred <- lavaan::lavPredict(object, newdata = chunk, ...)

      # Remove dummy rows from results
      if (is.list(pred) && !is.null(group_id)) {
        # Multi-group output, extract relevant group
        pred <- pred[[group_id]]
      }

      if (is.matrix(pred) || is.vector(pred)) {
        if (n_dummy > 0L) {
          pred[seq_len(nrow(pred))[-seq_len(n_dummy)], , drop = FALSE]
        } else {
          pred
        }
      } else {
        pred
      }
    }, .options = furrr::furrr_options(seed = TRUE))

  predictions
}

#' Combine Prediction Results
#' @param predictions list of prediction results
#' @keywords internal
combine_predictions <- function(predictions) {

  if (length(predictions) == 0L) {
    return(NULL)
  }

  first_pred <- predictions[[1]]

  if (is.matrix(first_pred) || is.vector(first_pred)) {
    do.call(rbind, predictions)
  } else if (is.list(first_pred)) {
    # Multi-group case
    predictions %>%
      purrr::transpose() %>%
      purrr::map(~ do.call(rbind, .x))
  } else {
    unlist(predictions)
  }
}
