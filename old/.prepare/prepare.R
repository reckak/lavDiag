#' Prepare: thin wrapper that works for continuous, ordinal, or mixed models
#'
#' @description
#' Minimal and robust unified wrapper:
#' - Tries to run `prepare_continuous(fit, ...)` and `prepare_ordinal(fit, ...)`.
#' - If only one succeeds, it returns that one.
#' - If both succeed, it merges them using `full_join()` over a stable ID key when possible.
#'
#' The goal is to avoid any manual manipulation of factor-score data or group columns.
#' All such handling is delegated to the respective sub-functions to minimize edge-case issues.
#'
#' @param fit A fitted lavaan/blavaan object.
#' @param data Optional pre-computed factor-score data; passed to sub-functions unchanged.
#' @param info Optional `model_info()` list; if `NULL`, it will be computed and forwarded.
#' @param plan Future plan for the ordinal branch; forwarded to `prepare_ordinal()`.
#' @param ... Additional arguments forwarded unchanged to both sub-functions (e.g.,
#'   `level`, `vcov_type`, `length.out`, `other_latents`, `latent_var_as_factor`,
#'   `se`, `se_summary`, etc.).
#' @return A tibble: either a single output (if the other is not applicable),
#'   or their `dplyr::full_join()` based on shared ID columns when available, otherwise
#'   on the full intersection of columns.
#' @export
prepare_old <- function(fit,
                    data = NULL,
                    info = NULL,
                    plan = c("auto","none","multisession","multicore","sequential","cluster"),
                    workers = NULL,
                    cluster = NULL,
                    ...) {
  # -- Match plan -------------------------------------------------------------
  plan <- match.arg(plan)

  # -- Compute/forward model_info --------------------------------------------
  if (is.null(info)) {
    info <- model_info(fit)
  }

  # -- Try both branches; fail-soft with NULL --------------------------------
  p_cont <- tryCatch(
    .prepare_continuous(fit, data = data, info = info, ...),
    error = function(e) NULL
  )
  p_ord <- tryCatch(
    .prepare_ordinal(fit, data = data, info = info, plan = plan,
                    workers = workers, cluster = cluster, ...),
    error = function(e) NULL
  )

  # -- If neither worked ------------------------------------------------------
  if (is.null(p_cont) && is.null(p_ord)) {
    rlang::abort("Neither continuous nor ordinal branch succeeded – check model/functions.")
  }

  # -- If only one worked, return it -----------------------------------------
  if (is.null(p_cont)) return(p_ord)
  if (is.null(p_ord))  return(p_cont)

  # -- Prefer a stable ID-based join when possible ---------------------------
  # Candidate ID columns commonly present across lavDiag outputs
  id_candidates <- c(".rid", ".gid", ".group", ".latent_var")
  ids_cont <- intersect(id_candidates, colnames(p_cont))
  ids_ord  <- intersect(id_candidates, colnames(p_ord))
  by_ids   <- intersect(ids_cont, ids_ord)

  # Fallback: use full intersection of columns if no stable IDs are shared
  by <- if (length(by_ids)) by_ids else intersect(colnames(p_cont), colnames(p_ord))

  if (!length(by)) {
    # Edge-case: no shared columns at all → row-bind with source tag
    rlang::warn("Continuous and ordinal outputs have no shared columns – returning row bind.")
    out <- dplyr::bind_rows(p_cont, p_ord, .id = ".source")
  } else {
    out <- dplyr::full_join(p_cont, p_ord, by = by)
  }

  # -- Relocate ID columns to the front for readability ----------------------
  id_front <- intersect(c(".rid", ".gid", ".group", ".latent_var"), colnames(out))
  if (length(id_front)) {
    out <- dplyr::relocate(out, dplyr::all_of(id_front))
  }

  out
}
