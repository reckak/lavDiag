#' Ensure object is a fitted lavaan model
#'
#' Stops with an error if `fit` is not a lavaan object.
#'
#' @noRd
#' @keywords internal
is_not_lavaan_fit <- function(fit) {
  if (!inherits(fit, "lavaan")) {
    cls <- paste(class(fit), collapse = ", ")
    stop("`fit` must be a fitted lavaan model (class 'lavaan'); received: ",
         cls, call. = FALSE)
  }
  invisible(FALSE)
}
