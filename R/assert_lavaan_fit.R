#' Ensure object is a fitted lavaan model
#' @noRd
#' @keywords internal
.assert_lavaan_fit <- function(fit) {
  if (!inherits(fit, "lavaan")) {
    cls <- paste(class(fit), collapse = ", ")
    stop("`fit` must be a fitted lavaan model (class 'lavaan'); received: ",
         cls, call. = FALSE)
  }
  invisible(NULL)  # return nothing if ok
}
