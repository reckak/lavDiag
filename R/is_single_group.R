#' Check if a lavaan model is single-group
#'
#' This helper function checks whether a fitted \code{lavaan} model
#' contains only a single group.
#'
#' @param fit An object of class \code{lavaan}.
#'
#' @return Logical scalar. Returns \code{TRUE} if the model has only
#'   one group, \code{FALSE} otherwise.
#'
#' @noRd
#' @keywords internal
.is_single_group <- function(fit) {
  length(unique(fit@ParTable$group)) == 1
}
