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
#' @examples
#' HS.model <- ' visual  =~ x1 + x2 + x3 '
#' fit <- lavaan::cfa(HS.model, data = lavaan::HolzingerSwineford1939)
#' is_single_group(fit)
#'
#' @export
is_single_group <- function(fit) {
  length(unique(fit@ParTable$group)) == 1
}
