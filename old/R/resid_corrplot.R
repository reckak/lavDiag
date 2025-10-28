#' Corrplot of residual correlations (Bentler-type "cor")
#'
#' Draw corrplot(s) of residual correlations from a fitted `lavaan` model.
#' For multi-group models you can harmonize the color scale across groups via
#' `common_scale = TRUE`. Produces base plots.
#'
#' If `record = TRUE`, returns a recorded plot object (single-group) or a named
#' list of recorded plots (multi-group) created with `grDevices::recordPlot()`.
#' You can later replay them with `grDevices::replayPlot()`.
#'
#' @param fit A fitted `lavaan` object.
#' @param order One of `c("original","AOE","FPC","hclust","alphabet")`.
#' @param hclust.method One of
#'   `c("complete","ward","ward.D","ward.D2","single","average","mcquitty","median","centroid")`.
#' @param common_scale Logical; use common symmetric color range across groups? Default `TRUE`.
#' @param title_prefix Optional character prefix for multi-group plot titles.
#' @param record Logical; if `TRUE`, return recorded plot(s) via `recordPlot()`. Default `FALSE`.
#'
#' @return
#' - If `record = FALSE` (default): invisibly returns `NULL` (plots are drawn as a side-effect).
#' - If `record = TRUE`: a `recordedplot` (single-group) or a named `list` of `recordedplot`s (multi-group).
#'
#' @examples
#' HS.model <- '
#'   visual  =~ x1 + x2 + x3
#'   textual =~ x4 + x5 + x6
#'   speed   =~ x7 + x8 + x9
#' '
#' fit <- lavaan::cfa(HS.model, data = lavaan::HolzingerSwineford1939)
#' # draw only:
#' resid_corrplot(fit, order = "hclust")
#' # capture plot object:
#' rec <- resid_corrplot(fit, order = "hclust", record = TRUE)
#' # later:
#' if (interactive()) grDevices::replayPlot(rec)
#'
#' @importFrom lavaan lavResiduals lavInspect
#' @importFrom corrplot corrplot.mixed
#' @importFrom grDevices recordPlot
#' @export
resid_corrplot <- function(fit,
                           order = c("original", "AOE", "FPC", "hclust", "alphabet"),
                           hclust.method = c("complete", "ward", "ward.D", "ward.D2",
                                             "single", "average", "mcquitty", "median", "centroid"),
                           common_scale = TRUE,
                           title_prefix = NULL,
                           record = FALSE) {
  .assert_lavaan_fit(fit)

  # tvrdý požadavek na verzi s `col.lim`
  if (utils::packageVersion("corrplot") < "0.90") {
    stop("`resid_corrplot()` requires corrplot >= 0.90 (supports `col.lim`). ",
         "Installed: ", as.character(utils::packageVersion("corrplot")), call. = FALSE)
  }

  order         <- match.arg(order)
  hclust.method <- match.arg(hclust.method)

  rs <- lavaan::lavResiduals(fit, type = "cor")
  ng <- lavaan::lavInspect(fit, "ngroups")

  # vykresli jeden panel; vrátí buď NULL, nebo recorded plot
  draw_one <- function(mat, ttl = NULL, lim = NULL, record = FALSE) {
    mat <- as.matrix(mat)
    corrplot::corrplot.mixed(
      mat,
      order = order,
      hclust.method = hclust.method,
      is.corr = FALSE,
      cl.ratio = 0.1, tl.srt = 90, mar = c(0, 0, 2, 0),
      lower.col = "black",
      col.lim = if (!is.null(lim)) c(-lim, lim) else NULL,
      title = ttl
    )
    if (record) grDevices::recordPlot() else invisible(NULL)
  }

  if (identical(ng, 1L)) {
    mat <- rs$cov
    lim <- if (isTRUE(common_scale)) max(abs(mat), na.rm = TRUE) else NULL
    res <- draw_one(mat, ttl = NULL, lim = lim, record = record)
    if (record) return(res)
    return(invisible(NULL))
  }

  # Multi-group
  mats <- lapply(rs, function(lst) lst$cov)
  labs <- tryCatch(lavaan::lavInspect(fit, "group.label"), error = function(e) NULL)
  if (is.null(labs) || length(labs) != length(mats)) {
    labs <- names(mats)
    if (is.null(labs) || any(!nzchar(labs))) labs <- as.character(seq_along(mats))
  }

  lim <- if (isTRUE(common_scale)) max(abs(unlist(mats, use.names = FALSE)), na.rm = TRUE) else NULL

  if (!record) {
    for (i in seq_along(mats)) {
      ttl <- if (!is.null(title_prefix)) paste0(title_prefix, labs[i]) else labs[i]
      draw_one(mats[[i]], ttl = ttl, lim = lim, record = FALSE)
    }
    return(invisible(NULL))
  } else {
    recs <- vector("list", length(mats))
    names(recs) <- labs
    for (i in seq_along(mats)) {
      ttl <- if (!is.null(title_prefix)) paste0(title_prefix, labs[i]) else labs[i]
      recs[[i]] <- draw_one(mats[[i]], ttl = ttl, lim = lim, record = TRUE)
    }
    return(recs)
  }
}
