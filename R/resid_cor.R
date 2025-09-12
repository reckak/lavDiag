if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("v1","v2","v1c","v2c","pair","cor","se","z"))
}
#' Residual correlations (Bentler) as a tidy tibble
#'
#' Creates a tidy tibble of residual **correlations** (Bentler type) from a
#' fitted \code{lavaan} model, including standard errors and z-statistics when
#' available. Supports single- and multi-group models.
#'
#' @param fit A fitted \code{lavaan} object.
#'
#' @details
#' Internally uses \code{lavaan::lavResiduals(type = "cor.bentler", se = TRUE)}.
#' For multi-group models, a \code{group} column is added (using
#' \code{lavaan::lavInspect(fit, "group.label")} when available).
#' Only unique variable pairs are kept (upper triangle without the diagonal).
#'
#' @return
#' A tibble with columns:
#' \itemize{
#'   \item \code{v1}, \code{v2} – variable names in the pair
#'   \item \code{pair} – canonical pair label \code{"v1-v2"} with alphabetical ordering via \code{pmin/pmax}
#'   \item \code{cor} – residual correlation (Bentler)
#'   \item \code{abs_cor} – absolute value of \code{cor}
#'   \item \code{se} – standard error (if available from lavaan)
#'   \item \code{z} – z-statistic (if available)
#'   \item \code{group} – group label (multi-group models only)
#' }
#'
#' @examples
#' HS.model <- '
#'   visual  =~ x1 + x2 + x3
#'   textual =~ x4 + x5 + x6
#'   speed   =~ x7 + x8 + x9
#' '
#' fit <- lavaan::cfa(HS.model, data = lavaan::HolzingerSwineford1939)
#' resid_cor(fit)
#'
#' @importFrom lavaan lavResiduals lavInspect
#' @importFrom tibble as_tibble tibble
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom dplyr mutate across everything distinct rename relocate filter
#' @importFrom dplyr bind_rows arrange select
#' @importFrom rlang .data
#' @export
resid_cor <- function(fit) {

  is_not_lavaan_fit(fit)

  x <- lavaan::lavResiduals(fit, type = "cor.bentler", se = TRUE)

  tidy_one <- function(lst) {
    have_se <- !is.null(lst$cov.se)
    have_z  <- !is.null(lst$cov.z)

    # shodíme třídu lvn.mtr. -> čistý double
    cov_tbl <- tibble::as_tibble(as.matrix(lst$cov), .name_repair = "minimal") |>
      dplyr::mutate(v1 = colnames(lst$cov))

    long <- tidyr::pivot_longer(cov_tbl, -v1, names_to = "v2", values_to = "cor") |>
      dplyr::mutate(
        v1c = pmin(.data$v1, .data$v2),
        v2c = pmax(.data$v1, .data$v2),
        pair = paste0(.data$v1c, "-", .data$v2c)
      ) |>
      dplyr::filter(.data$v1 != .data$v2) |>
      dplyr::distinct(.data$pair, .keep_all = TRUE) |>
      dplyr::select(v1 = .data$v1c, v2 = .data$v2c, .data$pair, .data$cor)

    if (have_se) {
      se_tbl <- tibble::as_tibble(as.matrix(lst$cov.se), .name_repair = "minimal") |>
        dplyr::mutate(v1 = colnames(lst$cov.se)) |>
        tidyr::pivot_longer(-v1, names_to = "v2", values_to = "se") |>
        dplyr::mutate(
          v1 = pmin(.data$v1, .data$v2),
          v2 = pmax(.data$v1, .data$v2),
          pair = paste0(.data$v1, "-", .data$v2)
        ) |>
        dplyr::distinct(.data$pair, .keep_all = TRUE) |>
        dplyr::select(.data$pair, .data$se)
      long <- dplyr::left_join(long, se_tbl, by = "pair")
    }

    if (have_z) {
      z_tbl <- tibble::as_tibble(as.matrix(lst$cov.z), .name_repair = "minimal") |>
        dplyr::mutate(v1 = colnames(lst$cov.z)) |>
        tidyr::pivot_longer(-v1, names_to = "v2", values_to = "z") |>
        dplyr::mutate(
          v1 = pmin(.data$v1, .data$v2),
          v2 = pmax(.data$v1, .data$v2),
          pair = paste0(.data$v1, "-", .data$v2)
        ) |>
        dplyr::distinct(.data$pair, .keep_all = TRUE) |>
        dplyr::select(.data$pair, .data$z)
      long <- dplyr::left_join(long, z_tbl, by = "pair")
    }

    long |>
      dplyr::mutate(
        # jistota: všechny cílové sloupce jako double
        dplyr::across(dplyr::any_of(c("cor", "se", "z")), ~ as.double(.)),
        abs_cor = as.double(abs(.data$cor))
      ) |>
      dplyr::relocate(.data$abs_cor, .after = .data$cor)
  }

  ng <- lavaan::lavInspect(fit, "ngroups")
  if (identical(ng, 1L)) {
    tidy_one(x) |>
      dplyr::arrange(dplyr::desc(.data$abs_cor))
  } else {
    labs <- tryCatch(lavaan::lavInspect(fit, "group.label"), error = function(e) NULL)
    if (is.null(labs) || length(labs) != ng) labs <- as.character(seq_len(ng))

    lapply(seq_len(ng), function(g) {
      tg <- tidy_one(x[[g]])
      tg$group <- labs[g]
      tg
    }) |>
      dplyr::bind_rows() |>
      dplyr::arrange(.data$group, dplyr::desc(.data$abs_cor))
  }
}
