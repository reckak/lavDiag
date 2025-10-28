#' Plot model-implied vs empirical item curves by latent factor (single- or multi-group)
#'
#' @description
#' Given the output list of `item_fit_data()`, draw smooth curves of the
#' model-implied relation and the empirical GAM-based relation between a chosen
#' latent factor (x-axis) and all items that load on it (y-axis). Works for
#' continuous, ordinal, and mixed models; supports single- and multi-group fits.
#'
#' @param x A list returned by `item_fit_data()`.
#' @param latent Character scalar: latent factor ID shown on x-axis.
#' @param items Optional character vector of item names to include; inferred if `NULL`.
#' @param show_model,show_empirical,show_points Logical toggles to draw the model curve,
#'   the empirical curve, and raw datapoints (all `TRUE` by default).
#' @param point_color,point_alpha,jitter_sd Appearance of datapoints (color, alpha, vertical jitter SD).
#' @param jitter_seed Optional integer seed for deterministic vertical jitter.
#' @param sample_frac Optional fraction in (0,1] to thin raw datapoints before plotting (default 1 = no thinning).
#' @param model_color,empirical_color Colors for model/empirical curves & ribbons.
#' @param ribbons Logical: draw CI ribbons when available (default `TRUE`).
#' @param penalized Logical: use penalized metrics in facet captions (default `FALSE`).
#' @param facet One of `"auto"`, `"wrap"`, `"grid"`. For multi-group inputs, a group grid is used
#'   regardless. For single-group, `"auto"` behaves like `"wrap"`.
#' @param ncol,nrow Optional layout hints used when `facet = "wrap"`.
#' @param scales Facet scales (passed to ggplot2 facets). Default `"free_y"`.
#'
#' @import ggplot2 dplyr tidyr stringr purrr rlang
#' @importFrom tidyselect any_of all_of
#' @export
item_fit_plot_old <- function(
    x,
    latent,
    items = NULL,
    show_model = TRUE,
    show_empirical = TRUE,
    show_points = TRUE,
    point_color = "black",
    point_alpha = 0.2,
    jitter_sd = 0.05,
    jitter_seed = NULL,
    sample_frac = 1,
    model_color = "blue",
    empirical_color = "red",
    ribbons = TRUE,
    penalized = FALSE,
    facet = c("auto","wrap","grid"),
    ncol = NULL,
    nrow = NULL,
    scales = c("fixed", "free_y","free", "free_x")
) {
  facet  <- match.arg(facet)
  scales <- match.arg(scales)
  stopifnot(is.null(jitter_seed) || (is.numeric(jitter_seed) && length(jitter_seed) == 1L && is.finite(jitter_seed)))
  if (!is.numeric(sample_frac) || length(sample_frac) != 1L || !is.finite(sample_frac) || sample_frac <= 0) sample_frac <- 1
  sample_frac <- min(sample_frac, 1)

  stopifnot(is.list(x), "original_data" %in% names(x))
  orig <- x$original_data
  new  <- x$new_data
  mets <- x$metrics
  if (is.null(new)) stop("x$new_data is NULL.")
  if (!latent %in% names(new)) stop(sprintf("Latent '%s' not found in new_data.", latent))

  if (is.null(items)) {
    rx <- paste0(
      "^(?:m|e)_(?:est|lwr|upr)_([^_]+)_",
      stringr::str_replace_all(latent, "([\\.^$|()\\[\\]{}*+?])", "\\\\\\1"),
      "$")
    cand  <- names(new)[stringr::str_detect(names(new), rx)]
    items <- sort(unique(stats::na.omit(stringr::str_match(cand, rx)[, 2])))
    if (!length(items)) stop("No <item>_<latent> prediction columns found in new_data for the requested latent.")
  }

  new_lat <- dplyr::filter(new, .latent_var %in% latent)
  if (!nrow(new_lat)) stop("No rows for the requested latent.")

  escaped_lat <- stringr::str_replace_all(latent, "([\\.^$|()\\[\\]{}*+?])", "\\\\\\1")
  cols_rx  <- paste0("^(?:m|e)_(?:est|lwr|upr)_[^_]+_", escaped_lat, "$")
  name_rx  <- paste0("^(m|e)_(est|lwr|upr)_([^_]+)_", escaped_lat, "$")

  core_cols <- c(".rid", ".gid", ".group", ".latent_var")
  pred_cols <- grep(cols_rx, names(new_lat), value = TRUE, perl = TRUE)
  if (!length(pred_cols)) stop("No prediction columns found for the chosen latent in new_data.")

  curves <- new_lat |>
    dplyr::mutate(x = .data[[latent]]) |>
    dplyr::select(dplyr::any_of(c(core_cols, "x")), tidyselect::all_of(pred_cols)) |>
    tidyr::pivot_longer(
      cols = tidyselect::all_of(pred_cols),
      names_to = c("src", "stat", "item"),
      names_pattern = name_rx,
      values_to = "val"
    ) |>
    tidyr::pivot_wider(
      names_from  = c("src", "stat"),
      values_from = "val",
      names_glue  = "{src}_{stat}"
    ) |>
    dplyr::mutate(
      .group = dplyr::if_else(is.na(.group) & is.finite(.gid), as.character(.gid), .group)
    )

  if (".group" %in% names(curves)) {
    ug <- unique(stats::na.omit(curves$.group))
    if (length(ug) <= 1L) curves$.group <- NULL
  }

  keep_items <- intersect(items, names(orig))
  if (length(keep_items) == 0L) show_points <- FALSE

  if (isTRUE(show_points)) {
    pts_long <- orig |>
      dplyr::select(dplyr::any_of(c(".gid", ".group", latent, keep_items))) |>
      tidyr::pivot_longer(cols = tidyselect::all_of(keep_items), names_to = "item", values_to = "y")

    if (sample_frac < 1 && nrow(pts_long) > 0) {
      pts_long <- dplyr::slice_sample(pts_long, prop = sample_frac)
    }

    if (".group" %in% names(curves)) {
      map_tbl <- dplyr::distinct(curves, .gid, .group)
      if (".group" %in% names(pts_long)) pts_long <- dplyr::select(pts_long, - .group)
      pts_long <- dplyr::left_join(pts_long, map_tbl, by = dplyr::join_by(.gid))
    }
  } else {
    pts_long <- NULL
  }

  p <- ggplot2::ggplot(curves, ggplot2::aes(x = x))

  if (isTRUE(show_points) && !is.null(pts_long) && nrow(pts_long) > 0) {
    p <- p + ggplot2::geom_point(
      data = pts_long,
      ggplot2::aes(x = !!rlang::sym(latent), y = y),
      inherit.aes = FALSE,
      color = point_color,
      alpha = point_alpha,
      shape = 16,
      show.legend = FALSE,
      position = ggplot2::position_jitter(height = jitter_sd, width = 0, seed = jitter_seed)
    )
  }

  if (isTRUE(show_model)) {
    if (isTRUE(ribbons) && any(!is.na(curves$m_lwr)))
      p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = m_lwr, ymax = m_upr, fill = "Model"), alpha = 0.2)
    p <- p + ggplot2::geom_line(ggplot2::aes(y = m_est, color = "Model"), linewidth = 0.9)
  }

  if (isTRUE(show_empirical)) {
    if (isTRUE(ribbons) && any(!is.na(curves$e_lwr)))
      p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = e_lwr, ymax = e_upr, fill = "Empirical"), alpha = 0.2)
    p <- p + ggplot2::geom_line(ggplot2::aes(y = e_est, color = "Empirical"), linewidth = 0.9)
  }

  if (!is.null(mets) && nrow(mets)) {
    met_cols <- if (isTRUE(penalized)) c("r2_pen", "rmse_pen", "mae_pen") else c("r2", "rmse", "mae")
    lab_tbl <- mets |>
      dplyr::filter(.data$item %in% items) |>
      dplyr::mutate(dplyr::across(dplyr::all_of(met_cols), ~ suppressWarnings(format(round(.x, 3), digits = 3, nsmall = 3)))) |>
      dplyr::transmute(
        .gid = .data$.gid,
        .group = .data$.group,
        item = .data$item,
        cap = paste0(
          "italic(R)^2==", .data[[met_cols[1]]], "*','~",
          "plain('RMSE')==", .data[[met_cols[2]]], "*','~",
          "plain('MAE')==", .data[[met_cols[3]]], "*'.'"
        )
      )

    if (".group" %in% names(curves)) {
      xr  <- curves |>
        dplyr::group_by(.group, item) |>
        dplyr::summarise(x_mid = mean(range(.data$x, na.rm = TRUE)), .groups = "drop")
      ann <- dplyr::left_join(lab_tbl, xr, by = dplyr::join_by(.group, item))
    } else {
      xr  <- curves |>
        dplyr::group_by(item) |>
        dplyr::summarise(x_mid = mean(range(.data$x, na.rm = TRUE)), .groups = "drop")
      ann <- dplyr::left_join(lab_tbl, xr, by = dplyr::join_by(item))
    }

    p <- p + ggplot2::geom_text(
      data = ann,
      ggplot2::aes(x = x_mid, y = Inf, label = cap),
      parse = TRUE,
      inherit.aes = FALSE,
      hjust = 0.5,
      vjust = 1.3
    ) + ggplot2::coord_cartesian(clip = "off")
  }

  p <- p +
    ggplot2::scale_color_manual(name = "Predictions",
                                values = c(Model = model_color, Empirical = empirical_color)) +
    ggplot2::scale_fill_manual(name = "Predictions",
                               values = c(Model = model_color, Empirical = empirical_color)) +
    ggplot2::labs(x = latent, y = NULL) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
    ggplot2::theme(legend.position = "top",
                   plot.margin = ggplot2::margin(10, 10, 20, 10))

  if (".group" %in% names(curves)) {
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(item), cols = ggplot2::vars(.group), scales = scales)
  } else {
    mode <- if (identical(facet, "auto")) "wrap" else facet
    if (identical(mode, "grid")) {
      p <- p + ggplot2::facet_grid(rows = ggplot2::vars(item), scales = scales)
    } else {
      if (is.null(ncol) && is.null(nrow)) {
        n_items <- length(unique(curves$item))
        ncol <- ceiling(sqrt(n_items))
      }
      p <- p + ggplot2::facet_wrap(ggplot2::vars(item), scales = scales, ncol = ncol, nrow = nrow)
    }
  }

  return(p)
}
