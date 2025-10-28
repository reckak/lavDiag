#' Interactive CFA/SEM diagram (single estimate type) using visNetwork
#'
#' Builds interactive diagram(s) for a selected estimate type:
#'   - "none"   → raw (unstandardized)
#'   - "std.all"
#'   - "std.lv"
#' Returns one widget per group (or a single widget if single-group).
#'
#' @param fit A fitted lavaan/blavaan object.
#' @param standardized Either logical (TRUE → "std.all") or one of
#'   c("none","std.all","std.lv"). Default "none".
#' @param include_r2 Logical; include R-squared in node tooltips (computed
#'   from raw estimates via parameterEstimates with rsquare = TRUE). Default TRUE.
#'
#' @return A list of visNetwork htmlwidgets (one per group). For single-group,
#'   the list has length 1.
#' @export
#' @importFrom dplyr mutate transmute select filter arrange rename left_join
#' @importFrom dplyr group_by ungroup summarise across if_else case_when distinct
#' @importFrom dplyr where
#' @importFrom tidyr expand_grid nest drop_na
#' @importFrom purrr map map2_chr pmap
#' @importFrom tibble as_tibble tibble
#' @importFrom stringr str_c str_replace
plot_cfa <- function(fit,
                     standardized = "none",
                     include_r2 = TRUE) {

  # -- Basic checks ------------------------------------------------------------
  if (!inherits(fit, "lavaan")) {
    stop("`fit` must be a lavaan-compatible object.")
  }

  # Normalize standardized argument
  if (is.logical(standardized)) {
    std_type <- if (isTRUE(standardized)) "std.all" else "none"
  } else {
    std_type <- match.arg(standardized, c("none", "std.all", "std.lv"))
  }

  # -- Model info --------------------------------------------------------------
  info <- model_info(fit)

  # -- Meanstructure check (intercepts availability) ---------------------------
  meanstr <- tryCatch(lavaan::lavInspect(fit, "meanstructure"),
                      error = function(e) NA)
  has_intercepts <- isTRUE(meanstr)
  if (!has_intercepts) {
    warning("Model was not fitted with meanstructure = TRUE; intercepts will not be shown. ",
            "Refit the model with meanstructure = TRUE if you need intercepts.")
  }

  # -- Estimates: only one type (raw or standardized) --------------------------
  pe <- parameter_estimates(fit, standardized = std_type)

  # -- R2 (compute from raw PE regardless of chosen display type) --------------
  r2 <- NULL
  if (isTRUE(include_r2)) {
    pe_raw <- parameter_estimates(fit, standardized = "none")
    r2 <- pe_raw |>
      dplyr::filter(.data$op == "r2") |>
      dplyr::transmute(group = .data$group,
                       id    = .data$lhs,
                       r2    = stringr::str_c("<i>R\u00B2</i> = ", .fnum(.data$est))) |>
      dplyr::distinct()
  }

  # If single group, enforce group column to 1
  if (isTRUE(info$is_single_group)) {
    pe$group <- 1L
    if (!is.null(r2)) r2$group <- 1L
  }

  # -- Node labels -------------------------------------------------------------
  node_labels <- c(info$latent_variables, info$observed_variables)

  nodes <- tidyr::expand_grid(
    group = seq_len(info$n_groups),
    id    = node_labels
  ) |>
    dplyr::arrange(.data$group, .data$id) |>
    dplyr::mutate(
      latent = .data$id %in% info$latent_variables,
      label  = purrr::map2_chr(.data$id, .data$latent, .generate_node_labels),
      shape  = dplyr::if_else(.data$latent, "circle", "box"),
      `font.size` = dplyr::if_else(.data$latent, 40L, 30L),
      `color.background` = "white",
      `color.border`     = "black"
    )

  # -- (Residual) variances ----------------------------------------------------
  sigma <- pe |>
    dplyr::filter(.data$op == "~~", .data$lhs == .data$rhs) |>
    dplyr::mutate(sigma = stringr::str_c(
      "<i>(Resid.) Var</i> = ", .ci(.data$est, .data$ci.lower, .data$ci.upper)
    )) |>
    dplyr::select(.data$group, id = .data$lhs, .data$sigma)

  # -- Intercepts (if available) -----------------------------------------------
  nu <- NULL
  if (has_intercepts) {
    nu <- pe |>
      dplyr::filter(.data$op == "~1") |>
      dplyr::mutate(nu = stringr::str_c(
        "<i>Int</i> = ", .ci(.data$est, .data$ci.lower, .data$ci.upper)
      )) |>
      dplyr::select(.data$group, id = .data$lhs, .data$nu)
  }

  # -- Thresholds (categorical) -----------------------------------------------
  tau <- NULL
  if (isTRUE(info$is_categorical)) {
    tau <- pe |>
      dplyr::filter(.data$op == "|") |>
      dplyr::group_by(.data$group, id = .data$lhs) |>
      dplyr::summarise(
        tau = stringr::str_c(.fnum(.data$est), collapse = ", "),
        .groups = "drop"
      ) |>
      dplyr::mutate(tau = stringr::str_c("Thresholds: ", .data$tau))
  }

  # -- Join node info & build tooltips -----------------------------------------
  nodes <- nodes

  if (!is.null(r2)) {
    nodes <- dplyr::left_join(nodes, r2, by = c("group", "id"))
  }
  nodes <- dplyr::left_join(nodes, sigma, by = c("group", "id"))

  if (!is.null(nu)) {
    nodes <- dplyr::left_join(nodes, nu, by = c("group", "id"))
  }
  if (!is.null(tau)) {
    nodes <- dplyr::left_join(nodes, tau, by = c("group", "id"))
  }

  nodes <- nodes |>
    dplyr::rowwise() |>
    dplyr::mutate(
      title = {
        # Pick whatever of these columns exist:
        parts <- dplyr::c_across(dplyr::any_of(c("sigma", "nu", "tau", "r2")))
        parts <- parts[!is.na(parts) & nzchar(parts)]
        paste(parts, collapse = "<br>")
      }
    ) |>
    dplyr::ungroup()

  # -- Edges -------------------------------------------------------------------
  edges <- pe |>
    dplyr::transmute(from = .data$lhs,
                     to   = .data$rhs,
                     group = .data$group,
                     op    = .data$op,
                     est   = .data$est,
                     ci.lower = .data$ci.lower,
                     ci.upper = .data$ci.upper,
                     se = .data$se,
                     z  = .data$z,
                     pvalue = .data$pvalue) |>
    dplyr::filter(.data$from != .data$to, .data$op %in% c("=~", "~", "~~")) |>
    dplyr::mutate(
      arrows = ifelse(.data$op %in% c("=~", "~"), "to", "to;from"),
      hidden = FALSE,
      smooth = dplyr::if_else(.data$op == "~~", TRUE, FALSE),
      width  = .rescale(.data$est^2, endpoints = c(1, 5)),
      `color.opacity` = .rescale(.data$width, endpoints = c(0.2, 1)),
      dashes = dplyr::if_else(.data$op == "~~", TRUE, FALSE),
      label  = {
        z <- round(.data$est, 2)
        fmt <- format(z, nsmall = 2, trim = TRUE)
        stringr::str_replace(fmt, "-", "\u2013")
      },
      length = dplyr::if_else(.data$op == "=~", 150L, 50L),
      title  = stringr::str_c(
        "<i>Coef</i> = ", .fnum(.data$est), " [",
        .fnum(.data$ci.lower), ", ", .fnum(.data$ci.upper), "]<br>",
        "<i>SE</i> = ", .fnum(.data$se, digits = 3), "<br>",
        "<i>z</i> = ", .fnum(.data$z), ", ",
        "<i>p </i>", .fnum(.data$pvalue, pvalue = TRUE)
      )
    ) |>
    tibble::as_tibble()

  # -- Hidden helper edges to keep adjacent indicators aligned -----------------
  hidden <- edges |>
    dplyr::filter(.data$op == "=~") |>
    dplyr::group_by(.data$from, .data$group) |>
    dplyr::transmute(
      from   = .data$to,
      to     = dplyr::lead(.data$from),
      group  = .data$group,
      hidden = TRUE,
      width  = 1,
      smooth = FALSE
    ) |>
    tidyr::drop_na() |>
    dplyr::ungroup()

  edges <- dplyr::bind_rows(edges, hidden)

  # -- Nest per-group and prepare panel labels ---------------------------------
  nested_nodes <- nodes |>
    dplyr::group_by(.data$group) |>
    tidyr::nest(.key = "nodes")

  nested_edges <- edges |>
    dplyr::group_by(.data$group) |>
    tidyr::nest(.key = "edges")

  nested <- dplyr::left_join(nested_nodes, nested_edges, by = "group") |>
    dplyr::mutate(
      main = {
        tlabel <- switch(std_type,
                         "none"    = "Raw estimates",
                         "std.all" = "Standardized (std.all)",
                         "std.lv"  = "Standardized (std.lv)")
        if (isTRUE(info$is_single_group)) {
          tlabel
        } else {
          grp <- factor(.data$group,
                        levels = seq_len(info$n_groups),
                        labels = info$group_labels)
          stringr::str_c(tlabel, ", ", grp)
        }
      }
    ) |>
    dplyr::ungroup()

  # --- Compute coordinates in R (once on the first panel) ----------------------
  # vezmeme uzly/hrany z první skupiny (stejná sada id pro všechny skupiny)
  first_nodes <- nested$nodes[[1]]
  first_edges <- nested$edges[[1]]

  coord_nodes <- .compute_xy(first_nodes, first_edges)  # returns id,x,y

  # Reuse coordinates across all panels by joining on id
  nested <- nested |>
    dplyr::mutate(
      nodes = purrr::map(
        .data$nodes,
        ~ dplyr::left_join(.x, coord_nodes, by = "id")
      )
    )


  # -- Build and return all panels --------------------------------------------
  purrr::pmap(nested[c("nodes","edges","main")], .vis_network)
}
nodes = purrr::map
