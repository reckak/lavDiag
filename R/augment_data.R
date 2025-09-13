#' Add predicted indicator values and residuals to factor-score output
#'
#' Computes \eqn{\hat{y} = \eta \Lambda^\top + \nu} for each observed indicator
#' and returns a data frame that contains (a) original lavaan factor scores
#' and observed variables, (b) predicted indicator values prefixed by
#' \code{prefix_yhat}, and (c) residuals (observed - predicted) prefixed by
#' \code{prefix_resid}. Works for single- and multi-group lavaan models with
#' continuous indicators.
#'
#' @param fit A fitted \code{lavaan} object.
#' @param prefix_yhat Character scalar, prefix for predicted columns. Default \code{".yhat_"}.
#' @param prefix_resid Character scalar, prefix for residual columns. Default \code{".resid_"}.
#'
#' @return A \code{tibble} with original factor-score output (including observed
#'   variables) plus columns of predicted indicators and their residuals. For
#'   multi-group models, a \code{group} column is included with group labels.
#'
#' @details
#' The function relies on \code{lavaan::lavPredict(..., transform = FALSE,
#' append.data = TRUE)} to obtain factor scores and observed variables.
#' Parameter estimates are taken from \code{lavInspect(fit, "est")} to extract
#' \eqn{\Lambda} (lambda) and intercepts \eqn{\nu} (nu). If the model was not
#' fitted with a meanstructure, intercepts are treated as zeros and a warning
#' is issued.
#'
#' @examples
#' # fit <- lavaan::cfa(model, data = dat, group = "grp", meanstructure = TRUE)
#' # out <- add_indicator_predictions(fit)
#' # dplyr::select(out, dplyr::starts_with(".yhat_"), dplyr::starts_with(".resid_"))
#'
#' @export
#' Add predicted indicator values and residuals for continuous indicators
#'
#' Computes \eqn{\hat{y} = \nu + \Lambda \eta} for observed **continuous**
#' indicators and returns a tibble with (a) original factor-score output
#' (including observed variables), (b) predicted values prefixed by
#' \code{prefix_yhat}, and (c) residuals (observed - predicted) prefixed by
#' \code{prefix_resid}. Works for single- and multi-group lavaan models.
#'
#' @param fit A fitted \code{lavaan} object with \code{meanstructure = TRUE}.
#' @param prefix_yhat Character scalar, prefix for predicted columns. Default \code{".yhat_"}.
#' @param prefix_resid Character scalar, prefix for residual columns. Default \code{".resid_"}.
#'
#' @return A \code{tibble}. For multi-group models, includes a \code{group}
#'   column with group labels. Predicted and residual columns are added **only**
#'   for continuous indicators; ordinal indicators are skipped.
#'
#' @details
#' The function requires intercepts (\eqn{\nu}) and loadings (\eqn{\Lambda})
#' to be available; it therefore enforces \code{meanstructure = TRUE}. It uses
#' \code{lavaan::lavPredict(..., transform = FALSE, append.data = TRUE)} to obtain
#' factor scores and observed variables, and \code{lavInspect(fit, "lambda")} and
#' \code{lavInspect(fit, "nu")} for measurement parameters.
#'
#' @examples
#' # fit <- lavaan::cfa('F =~ y1 + y2 + y3', data = dat, meanstructure = TRUE)
#' # out <- augment_indicator_predictions(fit)
#' # dplyr::select(out, dplyr::starts_with(".yhat_"), dplyr::starts_with(".resid_"))
#'
#' @export
augment <- function(fit,
                    prefix_yhat  = ".yhat_",
                    prefix_resid = ".resid_") {
  # Strict assertions (all backed by lavInspect(..., "est") only)
  .assert_lavaan_fit(
    fit,
    require_converged     = TRUE,
    require_meanstructure = TRUE,
    require_latent        = TRUE,
    require_measurement   = "lambda+nu",
    forbid_multilevel     = TRUE
  )

  # --- Model metadata ---------------------------------------------------------
  info    <- model_info(fit)
  ov_all  <- info$observed_variables
  ov_ord  <- lavaan::lavNames(fit, type = "ov.ord")
  ov_cont <- setdiff(ov_all, ov_ord)
  if (length(ov_cont) == 0L) stop("No continuous observed indicators detected; nothing to predict.", call. = FALSE)
  if (length(ov_ord) > 0L) {
    warning("Skipping ordinal indicators: ", paste(ov_ord, collapse = ", "),
            ". Predicted values and residuals are computed only for continuous indicators.")
  }

  n_groups     <- info$n_groups
  group_labels <- if (!is.null(info$group_labels)) info$group_labels else as.character(seq_len(n_groups))

  # --- Factor scores + observed variables per group ---------------------------
  fs_and_ov <- lavaan::lavPredict(
    fit, transform = FALSE, append.data = TRUE,
    assemble = FALSE, drop.list.single.group = FALSE
  )
  fs_list <- if (is.data.frame(fs_and_ov)) list(fs_and_ov) else fs_and_ov
  if (length(fs_list) != n_groups) {
    stop("lavPredict returned ", length(fs_list), " group table(s), expected ", n_groups, ".", call. = FALSE)
  }

  # --- Measurement parameters ONLY from lavInspect(..., "est") ----------------
  est <- lavaan::lavInspect(fit, "est")
  est_list <- if (is.list(est) && "lambda" %in% names(est)) list(est) else est

  # --- Extract eta-hat per group ---------------------------------------------
  eta_list <- lapply(fs_list, function(df) as.matrix(df[, info$latent_variables, drop = FALSE]))

  # --- Compute predictions for continuous indicators --------------------------
  yhat_list <- vector("list", n_groups)
  names(yhat_list) <- group_labels

  for (g in seq_len(n_groups)) {
    lam_g <- est_list[[g]]$lambda
    nu_g  <- est_list[[g]]$nu

    # Align Lambda: columns = latent order, rows = observed continuous order
    if (!is.null(colnames(lam_g))) lam_g <- lam_g[, info$latent_variables, drop = FALSE]
    if (!is.null(rownames(lam_g))) {
      lam_g <- lam_g[ov_cont, , drop = FALSE]
    } else {
      row_idx <- match(ov_cont, ov_all)
      lam_g <- lam_g[row_idx, , drop = FALSE]
    }

    # Align nu to continuous indicators
    if (!is.null(names(nu_g))) {
      nu_g <- nu_g[ov_cont]
    } else {
      row_idx <- match(ov_cont, ov_all)
      nu_g <- nu_g[row_idx]
    }

    eta_g  <- eta_list[[g]]                  # N x n_latent
    yhat_g <- eta_g %*% t(lam_g)             # N x n_cont
    yhat_g <- sweep(yhat_g, 2, nu_g, `+`)    # add intercepts
    colnames(yhat_g) <- paste0(prefix_yhat, ov_cont)
    yhat_list[[g]] <- tibble::as_tibble(yhat_g)
  }

  # --- Bind and compute residuals --------------------------------------------
  fs_tbl <- lapply(fs_list, tibble::as_tibble) |> dplyr::bind_rows(.id = "group")
  idx <- suppressWarnings(as.integer(fs_tbl$group))
  if (all(!is.na(idx)) && length(group_labels) >= max(idx)) fs_tbl$group <- group_labels[idx]

  yhat_tbl <- dplyr::bind_rows(yhat_list)
  out <- dplyr::bind_cols(fs_tbl, yhat_tbl)

  for (nm in ov_cont) {
    obs   <- nm
    yhat  <- paste0(prefix_yhat, nm)
    resid <- paste0(prefix_resid, nm)
    if (!obs %in% names(out))  stop("Observed variable '", obs,  "' not found in lavPredict output.", call. = FALSE)
    if (!yhat %in% names(out)) stop("Predicted column '", yhat, "' is missing; please report a bug.", call. = FALSE)
    out[[resid]] <- as.numeric(out[[obs]]) - as.numeric(out[[yhat]])
  }

  tibble::as_tibble(out)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
