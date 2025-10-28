#' Augment ordinal SEM output with factor scores, latent linear predictors (y*),
#' per–category probabilities, and expected observed scores (with optional CIs)
#'
#' @description
#' Extends the factor-score output of an ordinal \pkg{lavaan} model by adding:
#' (i) latent linear predictors \eqn{y^\*} for each ordinal indicator,
#' (ii) per–category probabilities for every observed (or threshold-inferred) level, and
#' (iii) expected observed scores computed from those probabilities; optionally
#' with delta-method confidence intervals that reflect **parameter** uncertainty.
#'
#' Supports single- and multi-group models with **only ordinal indicators**.
#'
#' @details
#' For each ordinal indicator \eqn{j}, the function computes
#' \deqn{y^\*_j = \nu_j + \lambda_{j\cdot}\,\eta,}
#' using estimates from `lavInspect(fit, "est")`. Under the `theta`
#' parameterization, residual variances are used to scale the probit links.
#'
#' Category probabilities are obtained from thresholded normal CDFs and are named
#' \code{<prefix_pr><level><sep><item>}. Expected observed scores
#' \code{<prefix_hat><item>} are probability-weighted means of category scores.
#' If factor levels are numeric-like (e.g., "1","2", ...), they are used
#' directly; otherwise a \code{1..K} sequence is used.
#'
#' If \code{ci = "delta"}, pointwise CIs for expected scores are computed via
#' the multivariate delta method with \code{lavInspect(fit, "vcov")}. These
#' intervals reflect uncertainty in thresholds, loadings, intercepts (and
#' residual variances in `theta`), but **not** factor-score indeterminacy or
#' sampling variability of new cases.
#'
#' Multi-group models are handled by aligning parameter blocks per group and
#' stacking outputs. Factor scores come from a parallelized wrapper around
#' \code{lavaan::lavPredict()}.
#'
#' @section Naming:
#' \itemize{
#'   \item \code{<prefix_ystar><item>} — latent linear predictor \eqn{y^\*}.
#'   \item \code{<prefix_pr><level><sep><item>} — per-category probability.
#'   \item \code{<prefix_hat><item>} — expected observed score.
#'   \item With \code{ci = "delta"} also:
#'         \code{<prefix_ci[1]><item>} (lower CI),
#'         \code{<prefix_ci[2]><item>} (upper CI).
#' }
#'
#' @param fit            A fitted \pkg{lavaan} object for an ordinal (ordered-categorical) SEM.
#' @param prefix_ystar   Prefix for latent linear predictors \eqn{y^\*}. Default \code{".ystar_"}.
#' @param prefix_pr      Prefix for per-category probabilities. Default \code{".pr_"}.
#' @param prefix_hat     Prefix for expected observed scores. Default \code{".yhat_"}.
#' @param prefix_ci      A length-2 character vector of prefixes for lower and upper
#'   CIs of expected scores (used when \code{ci = "delta"}). Default
#'   \code{c(".yhat_lwr_", ".yhat_upr_")}.
#' @param sep            Separator between \code{<level>} and \code{<item>} in probability names.
#'   Default \code{"__"}.
#' @param ci             Confidence-interval method for expected scores:
#'   \code{"none"} (default) or \code{"delta"}.
#' @param level          Confidence level for \code{ci = "delta"}. Default \code{0.95}.
#'
#' @returns
#' A tibble with:
#' \itemize{
#'   \item factor scores (as returned by \code{lavPredict()}),
#'   \item \code{.gid} (numeric group id) and \code{group} (group label),
#'   \item \code{<prefix_ystar><item>} for each ordinal indicator,
#'   \item \code{<prefix_hat><item>} and, if requested, its lower/upper CI via \code{prefix_ci},
#'   \item all per-category probability columns \code{<prefix_pr><level><sep><item>}.
#' }
#'
#' @seealso
#' \code{\link[lavaan]{lavPredict}}, internal checks \code{.assert_lavaan_fit()},
#' model summary helper \code{model_info()}, and the parallel wrapper
#' \code{lavPredict_parallel()}.
#'
#' @examples
#' \dontrun{
#' # fit <- lavaan::cfa(..., ordered = c("y1","y2","y3"))
#' aug <- augment_ordinal_3(fit, ci = "delta")
#' dplyr::glimpse(aug)
#' }
#' @export


augment_ordinal_3 <- function(fit,
                              prefix_ystar    = ".ystar_",
                              prefix_pr       = ".pr_",
                              prefix_yhat     = ".yhat_",
                              prefix_ci       = c(".yhat_lwr_", ".yhat_upr_"),
                              sep             = "__",
                              ci              = c("none","delta"),
                              level           = 0.95
                              ) {


  # -- Strict assertions -------------------------------------------------------
  # requires: converged fit, latent variables, and measurement (lambda + nu)

  .assert_lavaan_fit(
    fit,
    require_converged     = TRUE,
    require_latent        = TRUE,
    require_meanstructure = NA,
    require_measurement   = "lambda+nu",
    forbid_multilevel     = FALSE
  )

  ci <- match.arg(ci)

  # -- Model info --------------------------------------------------------------

  info        <- model_info(fit)
  ov_all      <- info$observed_variables
  ov_ord      <- info$ov_ordinal
  n_groups    <- info$n_groups %||% 1L
  group_labels<- info$group_labels %||% as.character(seq_len(n_groups))
  param       <- info$parameterization %||% "delta"

  fs_and_ov <- lavPredict_parallel(fit, return_type = "list")
  fs_list   <- if (is.data.frame(fs_and_ov)) list(fs_and_ov) else fs_and_ov

  est        <- lavaan::lavInspect(fit, "est")
  est_list   <- if (is.list(est) && "lambda" %in% names(est)) list(est) else est

  th_raw     <- tryCatch(lavaan::lavInspect(fit, "th"), error = function(e) NULL)
  th_list    <- if (is.null(th_raw)) replicate(n_groups, setNames(numeric(0), character(0)), simplify = FALSE)
  else if (is.list(th_raw)) th_raw else list(th_raw)
  if (length(th_list) == 1L && n_groups > 1L) th_list <- rep(th_list, n_groups)
  names(th_list) <- group_labels

  # -- Cache for delta method --------------------------------------------------
  coef_names <- names(lavaan::coef(fit))
  V_full     <- lavaan::lavInspect(fit, "vcov")
  zcrit      <- stats::qnorm(1 - (1 - level)/2)

  # -- Helpers -----------------------------------------------------------------
  # Align Lambda/nu to given indicators
  .align_Lambda_nu <- function(lam_g, nu_g, ov_subset, ov_all, latent_order) {
    # columns = latent variables: enforce order
    if (!is.null(colnames(lam_g))) lam_g <- lam_g[, latent_order, drop = FALSE]
    # rows = observed variables
    if (!is.null(rownames(lam_g))) {
      lam_g <- lam_g[ov_subset, , drop = FALSE]
    } else {
      row_idx <- match(ov_subset, ov_all); lam_g <- lam_g[row_idx, , drop = FALSE]
    }
    # intercepts/threshold centers (nu)
    if (is.null(nu_g) || length(nu_g) == 0L) {
      nu_g <- rep(0, length(ov_subset)); names(nu_g) <- ov_subset
    } else if (!is.null(names(nu_g))) {
      nu_g <- nu_g[ov_subset]
    } else {
      row_idx <- match(ov_subset, ov_all); nu_g <- nu_g[row_idx]
    }
    list(lambda = lam_g, nu = nu_g)
  }

  # Extract thresholds for variable j, ordered by numeric suffix t1 < t2 < ...
  .get_tau_for <- function(th_vec, j) {
    nm <- names(th_vec); if (is.null(nm)) return(numeric(0))
    sel <- startsWith(nm, paste0(j, "|"))
    tau <- th_vec[sel]
    if (!length(tau)) return(numeric(0))
    # order by numeric index after '|t'
    idx <- suppressWarnings(as.integer(gsub("^.*\\|t(\\d+)$", "\\1", names(tau))))
    idx[is.na(idx)] <- seq_along(tau)
    tau[order(idx)]
  }

  # Try to get levels from appended data; otherwise infer 1..m from thresholds
  .levels_for_item <- function(df, th_vec, j) {
    if (!is.null(df) && j %in% names(df)) {
      yo <- df[[j]]
      if (is.factor(yo)) {
        return(levels(yo))
      } else {
        lev <- sort(unique(yo))
        return(as.character(lev))
      }
    } else {
      tau_j <- .get_tau_for(th_vec, j)
      if (!length(tau_j)) return(character(0))
      m <- length(tau_j) + 1L
      return(as.character(seq_len(m)))
    }
  }

  # Derive numeric scores for levels (use numeric-like labels if possible)
  .scores_from_levels <- function(lev_chr) {
    if (!length(lev_chr)) return(numeric(0))
    sc <- suppressWarnings(as.numeric(lev_chr))
    if (anyNA(sc)) {
      # Fallback to 1..K when levels are non-numeric (e.g., "Never","Sometimes",...)
      sc <- seq_along(lev_chr)
    }
    sc
  }

  # Predicate for gentle numeric coercion (do not convert integers)
  .needs_double_from_special <- function(x) {
    inherits(x, "lvn.mtrx") || inherits(x, "lavaan.matrix")
  }

  ystar_list <- vector("list", n_groups)
  probs_list <- vector("list", n_groups)
  # we will build yexp+CI interleaved per item; keep as list-of-lists
  yexp_ci_list <- vector("list", n_groups)

  # latent scores per group
  eta_list <- lapply(fs_list, function(df) as.matrix(df[, info$latent_variables, drop = FALSE]))

  for (g in seq_len(n_groups)) {
    est_g <- est_list[[g]]
    lam_g <- est_g$lambda
    nu_g  <- est_g$nu %||% numeric(0)

    al <- .align_Lambda_nu(lam_g, nu_g,
                           ov_subset    = ov_ord,
                           ov_all       = ov_all,
                           latent_order = info$latent_variables)

    eta_g    <- eta_list[[g]]                                 # N x q
    ystar_g  <- sweep(eta_g %*% t(al$lambda), 2, al$nu, `+`)  # N x p_ord
    colnames(ystar_g) <- ov_ord

    # residual SDs in theta parametrization
    sig_eps <- NULL
    if (identical(param, "theta")) {
      theta_g <- est_g$theta
      if (!is.null(theta_g)) {
        if (!is.null(rownames(theta_g))) {
          sig_eps <- sqrt(diag(theta_g)[ov_ord])
        } else {
          idx <- match(ov_ord, ov_all); sig_eps <- sqrt(diag(theta_g)[idx])
        }
        sig_eps[!is.finite(sig_eps) | sig_eps <= 0] <- 1
      } else sig_eps <- rep(1, length(ov_ord))
    }

    N_g <- nrow(ystar_g)

    probs_g    <- NULL
    ystar_df   <- tibble::as_tibble(ystar_g)
    names(ystar_df) <- paste0(prefix_ystar, names(ystar_df))

    # interleaved yexp + CI accumulator for this group (bind by columns at the end)
    yexp_ci_cols <- list()

    # factor names aligned to lambda columns
    lam_cols <- if (!is.null(colnames(al$lambda))) colnames(al$lambda) else info$latent_variables

    for (jj in seq_along(ov_ord)) {
      j       <- ov_ord[jj]
      ystar_j <- ystar_g[, jj]

      # levels & thresholds
      lev   <- .levels_for_item(fs_list[[g]], th_list[[g]], j)
      tau_j <- as.numeric(.get_tau_for(th_list[[g]], j))
      if (length(tau_j) && length(lev) != (length(tau_j) + 1L)) lev <- as.character(seq_len(length(tau_j) + 1L))
      if (!length(lev) || !length(tau_j)) next

      m       <- length(lev)
      tau_pad <- c(-Inf, tau_j, Inf)

      # z-matrix: N x (m+2)
      if (identical(param, "theta")) {
        se_j <- if (is.null(sig_eps)) 1 else sig_eps[jj]
        zmat <- outer(ystar_j, tau_pad, function(y, t) (t - y) / se_j)
      } else {
        zmat <- outer(ystar_j, tau_pad, function(y, t) (t - y))
      }

      Phi     <- stats::pnorm(zmat)
      probs_j <- Phi[, -1, drop = FALSE] - Phi[, -ncol(Phi), drop = FALSE]  # N x m
      colnames(probs_j) <- paste0(prefix_pr, make.names(lev, allow_ = TRUE), sep, j)
      probs_g <- if (is.null(probs_g)) as.data.frame(probs_j, check.names = FALSE) else
        cbind(probs_g, as.data.frame(probs_j, check.names = FALSE))

      # expected scores
      scores_j <- { sc <- suppressWarnings(as.numeric(lev)); if (anyNA(sc)) seq_len(m) else sc }
      yexp_j   <- as.numeric(probs_j %*% scores_j)  # N-vector

      if (ci == "delta") {
        # --- Vectorized gradients --------------------------------------------
        phi_mat    <- stats::dnorm(zmat)                  # N x (m+2)
        phi_inner  <- phi_mat[, -1, drop = FALSE]         # z_k   (cols 2..m+1)
        phi_prev   <- phi_mat[, -ncol(phi_mat), drop = FALSE]  # z_{k-1} (cols 1..m)

        # Base term B = sum_k s_k * (phi(z_k) - phi(z_{k-1}))  -> N-vector
        # (matrix multiply N x m by m x 1)
        B <- as.numeric((phi_inner[, 1:m, drop = FALSE] - phi_prev[, 1:m, drop = FALSE]) %*% matrix(scores_j, ncol = 1))

        # dmu/d nu_j  (constant dz_nu factor across k)
        if (identical(param, "theta")) {
          se_j <- if (is.null(sig_eps)) 1 else sig_eps[jj]
          d_nu <- (-1 / se_j) * B
        } else {
          d_nu <- (-1) * B
        }

        # dmu/d lambda_jr for each factor r: outer(B, dz_lam_row)
        eta_mat    <- eta_g[, lam_cols, drop = FALSE]      # N x q
        if (identical(param, "theta")) {
          dz_lam_mat <- -eta_mat / (if (is.null(sig_eps)) 1 else sig_eps[jj])  # N x q
        } else {
          dz_lam_mat <- -eta_mat
        }
        d_lam_mat <- dz_lam_mat * B  # broadcasting by rows (N x q)

        # dmu/d tau_{jr}, r=1..m-1:  dz_tau * ( (s_r - s_{r+1}) * phi(z_r) )
        s_diff    <- scores_j[1:(m - 1)] - scores_j[2:m]   # length m-1
        phi_r_mat <- phi_mat[, 1 + seq_len(m - 1), drop = FALSE]  # columns for z_r (2..m)
        if (identical(param, "theta")) {
          se_j <- if (is.null(sig_eps)) 1 else sig_eps[jj]
          d_tau_mat <- sweep(phi_r_mat, 2, (s_diff / se_j), `*`)  # N x (m-1)
        } else {
          d_tau_mat <- sweep(phi_r_mat, 2, s_diff, `*`)           # N x (m-1)
        }

        # dmu/d theta_jj (theta only) via sigma:
        if (identical(param, "theta")) {
          se_j_local <- if (is.null(sig_eps)) 1 else sig_eps[jj]
          dz_sigma   <- -(zmat / se_j_local)                       # N x (m+2)
          A <- phi_inner * dz_sigma[, -1, drop = FALSE]            # align z_k columns
          C <- phi_prev  * dz_sigma[, -ncol(dz_sigma), drop = FALSE] # align z_{k-1}
          d_sigma <- as.numeric((A[, 1:m, drop = FALSE] - C[, 1:m, drop = FALSE]) %*% matrix(scores_j, ncol = 1)) # N-vector
          d_var   <- d_sigma * (1 / (2 * se_j_local))              # chain sigma->variance
        }

        # --- Pack gradients into G (N x P_item) in the same order as keep -----
        th_names <- paste0(j, "|t", seq_len(m - 1))
        nu_name  <- { nm <- paste0("nu[", j, "]"); if (!nm %in% coef_names) { alt <- paste0(j, "~1"); if (alt %in% coef_names) nm <- alt else nm <- NA_character_ }; nm }
        lam_names <- paste0(j, "=~", lam_cols)
        var_name  <- if (identical(param, "theta")) paste0(j, "~~", j) else NA_character_

        keep <- c(th_names, nu_name, lam_names, var_name)
        keep <- keep[!is.na(keep) & keep %in% coef_names]
        idx  <- match(keep, coef_names)
        Vsub <- V_full[idx, idx, drop = FALSE]

        # Build G:
        # thresholds block: d_tau_mat (N x (m-1)) but only keep those present
        G_parts <- list()
        if (length(th_names)) {
          sel <- th_names %in% keep
          if (any(sel)) G_parts[["tau"]] <- d_tau_mat[, sel, drop = FALSE]
        }
        if (!is.na(nu_name) && nu_name %in% keep) {
          G_parts[["nu"]] <- matrix(d_nu, ncol = 1)
          colnames(G_parts[["nu"]]) <- nu_name
        }
        # loadings: keep only existing (lavaan may drop zeros)
        lam_keep <- lam_names[lam_names %in% keep]
        if (length(lam_keep)) {
          # d_lam_mat columns in lam_cols order; select those present
          sel <- lam_cols %in% sub("^.*=~", "", lam_keep)
          G_parts[["lam"]] <- d_lam_mat[, sel, drop = FALSE]
          colnames(G_parts[["lam"]]) <- lam_keep
        }
        if (identical(param, "theta") && !is.na(var_name) && var_name %in% keep) {
          G_parts[["var"]] <- matrix(d_var, ncol = 1)
          colnames(G_parts[["var"]]) <- var_name
        }

        if (length(G_parts)) {
          G <- do.call(cbind, G_parts)  # N x P_item, columns in the same order as `keep`
          # SE_i = sqrt( diag( G Vsub G^T ) )  -> efficient via rowSums((G %*% Vsub) * G)
          GV <- G %*% Vsub
          se <- sqrt(rowSums(GV * G))
        } else {
          se <- rep(NA_real_, N_g)
        }

        # CI vectors
        yexp_lwr_j <- yexp_j - zcrit * se
        yexp_upr_j <- yexp_j + zcrit * se

        # Interleaved per-item tibble: .yexp_j, .yexp_lwr_j, .yexp_upr_j
        col_core <- setNames(tibble::tibble(
          !!paste0(prefix_yhat, j)     := yexp_j,
          !!paste0(prefix_ci[[1]], j)  := yexp_lwr_j,
          !!paste0(prefix_ci[[2]], j)  := yexp_upr_j
        ), c(paste0(prefix_yhat, j), paste0(prefix_ci[[1]], j), paste0(prefix_ci[[2]], j)))

      } else {
        # No CI -> only expected score column
        col_core <- setNames(tibble::tibble(!!paste0(prefix_yhat, j) := yexp_j),
                             paste0(prefix_yhat, j))
      }

      yexp_ci_cols[[length(yexp_ci_cols) + 1L]] <- col_core
    } # end items

    ystar_list[[g]] <- ystar_df
    probs_list[[g]] <- if (is.null(probs_g)) tibble::tibble() else tibble::as_tibble(probs_g)
    # bind interleaved yexp (+/- CI) columns for this group
    yexp_ci_list[[g]] <- if (length(yexp_ci_cols)) dplyr::bind_cols(yexp_ci_cols) else tibble::tibble()
  } # end groups

  # -- Final bind: keep order as FS | ystar | [ yexp_j, lwr_j, upr_j ... ] | probs
  fs_tbl_parts <- vector("list", n_groups)
  for (g in seq_len(n_groups)) {
    fs_g <- tibble::as_tibble(fs_list[[g]])
    fs_g <- tibble::add_column(fs_g, .gid = g, .before = 1)
    fs_g <- tibble::add_column(fs_g, group = group_labels[g], .after = ".gid")
    fs_tbl_parts[[g]] <- fs_g
  }
  fs_tbl    <- dplyr::bind_rows(fs_tbl_parts)
  ystar_tbl <- dplyr::bind_rows(ystar_list)
  yexp_tbl  <- dplyr::bind_rows(yexp_ci_list)
  probs_tbl <- dplyr::bind_rows(probs_list)

  out <- dplyr::bind_cols(fs_tbl, ystar_tbl, yexp_tbl, probs_tbl)

  # Coerce special classes -> double
  .needs_double_from_special <- function(x) inherits(x, "lvn.mtrx") || inherits(x, "lavaan.matrix")
  out <- dplyr::mutate(out, dplyr::across(.cols = where(.needs_double_from_special), .fns = ~ as.double(.x)))

  tibble::as_tibble(out)
}
