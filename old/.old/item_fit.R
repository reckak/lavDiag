#' Prepare original data + model-based and empirical item curves (continuous & ordinal)
#'
#' Returns `.rid`, `.gid`, `.group`, original indicators, factor scores
#' (including SEs for continuous-only models), and two kinds of item curves:
#' (i) model-based expected scores + CIs from `augment()` with prefixes `m_*`,
#' and (ii) empirical expected-score curves + CIs from `mgcv::gam()` with prefixes `e_*`.
#' Additionally returns a per-item/per-group metrics tibble that quantifies the
#' *agreement between the model-based and empirical curves* via R^2, RMSE, MAE,
#' and their complexity-penalized variants.
#'
#' @param fit        Fitted `lavaan` object; must be converged and contain latent variables.
#' @param data       Optional `lavPredict` data/list; computed if `NULL`.
#' @param info       Optional `model_info(fit)`; computed if `NULL`.
#' @param level      Confidence level for empirical CIs (default `0.95`).
#' @param fam_cont   `mgcv` family for continuous items (default `stats::gaussian()`).
#' @param fam_ord    `mgcv` family for ordinal items (default `mgcv::betar(link = "logit")`).
#' @param gam_args_cont List of extra args for `mgcv::gam()` **for continuous items** (formula/data/family locked).
#' @param gam_args_ord  List of extra args for `mgcv::gam()` **for ordinal items** (formula/data/family locked).
#' @param plan       Parallel plan: one of `c("auto","multisession","multicore","sequential","cluster","none")`.
#' @param workers    Integer number of workers (defaults inside `.set_future_plan()` when `NULL`).
#' @param cluster    Optional cluster object/host spec if `plan = "cluster"`.
#' @param progress   Logical; show `furrr` progress bar (default `FALSE`).
#'
#' @details
#' Empirical curves are fit per item and per group using `mgcv::gam()` on the link scale,
#' then transformed back using the family's inverse link. For ordinal outcomes the response
#' is rescaled to (0,1) with the Smithson–Verkuilen adjustment and modeled with `betar`;
#' predictions are clamped in (0,1) on the probability scale and then mapped back to the
#' original item scale.
#'
#' Parallelization: groups are processed sequentially while items within a group are
#' processed in parallel via `furrr::future_map()`. The helper `.set_future_plan()` is used
#' to safely set and restore `future::plan()`.
#'
#' Factor loadings: relevant latent predictors for each item can vary by group. The function
#' first tries group-specific loadings; if missing, it falls back to overall loadings and then
#' to all available factor-score columns. For speed, if loadings do not vary across groups,
#' group-specific lookup is skipped.
#'
#' Metrics: For each item and group, the function computes
#'   - `R^2` as Pearson `cor(m_est_yhat_*, e_est_yhat_*)^2`,
#'   - `RMSE` and `MAE` between the two curve estimates on the original item scale.
#' Penalized variants account for effective complexity with
#' **`k_eff = c_m + sum(edf)`**, where `c_m` is model-side complexity
#'   (number of relevant latents + intercept for continuous items, or number of thresholds
#'   for ordinal items) and `sum(edf)` is the sum of smoothing edf from the fitted GAM.
#' Formulas:
#' `R2_pen = 1 - (1 - R2) * (n_eff - 1) / pmax(1, n_eff - k_eff - 1)`,
#' `RMSE_pen = sqrt(MSE * n_eff / pmax(1, n_eff - k_eff))`,
#' `MAE_pen = MAE * sqrt(n_eff / pmax(1, n_eff - k_eff))`.
#'
#' @return A list with two tibbles:
#'   - `$data`: original augmented data with `m_*` and `e_*` columns,
#'   - `$metrics`: per-item/per-group tibble with columns
#'     `.gid`, `.group`, `item`, `type` ("continuous"/"ordinal"),
#'     `n_eff`, `c_m`, `c_e`, `k_eff`, `r2`, `rmse`, `mae`, `r2_pen`, `rmse_pen`, `mae_pen`.
#'
#' @export
item_fit_data_old <- function(fit,
                          data     = NULL,
                          info     = NULL,
                          level    = 0.95,
                          fam_cont = stats::gaussian(),
                          fam_ord  = mgcv::betar(link = "logit"),
                          gam_args_cont = list(method = "REML"),
                          gam_args_ord  = list(method = "REML"),
                          plan     = c("auto","multisession","multicore",
                                       "sequential","cluster","none"),
                          workers  = NULL,
                          cluster  = NULL,
                          progress = FALSE) {

  plan <- match.arg(plan)

  # -- Assertions --------------------------------------------------------------
  .assert_lavaan_fit(fit, require_converged = TRUE, require_latent = TRUE, forbid_multilevel = TRUE)
  if (is.null(info)) info <- model_info(fit)

  # Ensure mgcv (and betar) availability (betar introduced in mgcv 1.8-31)
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is required for empirical curves.")
  }
  if (isTRUE(inherits(fam_ord, "family")) && identical(fam_ord$family, "Beta regression")) {
    # user supplied their own betar(), skip version check
  } else {
    # check version when using default betar()
    vmgcv <- utils::packageVersion("mgcv")
    if (vmgcv < "1.8.31") {
      stop("mgcv >= 1.8-31 is required for mgcv::betar(). Installed: ", as.character(vmgcv))
    }
  }

  ov_all   <- info$observed_variables
  ov_cont  <- info$ov_continuous
  ov_ord   <- info$ov_ordinal
  lats     <- info$latent_variables
  n_groups <- info$n_groups %||% 1L
  gvar     <- info$group_var %||% ".group"

  # -- Utility: safe builder for mgcv::gam calls -------------------------------
  build_gam_call <- function(fml, dat, fam, user_args) {
    # Lock fields we control
    user_args$formula <- NULL
    user_args$data    <- NULL
    user_args$family  <- NULL
    args <- c(list(formula = fml, data = dat, family = fam), user_args)
    do.call(mgcv::gam, args)
  }

  # -- 1) Model-based curves via augment() ------------------------------------
  out <- augment(
    fit            = fit,
    data           = data,
    info           = info,
    yhat           = TRUE, ci = TRUE, resid = FALSE,
    se_yhat        = FALSE, ystar = FALSE, pr = FALSE,
    se_fs          = TRUE,
    prefix_yhat    = "m_est_yhat_",
    prefix_ci      = c("m_lwr_yhat_", "m_upr_yhat_"),
    col_layout     = "by_type"
  )

  # -- Helpers -----------------------------------------------------------------
  rescale_01 <- function(x) {
    r <- range(x, na.rm = TRUE, finite = TRUE)
    if (!is.finite(r[1]) || !is.finite(r[2]) || r[2] == r[1]) return(rep(0.5, length(x)))
    (x - r[1]) / (r[2] - r[1])
  }
  sv01 <- function(y01) {
    n <- sum(!is.na(y01))
    if (n <= 1L) return(y01)
    ((y01 * (n - 1)) + 0.5) / n
  }
  linkinv_fun <- function(fam) {
    if (!is.null(fam$linkinv)) fam$linkinv else stop("Family lacks linkinv.")
  }
  clamp01 <- function(p, eps = 1e-6) pmin(pmax(p, eps), 1 - eps)

  # -- Loadings & thresholds (including group-specific) -----------------------
  PE_raw <- lavaan::parameterEstimates(fit, standardized = FALSE)
  # loadings
  PE_sel <- PE_raw[PE_raw$op == "=~" & PE_raw$rhs %in% ov_all, , drop = FALSE]
  keep_cols <- intersect(c("group","lhs","rhs","est"), names(PE_sel))
  PE_L <- PE_sel[, keep_cols, drop = FALSE]
  if (!"group" %in% names(PE_L)) PE_L$group <- 1L

  tol <- 1e-10
  loads_by_item_overall <- tapply(seq_len(nrow(PE_L)), PE_L$rhs, function(idx) {
    lhs <- PE_L$lhs[idx]; est <- PE_L$est[idx]
    unique(lhs[is.finite(est) & abs(est) > tol])
  })
  loads_by_item_overall <- loads_by_item_overall[ov_all]

  PE_by_group <- split(PE_L, PE_L$group)
  loads_by_item_by_group <- lapply(PE_by_group, function(PEg) {
    if (!nrow(PEg)) return(setNames(vector("list", length(ov_all)), ov_all))
    li <- tapply(seq_len(nrow(PEg)), PEg$rhs, function(idx) {
      lhs <- PEg$lhs[idx]; est <- PEg$est[idx]
      unique(lhs[is.finite(est) & abs(est) > tol])
    })
    li[ov_all[!ov_all %in% names(li)]] <- list(NULL)
    li[ov_all]
  })

  # thresholds (ordinal complexity component): count number of FREE thresholds per item/group
  PE_thr <- PE_raw[PE_raw$op == "|" & PE_raw$lhs %in% ov_all, , drop = FALSE]
  if (!"group" %in% names(PE_thr)) PE_thr$group <- 1L

  thr_counts_by_group <- lapply(split(PE_thr, PE_thr$group), function(PEg) {
    if (!nrow(PEg)) {
      return(list(
        all  = setNames(integer(length(ov_all)), ov_all),
        free = setNames(integer(length(ov_all)), ov_all)
      ))
    }
    tabs_all <- tapply(PEg$lhs, PEg$lhs, length)
    tabs_free <- tapply(seq_len(nrow(PEg)), PEg$lhs, function(idx) {
      f <- PEg$free[idx]; length(unique(f[f > 0]))
    })

    all_vec  <- setNames(integer(length(ov_all)), ov_all)
    free_vec <- setNames(integer(length(ov_all)), ov_all)
    all_vec[names(tabs_all)]   <- as.integer(tabs_all)
    free_vec[names(tabs_free)] <- as.integer(tabs_free)

    list(all = all_vec, free = free_vec)
  })

  # detect if loadings vary by group
  set_equal <- function(a, b) {
    if (length(a) != length(b)) return(FALSE)
    all(sort(a) == sort(b))
  }
  vary_by_group <- any(vapply(ov_all, function(j) {
    overall <- loads_by_item_overall[[j]]
    if (is.null(overall)) return(TRUE)  # conservative
    any(!vapply(loads_by_item_by_group, function(li) {
      set_equal(li[[j]] %||% character(0L), overall %||% character(0L))
    }, logical(1L)))
  }, logical(1L)))

  # -- mgcv setup --------------------------------------------------------------
  zcrit <- stats::qnorm(1 - (1 - level) / 2)
  fam_for_item <- function(j) if (j %in% ov_cont) fam_cont else fam_ord
  inv_for_item <- function(j) linkinv_fun(fam_for_item(j))

  # original-scale ranges for back-transform of ordinal expected scores
  rng_by_item <- lapply(ov_all, function(j) {
    if (!j %in% names(out)) return(c(0,1))
    x <- out[[j]]
    if (is.factor(x)) c(1L, nlevels(x)) else range(as.numeric(x), na.rm = TRUE, finite = TRUE)
  })
  names(rng_by_item) <- ov_all

  # -- Parallel plan -----------------------------------------------------------
  reset_plan <- .set_future_plan(plan = plan, workers = workers, cluster = cluster)
  on.exit(reset_plan(), add = TRUE)

  # storage for empirical complexity (sum edf) per item/group
  ce_map <- vector("list", n_groups)
  names(ce_map) <- as.character(seq_len(n_groups))

  # -- 2) Empirical GAM curves (parallel over items within each group) --------
  for (g in seq_len(n_groups)) {
    rows_g <- if (".gid" %in% names(out)) which(out$.gid == g) else seq_len(nrow(out))
    if (!length(rows_g)) next

    eta_cols <- intersect(lats, names(out))
    if (!length(eta_cols)) stop("Factor scores missing in data; cannot build empirical curves.")

    out_g <- out[rows_g, c(ov_all, eta_cols), drop = FALSE]
    fml_cache <- new.env(parent = emptyenv())      # cache formulas per unique predictor set
    newdata_cache <- new.env(parent = emptyenv())  # cache newdata per unique predictor set

    item_results <- furrr::future_map(
      ov_all,
      function(j) {
        # choose relevant latents
        rel_lats <- if (vary_by_group) {
          li_g <- loads_by_item_by_group[[as.character(g)]] %||% loads_by_item_by_group[[g]]
          (if (!is.null(li_g)) li_g[[j]] else NULL) %||% loads_by_item_overall[[j]]
        } else {
          loads_by_item_overall[[j]]
        }
        rel_lats <- intersect(rel_lats %||% eta_cols, eta_cols)
        if (!length(rel_lats) || !j %in% names(out_g)) return(NULL)

        df_g <- out_g[, c(j, rel_lats), drop = FALSE]
        fam_j   <- fam_for_item(j)
        invlink <- inv_for_item(j)
        rj      <- rng_by_item[[j]]

        # build formula from cache
        key <- paste(rel_lats, collapse = "|")
        fml <- get0(key, envir = fml_cache, inherits = FALSE)
        if (is.null(fml)) {
          sm_terms <- paste(sprintf("s(%s)", rel_lats), collapse = " + ")
          fml <- stats::as.formula(paste("y_resp ~", sm_terms))
          assign(key, fml, envir = fml_cache)
        }

        ok    <- stats::complete.cases(df_g[, c(j, rel_lats), drop = FALSE])
        df_ok <- df_g[ok, , drop = FALSE]
        if (!nrow(df_ok)) {
          n <- nrow(df_g)
          return(list(
            j   = j,
            mu  = rep(NA_real_, n),
            lwr = rep(NA_real_, n),
            upr = rep(NA_real_, n),
            ce  = NA_real_
          ))
        }

        # outcome + back-transformers
        if (j %in% ov_cont) {
          y_fit <- df_ok[[j]]
          back_to_original <- function(mu_link) invlink(mu_link)
          back_ci <- function(mu_link) invlink(mu_link)
        } else {
          y_fit <- sv01(rescale_01(df_ok[[j]]))
          to_original_scale <- function(p01) { d <- rj[2] - rj[1]; p01 * d + rj[1] }
          back_to_original  <- function(mu_link) to_original_scale(clamp01(invlink(mu_link)))
          back_ci           <- function(mu_link) to_original_scale(clamp01(invlink(mu_link)))
        }
        X_fit <- df_ok[, rel_lats, drop = FALSE]

        if (nrow(X_fit) >= length(rel_lats) + 2L) {
          user_args <- if (j %in% ov_cont) gam_args_cont else gam_args_ord

          fit_gam <- build_gam_call(
            fml,
            dat = cbind(y_resp = y_fit, X_fit),
            fam = fam_j,
            user_args = user_args
          )

          # newdata from cache (same row order as df_g)
          nd <- get0(key, envir = newdata_cache, inherits = FALSE)
          if (is.null(nd)) {
            nd <- df_g[, rel_lats, drop = FALSE]
            assign(key, nd, envir = newdata_cache)
          }

          pr <- predict(
            fit_gam,
            newdata = nd,
            type = "link",
            se.fit = TRUE
          )

          mu  <- back_to_original(pr$fit)
          lwr <- back_ci(pr$fit - zcrit * pr$se.fit)
          upr <- back_ci(pr$fit + zcrit * pr$se.fit)

          # empirical complexity: sum of edf (includes intercept if present — that's fine)
          ce  <- sum(fit_gam$edf, na.rm = TRUE)
        } else {
          n <- nrow(df_g)
          mu <- lwr <- upr <- rep(NA_real_, n)
          ce <- NA_real_
        }

        list(j = j, mu = as.numeric(mu), lwr = as.numeric(lwr), upr = as.numeric(upr), ce = ce)
      },
      .options  = furrr::furrr_options(
        seed     = TRUE,
        packages = c("mgcv","stats")
      ),
      .progress = progress
    )

    # write back & store empirical complexity
    ce_vec <- setNames(numeric(length(ov_all)), ov_all)
    for (res in item_results) {
      if (is.null(res)) next
      nm_est <- paste0("e_est_yhat_", res$j)
      nm_lwr <- paste0("e_lwr_yhat_", res$j)
      nm_upr <- paste0("e_upr_yhat_", res$j)
      out[rows_g, nm_est] <- res$mu
      out[rows_g, nm_lwr] <- res$lwr
      out[rows_g, nm_upr] <- res$upr
      ce_vec[[res$j]] <- res$ce
    }
    ce_map[[as.character(g)]] <- ce_vec
  }

  # -- Numeric types cleanup ---------------------------------------------------
  is_num <- vapply(out, is.numeric, logical(1L))
  keep_int <- names(out) %in% c(".rid", ".gid")
  to_double <- is_num & !keep_int
  for (nm in names(out)[to_double]) out[[nm]] <- as.double(out[[nm]])

  # -- 3) Metrics per item & group --------------------------------------------
  metrics <- list()
  for (g in seq_len(n_groups)) {
    rows_g <- if (".gid" %in% names(out)) which(out$.gid == g) else seq_len(nrow(out))

    g_label <- if (n_groups > 1L) {
      if (gvar %in% names(out)) {
        gvals <- out[[gvar]][rows_g]
        idx   <- which(!is.na(gvals) & nzchar(as.character(gvals)))
        if (length(idx)) as.character(gvals[idx[1]]) else as.character(g)
      } else {
        if (!is.null(info$group_labels) && length(info$group_labels) >= g) {
          as.character(info$group_labels[g])
        } else {
          as.character(g)
        }
      }
    } else {
      if (!is.null(info$group_labels) && length(info$group_labels) >= 1L) {
        as.character(info$group_labels[1L])
      } else {
        "1"
      }
    }

    ce_vec  <- ce_map[[as.character(g)]]

    for (j in ov_all) {
      mcol <- paste0("m_est_yhat_", j)
      ecol <- paste0("e_est_yhat_", j)
      if (!all(c(mcol, ecol) %in% names(out))) next

      x <- out[rows_g, mcol, drop = TRUE]
      y <- out[rows_g, ecol, drop = TRUE]
      ok <- is.finite(x) & is.finite(y)
      n_eff <- sum(ok)
      if (n_eff < 2L || (sd(x[ok]) == 0 && sd(y[ok]) == 0)) {
        r2 <- rmse <- mae <- NA_real_
      } else {
        r  <- suppressWarnings(stats::cor(x[ok], y[ok]))
        r2 <- if (is.finite(r)) r^2 else NA_real_
        dif <- y[ok] - x[ok]
        rmse <- sqrt(mean(dif^2))
        mae  <- mean(abs(dif))
      }

      # model-side complexity c_m
      rel_lats <- if (vary_by_group) {
        li_g <- loads_by_item_by_group[[as.character(g)]] %||% loads_by_item_by_group[[g]]
        (if (!is.null(li_g)) li_g[[j]] else NULL) %||% loads_by_item_overall[[j]]
      } else {
        loads_by_item_overall[[j]]
      }
      rel_lats <- rel_lats %||% character(0L)
      n_rel <- length(rel_lats)

      if (j %in% ov_cont) {
        cm <- n_rel + 1L  # latent slopes + intercept (nu_j)
        type_j <- "continuous"
      } else {
        thr_gc <- thr_counts_by_group[[as.character(g)]]
        n_thr_all  <- if (!is.null(thr_gc)) as.integer(thr_gc$all[[j]]  %||% 0L) else 0L
        n_thr_free <- if (!is.null(thr_gc)) as.integer(thr_gc$free[[j]] %||% 0L) else 0L
        n_thr <- if (is.finite(n_thr_all) && n_thr_all > 0L) n_thr_all else n_thr_free
        if (!is.finite(n_thr) || n_thr <= 0L) {
          xj <- out[[j]]
          if (is.factor(xj)) n_thr <- max(0L, nlevels(xj) - 1L) else n_thr <- 0L
        }
        cm <- n_rel + n_thr
        type_j <- "ordinal"
      }

      ce <- ce_vec[[j]]
      # Consistent complexity: k_eff = c_m + sum(edf)
      if (is.finite(ce)) {
        k_eff <- cm + ce
      } else {
        k_eff <- cm
      }

      # Penalized variants (agreement-adjusted for complexity)
      if (is.finite(r2) && n_eff > (k_eff + 1)) {
        r2_pen <- 1 - (1 - r2) * (n_eff - 1) / max(1, n_eff - k_eff - 1)
      } else {
        r2_pen <- NA_real_
      }
      if (is.finite(rmse) && n_eff > k_eff) {
        mse <- rmse^2
        rmse_pen <- sqrt(mse * n_eff / max(1, n_eff - k_eff))
      } else {
        rmse_pen <- NA_real_
      }
      if (is.finite(mae) && n_eff > k_eff) {
        mae_pen <- mae * sqrt(n_eff / max(1, n_eff - k_eff))  # heuristic scaling
      } else {
        mae_pen <- NA_real_
      }

      metrics[[length(metrics) + 1L]] <- list(
        .gid      = g,
        .group    = g_label,
        item      = j,
        type      = type_j,
        n_eff     = n_eff,
        c_m       = as.numeric(cm),
        c_e       = as.numeric(if (is.finite(ce)) ce else NA_real_),
        k_eff     = as.numeric(k_eff),
        r2        = r2,
        rmse      = rmse,
        mae       = mae,
        r2_pen    = r2_pen,
        rmse_pen  = rmse_pen,
        mae_pen   = mae_pen
      )
    }
  }
  metrics_tbl <- if (length(metrics)) tibble::as_tibble(do.call(rbind, lapply(metrics, tibble::as_tibble_row)))
  else tibble::tibble(.gid = integer(),
                      .group = character(),
                      item = character(),
                      type = character(),
                      n_eff = integer(),
                      c_m = double(),
                      c_e = double(),
                      k_eff = double(),
                      r2 = double(),
                      rmse = double(),
                      mae = double(),
                      r2_pen = double(),
                      rmse_pen = double(),
                      mae_pen = double())

  # -- Return both data and metrics -------------------------------------------
  list(
    data    = tibble::as_tibble(out),
    metrics = metrics_tbl
  )
}
