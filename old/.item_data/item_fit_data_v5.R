item_fit_data_old <- function(fit,
                          data     = NULL,
                          info     = NULL,
                          level    = 0.95,
                          fam_cont = stats::gaussian(),
                          fam_ord  = "ocat",
                          gam_args_cont = list(method = "REML"),
                          gam_args_ord  = list(method = "REML"),
                          plan     = c("auto","multisession","multicore",
                                       "sequential","cluster","none"),
                          workers  = NULL,
                          cluster  = NULL,
                          progress = FALSE,
                          verbose  = TRUE) {

  plan <- match.arg(plan)

  # Verbose helper (prints only when verbose = TRUE)
  vmsg <- function(...) { if (isTRUE(verbose)) message(...) }

  # -- External helpers presence check ---------------------------------------
  required_helpers <- c(
    ".fit_gam_models", ".predict_gam_cols",
    ".rescale_01", ".sv01", ".clamp01", ".linkinv_fun"
  )
  # If user asked for ocat, also require the delta-CI predictor helper
  is_ocat <- (is.character(fam_ord) && identical(fam_ord, "ocat")) || (isTRUE(inherits(fam_ord, "family")) && identical(fam_ord$family, "ocat"))
  if (is_ocat) required_helpers <- c(required_helpers, ".predict_ocat_ci_delta")

  missing_helpers <- required_helpers[!vapply(required_helpers, function(nm) exists(nm, mode = "function", inherits = TRUE), logical(1))]
  if (length(missing_helpers)) {
    stop("Missing helper functions: ", paste(missing_helpers, collapse = ", "),
         ". Please source/load them before calling item_fit_data().")
  }
  # Bind helpers to local symbols so furrr/future serializes them explicitly to workers
  fit_fun    <- .fit_gam_models
  pred_fun   <- .predict_gam_cols
  rescale_fn <- .rescale_01
  sv_fn      <- .sv01
  clamp_fn   <- .clamp01
  if (is_ocat) pred_ocat <- .predict_ocat_ci_delta

  vmsg("Starting item_fit_data()...")

  # -- Assertions --------------------------------------------------------------
  .assert_lavaan_fit(fit, require_converged = TRUE, require_latent = TRUE, forbid_multilevel = TRUE)
  if (is.null(info)) info <- model_info(fit)

  # Ensure mgcv (and betar/ocat) availability (betar introduced in mgcv 1.8-31)
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is required for empirical curves.")
  }
  # If user explicitly passes betar(), check mgcv version
  if (isTRUE(inherits(fam_ord, "family")) && identical(fam_ord$family, "Beta regression")) {
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

  # -- 1) Model-based curves via augment() ------------------------------------
  vmsg("Computing model-based curves (augment)...")
  out <- augment(
    fit            = fit,
    data           = data,
    info           = info,
    yhat = TRUE, ci = TRUE, resid = TRUE,
    se_yhat        = FALSE, ystar = FALSE, pr = FALSE,
    se_fs          = TRUE,
    prefix_yhat    = "m_est_yhat_",
    prefix_ci      = c("m_lwr_yhat_", "m_upr_yhat_"),
    col_layout     = "by_type"
  )
  vmsg("...augment done.")
  # -- Loadings & thresholds (including group-specific) -----------------------
  vmsg("Extracting loadings/thresholds & preparing per-group maps...")
  PE_raw <- lavaan::parameterEstimates(fit, standardized = FALSE)

  # loadings
  PE_sel <- PE_raw[PE_raw$op == "=~" & PE_raw$rhs %in% ov_all, , drop = FALSE]
  keep_cols <- intersect(c("group","lhs","rhs","est"), names(PE_sel))
  PE_L <- PE_sel[, keep_cols, drop = FALSE]
  if (!"group" %in% names(PE_L)) PE_L$group <- rep_len(1L, nrow(PE_L))

  tol <- 1e-10
  if (nrow(PE_L) == 0L) {
    loads_by_item_overall <- setNames(vector("list", length(ov_all)), ov_all)
  } else {
    loads_by_item_overall <- tapply(seq_len(nrow(PE_L)), PE_L$rhs, function(idx) {
      lhs <- PE_L$lhs[idx]; est <- PE_L$est[idx]
      unique(lhs[is.finite(est) & abs(est) > tol])
    })
    loads_by_item_overall <- loads_by_item_overall[ov_all]
  }

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

  # thresholds (ordinal complexity component)
  PE_thr <- PE_raw[PE_raw$op == "|" & PE_raw$lhs %in% ov_all, , drop = FALSE]
  if (!"group" %in% names(PE_thr)) PE_thr$group <- rep_len(1L, nrow(PE_thr))

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
  vmsg("...loadings/thresholds processed.")

  # -- mgcv setup --------------------------------------------------------------
  zcrit <- stats::qnorm(1 - (1 - level) / 2)
  fam_for_item <- function(j) {
    if (j %in% ov_cont) return(fam_cont)
    # ordinal branch
    return(fam_ord)
  }
  # original-scale ranges for back-transform of ordinal expected scores
  rng_by_item <- lapply(ov_all, function(j) {
    if (!j %in% names(out)) return(c(0,1))
    x <- out[[j]]
    if (is.factor(x)) c(1L, nlevels(x)) else range(as.numeric(x), na.rm = TRUE, finite = TRUE)
  })
  names(rng_by_item) <- ov_all

  # Build per-item value maps for ordinal items (exact observed codes per category)
  ord_value_map <- lapply(ov_all, function(j) {
    x <- out[[j]]
    if (is.null(x)) return(NULL)
    if (is.factor(x) || is.ordered(x)) {
      levs <- levels(x)
      suppressWarnings({ lev_num <- as.numeric(levs) })
      if (all(is.finite(lev_num))) lev_num else seq_len(nlevels(x))
    } else if (j %in% ov_ord) {
      sort(unique(as.numeric(x)))
    } else {
      NULL
    }
  })
  names(ord_value_map) <- ov_all

  # Helper: expected value on original (possibly non-equidistant) codes
  .emap_expect <- function(E_index, vals) {
    R <- length(vals)
    lo <- pmax(1, pmin(R - 1, floor(E_index)))
    hi <- lo + 1
    w  <- E_index - lo
    vals[lo] + w * (vals[hi] - vals[lo])
  }

  # -- Parallel plan -----------------------------------------------------------
  vmsg(sprintf("Fitting empirical GAM curves (%d group%s)...", n_groups, ifelse(n_groups > 1, "s", "")))
  # Store and restore future plan
  reset_plan <- .set_future_plan(plan = plan, workers = workers, cluster = cluster)
  on.exit(reset_plan(), add = TRUE)

  # storage for empirical complexity (sum edf) per item/group
  ce_map   <- vector("list", n_groups)
  gam_fits <- vector("list", n_groups)
  names(ce_map) <- as.character(seq_len(n_groups))
  names(gam_fits) <- as.character(seq_len(n_groups))

  # -- 2) Empirical GAM curves (parallel over items within each group) --------
  for (g in seq_len(n_groups)) {
    vmsg(sprintf("  Group %d: fitting %d item GAMs...", g, length(ov_all)))
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

        # Two branches: continuous, ordinal-betar, or ordinal-ocat
        if (j %in% ov_cont) {
          # Continuous -> gaussian/beta-like on original scale
          y_fit <- df_ok[[j]]
          X_fit <- df_ok[, rel_lats, drop = FALSE]

          if (nrow(X_fit) >= length(rel_lats) + 2L) {
            user_args <- gam_args_cont
            fit_gam <- fit_fun(
              data    = cbind(y_resp = y_fit, X_fit),
              formula = fml,
              family  = fam_j,
              args    = user_args
            )

            nd <- get0(key, envir = newdata_cache, inherits = FALSE)
            if (is.null(nd)) { nd <- df_g[, rel_lats, drop = FALSE]; assign(key, nd, envir = newdata_cache) }

            pred_df <- pred_fun(dat = nd, fit = fit_gam,
                                prefixes = c(".tmp_est_", ".tmp_lwr_", ".tmp_upr_"),
                                suffix   = "", level = level)
            mu  <- pred_df[[1L]]; lwr <- pred_df[[2L]]; upr <- pred_df[[3L]]
            ce  <- sum(fit_gam$edf, na.rm = TRUE)
          } else {
            n <- nrow(df_g); mu <- lwr <- upr <- rep(NA_real_, n); ce <- NA_real_
          }
        } else if (is_ocat) {
          # Ordinal via mgcv::ocat with ordered outcome, then expected score using delta CI
          # Ensure ordered coding 1..R
          y_ord <- df_ok[[j]]
          if (!is.ordered(y_ord)) y_ord <- ordered(y_ord)
          # mgcv::ocat() expects INTEGER class labels 1..R, not a factor
          y_fit <- as.integer(y_ord)
          X_fit <- df_ok[, rel_lats, drop = FALSE]

          if (nrow(X_fit) >= length(rel_lats) + 2L) {
            # Ensure mgcv::ocat knows the number of categories (R)
            fam_use <- fam_j
            if ((is.character(fam_j) && identical(fam_j, "ocat")) ||
                (isTRUE(inherits(fam_j, "family")) && identical(fam_j$family, "ocat"))) {
              Rj <- nlevels(y_ord)
              if (!is.finite(Rj) || Rj < 2L) stop(sprintf("ocat needs at least 2 categories for item '%s'", j))
              fam_use <- mgcv::ocat(link = "identity", R = Rj)
            }
            fit_gam <- fit_fun(
              data    = cbind(y_resp = y_fit, X_fit),
              formula = fml,
              family  = fam_use,
              args    = gam_args_ord
            )

            nd <- get0(key, envir = newdata_cache, inherits = FALSE)
            if (is.null(nd)) { nd <- df_g[, rel_lats, drop = FALSE]; assign(key, nd, envir = newdata_cache) }

            oc <- pred_ocat(fit_gam, newdata = nd, level = level)
            pr_mat <- try(suppressWarnings(predict(fit_gam, newdata = nd, type = "response")), silent = TRUE)
            vals <- ord_value_map[[j]]
            if (is.matrix(pr_mat) && !is.null(vals) && ncol(pr_mat) == length(vals)) {
              mu  <- as.numeric(pr_mat %*% vals)
              lwr <- as.numeric(.emap_expect(oc$E_lwr, vals))
              upr <- as.numeric(.emap_expect(oc$E_upr, vals))
            } else {
              if (is.null(vals)) {
                to_original_from_scores <- function(E) {
                  R <- length(fit_gam$family$getTheta(trans = TRUE)) + 1L
                  if (!is.finite(rj[1]) || !is.finite(rj[2]) || R <= 1L) return(E)
                  rj[1] + (E - 1) * ((rj[2] - rj[1]) / (R - 1))
                }
                mu  <- to_original_from_scores(oc$E_est)
                lwr <- to_original_from_scores(oc$E_lwr)
                upr <- to_original_from_scores(oc$E_upr)
              } else {
                mu  <- .emap_expect(oc$E_est, vals)
                lwr <- .emap_expect(oc$E_lwr, vals)
                upr <- .emap_expect(oc$E_upr, vals)
              }
            }
            ce  <- sum(fit_gam$edf, na.rm = TRUE)
          } else {
            n <- nrow(df_g); mu <- lwr <- upr <- rep(NA_real_, n); ce <- NA_real_
          }
        } else {
          # Ordinal via Beta regression on scaled [0,1]
          y_fit <- sv_fn(rescale_fn(df_ok[[j]]))
          X_fit <- df_ok[, rel_lats, drop = FALSE]

          if (nrow(X_fit) >= length(rel_lats) + 2L) {
            fit_gam <- fit_fun(
              data    = cbind(y_resp = y_fit, X_fit),
              formula = fml,
              family  = fam_j,
              args    = gam_args_ord
            )

            nd <- get0(key, envir = newdata_cache, inherits = FALSE)
            if (is.null(nd)) { nd <- df_g[, rel_lats, drop = FALSE]; assign(key, nd, envir = newdata_cache) }

            pred_df <- pred_fun(dat = nd, fit = fit_gam,
                                prefixes = c(".tmp_est_", ".tmp_lwr_", ".tmp_upr_"),
                                suffix   = "", level = level)
            to_original_scale <- function(p01) { d <- rj[2] - rj[1]; p01 * d + rj[1] }
            mu  <- to_original_scale(clamp_fn(pred_df[[1L]]))
            lwr <- to_original_scale(clamp_fn(pred_df[[2L]]))
            upr <- to_original_scale(clamp_fn(pred_df[[3L]]))
            ce  <- sum(fit_gam$edf, na.rm = TRUE)
          } else {
            n <- nrow(df_g); mu <- lwr <- upr <- rep(NA_real_, n); ce <- NA_real_
          }
        }

        list(j = j, mu = as.numeric(mu), lwr = as.numeric(lwr), upr = as.numeric(upr), ce = ce, fit = if (exists("fit_gam")) fit_gam else NULL)
      },
      .options  = furrr::furrr_options(
        seed     = TRUE,
        packages = c("mgcv","stats"),
        globals  = list(fit_fun = fit_fun, pred_fun = pred_fun,
                        rescale_fn = rescale_fn, sv_fn = sv_fn, clamp_fn = clamp_fn,
                        pred_ocat = if (is_ocat) pred_ocat else NULL,
                        is_ocat = is_ocat)
      ),
      .progress = progress
    )

    # write back & store empirical complexity and fitted GAMs
    ce_vec <- setNames(numeric(length(ov_all)), ov_all)
    gf     <- vector("list", length(ov_all)); names(gf) <- ov_all
    for (res in item_results) {
      if (is.null(res)) next
      nm_est <- paste0("e_est_yhat_", res$j)
      nm_lwr <- paste0("e_lwr_yhat_", res$j)
      nm_upr <- paste0("e_upr_yhat_", res$j)
      out[rows_g, nm_est] <- res$mu
      out[rows_g, nm_lwr] <- res$lwr
      out[rows_g, nm_upr] <- res$upr
      ce_vec[[res$j]] <- res$ce
      gf[[res$j]]     <- res$fit
    }
    ce_map[[as.character(g)]] <- ce_vec
    gam_fits[[as.character(g)]] <- gf
    vmsg(sprintf("  Group %d: done.", g))
  }
  vmsg("Empirical GAM curves done.")

  # -- Numeric types cleanup ---------------------------------------------------
  is_num <- vapply(out, is.numeric, logical(1L))
  keep_int <- names(out) %in% c(".rid", ".gid")
  to_double <- is_num & !keep_int
  for (nm in names(out)[to_double]) out[[nm]] <- as.double(out[[nm]])

  # -- 3) Metrics per item & group --------------------------------------------
  vmsg("Computing agreement metrics (r2, rmse, mae, penalized variants)...")
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
  else tibble::tibble(
    .gid   = integer(),
    .group = character(),
    item   = character(),
    type   = character(),
    n_eff  = integer(),
    c_m    = double(),
    c_e    = double(),
    k_eff  = double(),
    r2     = double(),
    rmse   = double(),
    mae    = double(),
    r2_pen   = double(),
    rmse_pen = double(),
    mae_pen  = double()
  )
  vmsg("...metrics done.")

  # -- Prepare smoothed latent-grid curves as well ------------------------------
  # Speed-ups:
  # 1) Do NOT check for prepare() existence — it is part of the package.
  # 2) Reuse already computed original_data (with factor scores etc.) as input to
  #    avoid recomputing lavPredict_parallel inside prepare().
  # 3) Forward parallelization knobs (plan, workers) to prepare() so that it can
  #    leverage the already chosen backend.
  {
    vmsg("Building new_data via prepare() using original_data as input...")
    new_data <- tryCatch(
      prepare(
        fit   = fit,
        info  = info,
        # Key change: feed 'out' (original_data with factor scores and predictions)
        # so prepare() can skip recomputing scores and simply build latent grids.
        data  = out,
        level = level,
        # Pass through parallelization knobs to keep a single planning locus
        plan    = plan,
        workers = workers,
        cluster = cluster
      ),
      error = function(e) {
        warning(sprintf("prepare() failed inside item_fit_data(): %s", conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(new_data)) vmsg("...prepare() done.") else vmsg("...prepare() failed (see warning).")
  }

  # -- Empirical predictions on new_data (latent-by-latent grid) --------------
  if (!is.null(new_data)) {
    vmsg("Scoring empirical predictions for new_data latent grids...")
    # Ensure needed columns
    if (!all(c(".gid", ".latent_var") %in% names(new_data))) {
      warning("new_data lacks .gid or .latent_var; skipping empirical predictions for new_data.")
    } else {
      for (g in seq_len(n_groups)) {
        rows_ng <- which(new_data$.gid == g)
        if (!length(rows_ng)) next
        nd_g <- new_data[rows_ng, , drop = FALSE]

        # pull fitted models for this group
        gf <- gam_fits[[as.character(g)]]
        if (is.null(gf)) next

        # factor being varied in each row
        lat_in_row <- as.character(nd_g$.latent_var)

        for (j in ov_all) {
          fit_gj <- gf[[j]]
          if (is.null(fit_gj)) next

          # determine relevant latents for this item in this group
          rel_lats <- if (vary_by_group) {
            li_g <- loads_by_item_by_group[[as.character(g)]] %||% loads_by_item_by_group[[g]]
            (if (!is.null(li_g)) li_g[[j]] else NULL) %||% loads_by_item_overall[[j]]
          } else {
            loads_by_item_overall[[j]]
          }
          rel_lats <- rel_lats %||% character(0L)
          if (!length(rel_lats)) next

          # We only predict when the currently varied latent is in rel_lats
          for (k in unique(lat_in_row)) {
            if (!k %in% rel_lats) next
            sel_k <- which(lat_in_row == k)
            if (!length(sel_k)) next

            nd_pred <- nd_g[sel_k, rel_lats, drop = FALSE]

            if (is_ocat && (j %in% ov_ord)) {
              oc <- pred_ocat(fit_gj, newdata = nd_pred, level = level)
              pr_mat <- try(suppressWarnings(predict(fit_gj, newdata = nd_pred, type = "response")), silent = TRUE)
              vals <- ord_value_map[[j]]
              if (is.matrix(pr_mat) && !is.null(vals) && ncol(pr_mat) == length(vals)) {
                mu  <- as.numeric(pr_mat %*% vals)
                lwr <- .emap_expect(oc$E_lwr, vals)
                upr <- .emap_expect(oc$E_upr, vals)
                pred_df <- data.frame(
                  setNames(list(mu),  paste0("e_est_", j, "_", k)),
                  setNames(list(lwr), paste0("e_lwr_", j, "_", k)),
                  setNames(list(upr), paste0("e_upr_", j, "_", k))
                )
              } else {
                if (is.null(vals)) {
                  R <- length(fit_gj$family$getTheta(trans = TRUE)) + 1L
                  rj <- rng_by_item[[j]]
                  to_original_from_scores <- function(E) {
                    if (!is.finite(rj[1]) || !is.finite(rj[2]) || R <= 1L) return(E)
                    rj[1] + (E - 1) * ((rj[2] - rj[1]) / (R - 1))
                  }
                  pred_df <- data.frame(
                    setNames(list(to_original_from_scores(oc$E_est)), paste0("e_est_", j, "_", k)),
                    setNames(list(to_original_from_scores(oc$E_lwr)), paste0("e_lwr_", j, "_", k)),
                    setNames(list(to_original_from_scores(oc$E_upr)), paste0("e_upr_", j, "_", k))
                  )
                } else {
                  pred_df <- data.frame(
                    setNames(list(.emap_expect(oc$E_est, vals)), paste0("e_est_", j, "_", k)),
                    setNames(list(.emap_expect(oc$E_lwr, vals)), paste0("e_lwr_", j, "_", k)),
                    setNames(list(.emap_expect(oc$E_upr, vals)), paste0("e_upr_", j, "_", k))
                  )
                }
              }
            } else {
              # default branch (continuous or beta)
              pred_df <- pred_fun(
                dat      = nd_pred,
                fit      = fit_gj,
                prefixes = c("e_est_", "e_lwr_", "e_upr_"),
                suffix   = paste0(j, "_", k),
                level    = level
              )

              if (j %in% ov_ord && !is_ocat) {
                # map [0,1] -> original scale
                rj <- rng_by_item[[j]]
                to_original_scale <- function(p01) { d <- rj[2] - rj[1]; p01 * d + rj[1] }
                cn <- names(pred_df)
                pred_df[[cn[1]]] <- to_original_scale(clamp_fn(pred_df[[cn[1]]]))
                pred_df[[cn[2]]] <- to_original_scale(clamp_fn(pred_df[[cn[2]]]))
                pred_df[[cn[3]]] <- to_original_scale(clamp_fn(pred_df[[cn[3]]]))
              }
            }

            # write back only into selected rows
            for (cn in names(pred_df)) {
              new_data[rows_ng[sel_k], cn] <- pred_df[[cn]]
            }
          }
        }
      }
    }
    vmsg("...empirical predictions for new_data done.")
  }

  # -- Return both data and metrics -------------------------------------------
  vmsg("Assembling final output list...")
  return(list(
    original_data = tibble::as_tibble(out),
    metrics       = metrics_tbl,
    new_data      = new_data
  ))
}

# ---- Helper: Ordered-category expected value + CI (E-only, delta method) -----
# Returns only E_est, E_lwr, E_upr for efficiency & stability.
.predict_ocat_ci_delta <- function(fit, newdata, level = 0.95) {
  stopifnot(inherits(fit, "gam"))

  pr <- predict(fit, newdata = newdata, type = "link", se.fit = TRUE)
  eta_hat <- as.numeric(pr$fit)
  se_hat  <- pmax(0, as.numeric(pr$se.fit))
  n <- length(eta_hat)

  theta <- fit$family$getTheta(trans = TRUE)
  if (is.null(theta)) theta <- numeric(0)
  Rm1 <- length(theta)
  if (Rm1 < 1L) stop("ocat family returned empty theta; ensure R was specified when fitting.")
  R <- Rm1 + 1L

  # helpers that always return matrices
  Fcum <- function(et) {
    mat <- plogis(outer(et, theta, function(e, th) th - e))
    dim(mat) <- c(length(et), Rm1)
    mat
  }
  flog <- function(et) {
    F <- Fcum(et)
    F * (1 - F)
  }
  make_probs <- function(Fc) {
    Fc <- as.matrix(Fc)
    out <- matrix(NA_real_, nrow = nrow(Fc), ncol = ncol(Fc) + 1L)
    out[, 1] <- Fc[, 1]
    if (ncol(Fc) > 1L) out[, 2:ncol(Fc)] <- Fc[, 2:ncol(Fc), drop = FALSE] - Fc[, 1:(ncol(Fc)-1L), drop = FALSE]
    out[, ncol(out)] <- 1 - Fc[, ncol(Fc)]
    rs <- rowSums(out)
    rs[!is.finite(rs) | rs <= 0] <- 1
    out / rs
  }

  # point expectation E[Y]
  P_est <- make_probs(Fcum(eta_hat))
  r_idx <- matrix(rep(seq_len(R), each = n), nrow = n)
  E_est <- rowSums(P_est * r_idx)

  # delta-method CI for E[Y]
  z <- stats::qnorm(1 - (1 - level)/2)
  f_vals <- flog(eta_hat)
  dP1 <- -f_vals[, 1, drop = TRUE]
  dPR <-  f_vals[, ncol(f_vals), drop = TRUE]
  dPmid <- if (R > 2L) f_vals[, 1:(R-2), drop = FALSE] - f_vals[, 2:(R-1), drop = FALSE] else NULL
  dP <- matrix(0, nrow = n, ncol = R)
  dP[, 1] <- dP1
  if (!is.null(dPmid)) dP[, 2:(R-1)] <- dPmid
  dP[, R] <- dPR
  dE <- rowSums(dP * r_idx)
  se_E <- abs(dE) * se_hat
  E_lwr <- pmax(1, pmin(R, E_est - z * se_E))
  E_upr <- pmax(1, pmin(R, E_est + z * se_E))

  out <- data.frame(E_est = E_est, E_lwr = E_lwr, E_upr = E_upr)
  out
}
