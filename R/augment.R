#' Unified augment(): one lavPredict pass, then type-specific add-ons
#'
#' Internals:
#'   - .augment_ordinal_dat(fit, dat, ...)    # expects FS+data from lavPredict_parallel()
#'   - .augment_continuous_dat(fit, dat, ...) # expects FS+data from lavPredict_parallel()
#'
#' Notes:
#'   * Factor-score SEs are only available for continuous-only models.
#'     We therefore request SEs from lavPredict_parallel() only when no ordinal ov is present.
#'   * Delta-method SEs / CIs for .yhat_* are available for both types (if vcov exists).
#'
#' @export
augment <- function(fit,
                    # shared toggles
                    yhat            = TRUE,
                    resid           = TRUE,
                    ci              = TRUE,
                    level           = 0.95,
                    # ordinal-only toggles
                    ystar           = TRUE,
                    pr              = TRUE,
                    # SE controls
                    se_fs           = TRUE,   # FS SEs (continuous-only)
                    se_yhat         = TRUE,   # delta-method SE(.yhat_*)
                    # prefixes
                    prefix_yhat     = ".yhat_",
                    prefix_ci       = c(".yhat_lwr_", ".yhat_upr_"),
                    prefix_resid    = ".resid_",
                    prefix_ystar    = ".ystar_",
                    prefix_pr       = ".pr_",
                    prefix_se       = ".se_",
                    prefix_se_yhat  = ".se_yhat_",
                    sep             = "__",
                    vcov_type       = NULL,
                    # lavPredict_parallel controls (optional)
                    workers         = NULL,
                    plan            = c("auto","multisession","multicore","sequential"),
                    chunk_size      = NULL,
                    progress        = FALSE,
                    ...) {

  # -- Basic checks ------------------------------------------------------------
  .assert_lavaan_fit(fit, require_converged = TRUE, require_latent = TRUE)

  plan <- match.arg(plan)
  info <- model_info(fit)

  ov_ord  <- info$ov_ordinal
  ov_cont <- info$ov_continuous
  has_ord  <- length(ov_ord)  > 0L
  has_cont <- length(ov_cont) > 0L

  # -- Request FS SEs only if continuous-only ---------------------------------
  request_fs_se <- isTRUE(se_fs) && !has_ord

  # One and only lavPredict run (returns observed vars + factor scores [+ FS SEs if possible])
  dat <- lavPredict_parallel(
    fit,
    workers     = workers,
    plan        = plan,
    chunk_size  = chunk_size,
    return_type = "data",
    progress    = progress,
    se          = request_fs_se,
    prefix_se   = prefix_se,
    ...
  )

  # Ensure .gid and .group exist and are ordered nicely
  if (!(".gid" %in% names(dat))) {
    # Single-group fallback
    dat$.gid <- 1L
  }
  # If .group missing, try to add labels from model_info
  if (!(".group" %in% names(dat))) {
    glabs <- info$group_labels %||% as.character(sort(unique(dat$.gid)))
    dat$.group <- glabs[dat$.gid]
  }
  # Relocate for stable order: .gid, .group, then factor scores
  anchor <- ".gid"
  dat <- dplyr::relocate(dat, .group, .after = dplyr::all_of(anchor))

  # Gate: nothing to do?
  if (!has_ord && !has_cont) return(tibble::as_tibble(dat))

  # --- Ordinal branch on top of the single lavPredict dataset -----------------
  if (has_ord) {
    dat <- .augment_ordinal_dat(
      fit            = fit,
      dat            = dat,
      ov_ord         = ov_ord,
      ystar          = ystar,
      yhat           = yhat,
      pr             = pr,
      ci             = ci,
      level          = level,
      resid          = resid,
      se_yhat        = se_yhat,
      prefix_ystar   = prefix_ystar,
      prefix_yhat    = prefix_yhat,
      prefix_pr      = prefix_pr,
      prefix_ci      = prefix_ci,
      prefix_resid   = prefix_resid,
      prefix_se_yhat = prefix_se_yhat,
      sep            = sep,
      vcov_type      = vcov_type
    )
  }

  # --- Continuous branch on top of the same dataset ---------------------------
  if (has_cont) {
    dat <- .augment_continuous_dat(
      fit            = fit,
      dat            = dat,
      ov_cont        = ov_cont,
      yhat           = yhat,
      resid          = resid,
      ci             = ci,
      level          = level,
      se_yhat        = se_yhat,
      prefix_yhat    = prefix_yhat,
      prefix_ci      = prefix_ci,
      prefix_resid   = prefix_resid,
      prefix_se_yhat = prefix_se_yhat,
      vcov_type      = vcov_type
    )
  }

  # Cast numerics robustly (preserve factors for ov_ord)
  dat <- dplyr::mutate(dat, dplyr::across(dplyr::where(is.numeric), as.double))
  tibble::as_tibble(dat)
}
