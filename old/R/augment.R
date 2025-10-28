#' Augment SEM data with latent predictions, residuals, and (optional) CIs
#'
#' Public user-facing wrapper that augments fitted `lavaan` models with predicted
#' values, residuals, optional delta-method SEs/CIs, and (for ordinal models) also
#' latent linear predictors (y*) and per-category probabilities. Supports
#' continuous-only, ordinal-only, and mixed (both types) models.
#'
#' Internally routes to `.augment_continuous()` and/or `.augment_ordinal()` and
#' reuses `data` (lavPredict output) and `info` (model_info) so they are computed
#' at most once.
#'
#' @param fit A fitted `lavaan` object.
#' @param data Optional lavPredict output to reuse. Either a data.frame (single-group)
#'   or a list of data.frames (per group) as returned by `lavPredict_parallel()`.
#'   If `NULL`, it will be computed once and reused.
#' @param info Optional `model_info(fit)` list to reuse. If `NULL`, it will be computed.
#' @param # ---- Common toggles ----
#' @param yhat Logical; include predicted observed values (default TRUE).
#' @param resid Logical; include residuals (obs - yhat) (default TRUE).
#' @param ci Logical; include delta-method CIs for `yhat` (default TRUE).
#' @param level CI level (default 0.95).
#' @param se_yhat Logical; include delta-method SEs of `yhat` (default TRUE).
#' @param # ---- Ordinal-only toggles ----
#' @param ystar Logical; include latent linear predictors y* for ordinal items (default TRUE).
#' @param pr Logical; include per-category probabilities for ordinal items (default TRUE).
#' @param # ---- Factor-score SEs (continuous-only) ----
#' @param se_fs Logical; request factor-score SEs in `lavPredict_parallel()` for continuous-only
#'   models (default TRUE). Ignored when any ordinal indicator is present.
#' @param vcov_type Optional `vcov()` type passed to the continuous branch for robust CIs.
#' @param # ---- Column prefixes ----
#' @param prefix_ystar,prefix_yhat,prefix_pr,prefix_ci,prefix_resid,
#'   prefix_se_fs,prefix_se_yhat Character prefixes for generated columns. Must be consistent
#'   with the internal augmenters' expectations (defaults mirror internal functions).
#'
#' @return A tibble-like `data.frame` with `.rid`, `.gid`, `.group`, original lavPredict
#'   columns (observed + factor scores [+ FS SEs if requested/available]) and appended
#'   augmentation columns from the relevant branches.
#'
#' @export
augment <- function(fit,
                    data            = NULL,
                    info            = NULL,
                    # common toggles
                    yhat            = TRUE,
                    resid           = TRUE,
                    ci              = TRUE,
                    level           = 0.95,
                    se_yhat         = TRUE,
                    # ordinal-only toggles
                    ystar           = TRUE,
                    pr              = TRUE,
                    # FS SEs (cont-only)
                    se_fs           = TRUE,
                    vcov_type       = NULL,
                    # prefixes (kept aligned with internals)
                    prefix_ystar    = ".ystar_",
                    prefix_yhat     = ".yhat_",
                    prefix_pr       = ".pr_",
                    prefix_ci       = c(".yhat_lwr_", ".yhat_upr_"),
                    prefix_resid    = ".resid_",
                    prefix_se_fs    = ".se_",
                    prefix_se_yhat  = ".se_yhat_",
                    sep             = "__") {

  # -- Assertions on fit (light; internals do strict checks) -------------------
  .assert_lavaan_fit(fit)

  # -- One-time model info -----------------------------------------------------
  if (is.null(info)) info <- model_info(fit)
  ov_cont <- info$ov_continuous
  ov_ord  <- info$ov_ordinal

  # -- One-time lavPredict data (reused) ---------------------------------------
  # We only compute here if user didn't provide `data`.
  if (is.null(data)) {
    # Ask for FS SEs only if purely continuous; internal cont-augmenter also guards this.
    request_fs_se <- isTRUE(se_fs) && length(ov_ord) == 0L
    data <- lavPredict_parallel(
      fit,
      return_type  = "list",
      se           = request_fs_se,
      prefix_se_fs = prefix_se_fs
    )
  }

  # -- Routing by measurement type --------------------------------------------
  has_cont <- length(ov_cont) > 0L
  has_ord  <- length(ov_ord)  > 0L

  if (!has_cont && !has_ord) {
    stop("The fitted model has neither continuous nor ordinal observed indicators.", call. = FALSE)
  }

  # Helper to extract "base" columns (rid/gid/group + original lavPredict table)
  # from an augmented table and to isolate augmentation-only columns.
  split_base_aug <- function(df) {
    # base columns are .rid, .gid, .group and everything that was already in lavPredict
    # Heuristic: any column that does NOT start with known augmentation prefixes is base.
    aug_prefixes <- c(prefix_ystar, prefix_yhat, prefix_pr, prefix_resid, prefix_se_yhat,
                      prefix_ci[1L], prefix_ci[2L])
    starts_with_any <- function(nm, pref) any(startsWith(nm, pref))
    is_aug <- vapply(names(df), starts_with_any, logical(1L), pref = aug_prefixes)
    list(base = df[, !is_aug, drop = FALSE],
         aug  = df[,  is_aug, drop = FALSE])
  }

  out <- NULL

  if (has_cont && !has_ord) {
    # Continuous-only branch
    out <- .augment_continuous(
      fit            = fit,
      data           = data,
      info           = info,
      yhat           = yhat,
      ci             = ci,
      level          = level,
      resid          = resid,
      se_fs          = se_fs,
      se_yhat        = se_yhat,
      prefix_yhat    = prefix_yhat,
      prefix_ci      = prefix_ci,
      prefix_resid   = prefix_resid,
      prefix_se_fs   = prefix_se_fs,
      prefix_se_yhat = prefix_se_yhat,
      vcov_type      = vcov_type
    )
    return(out)
  }

  if (!has_cont && has_ord) {
    # Ordinal-only branch
    out <- .augment_ordinal(
      fit            = fit,
      data           = data,
      info           = info,
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
      sep            = sep
    )
    return(out)
  }

  # Mixed (both continuous and ordinal): run both and merge augmentation columns
  aug_cont <- .augment_continuous(
    fit            = fit,
    data           = data,
    info           = info,
    yhat           = yhat,
    ci             = ci,
    level          = level,
    resid          = resid,
    se_fs          = se_fs,        # internally ignored if ord present
    se_yhat        = se_yhat,
    prefix_yhat    = prefix_yhat,
    prefix_ci      = prefix_ci,
    prefix_resid   = prefix_resid,
    prefix_se_fs   = prefix_se_fs,
    prefix_se_yhat = prefix_se_yhat,
    vcov_type      = vcov_type
  )

  aug_ord <- .augment_ordinal(
    fit            = fit,
    data           = data,
    info           = info,
    ystar          = ystar,
    yhat          = yhat,
    pr            = pr,
    ci            = ci,
    level         = level,
    resid         = resid,
    se_yhat       = se_yhat,
    prefix_ystar  = prefix_ystar,
    prefix_yhat   = prefix_yhat,
    prefix_pr     = prefix_pr,
    prefix_ci     = prefix_ci,
    prefix_resid  = prefix_resid,
    prefix_se_yhat = prefix_se_yhat,
    sep           = sep
  )

  # Split base vs augmentation parts
  sbc <- split_base_aug(aug_cont)
  sbo <- split_base_aug(aug_ord)

  # Use base from the continuous branch (arbitrary but consistent)
  base <- sbc$base

  # Sanity: ensure row alignment by .rid/.gid if both exist
  # (Assumes both augmenters preserve order; cheap alignment guard.)
  if (all(c(".rid", ".gid") %in% names(base)) &&
      all(c(".rid", ".gid") %in% names(sbo$base))) {
    # if necessary, align sbo to base by (.rid,.gid)
    # Build index map; if identical, match returns 1:nrow
    key_base <- paste0(base$.rid, "_", base$.gid)
    key_ord  <- paste0(sbo$base$.rid, "_", sbo$base$.gid)
    idx <- match(key_base, key_ord)
    if (!all(!is.na(idx) & idx == seq_along(idx))) {
      sbo$aug <- sbo$aug[idx, , drop = FALSE]
    }
  }

  # Drop any augmentation-name overlaps (should be none; guard anyway)
  dup_cols <- intersect(names(sbc$aug), names(sbo$aug))
  if (length(dup_cols)) {
    # keep continuous' versions, drop duplicates from ordinal aug
    sbo$aug <- sbo$aug[, setdiff(names(sbo$aug), dup_cols), drop = FALSE]
  }

  # Final bind: base + cont_aug + ord_aug
  out <- cbind(base, sbc$aug, sbo$aug, stringsAsFactors = FALSE)

  # Return class consistent with internals
  class(out) <- c("tbl_df", "tbl", "data.frame")
  out
}
