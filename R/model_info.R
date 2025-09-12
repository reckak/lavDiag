#' Summarize basic model properties from a lavaan/blavaan fit
#'
#' @param fit A fitted lavaan or blavaan object.
#' @return A named list with model properties (groups, categorical, multilevel, etc.).
#' @examples
#' # info <- model_info(fit)
model_info <- function(fit) {
  # -- Defensive check ---------------------------------------------------------
  .assert_lavaan_fit(fit)

  # -- Groups ------------------------------------------------------------------
  n_groups <- tryCatch(lavaan::lavInspect(fit, "ngroups"), error = function(e) NA_integer_)
  is_single_group <- isTRUE(n_groups == 1L)

  # -- Categorical -------------------------------------------------------------
  is_categorical <- tryCatch(lavaan::lavInspect(fit, "categorical"),
                             error = function(e) NA)

  # -- Variables ---------------------------------------------------------------
  observed_variables <- tryCatch(lavaan::lavNames(fit, type = "ov"),
                                 error = function(e) character())
  latent_variables   <- tryCatch(lavaan::lavNames(fit, type = "lv"),
                                 error = function(e) character())

  # -- Group meta --------------------------------------------------------------
  group_labels <- tryCatch(lavaan::lavInspect(fit, "group.label"),
                           error = function(e) NULL)
  n_obs        <- tryCatch(lavaan::lavInspect(fit, "nobs"),
                           error = function(e) NULL)

  # -- Multilevel detection ----------------------------------------------------
  # Use public inspectors documented in lavaan manual:
  # - "nlevels": number of levels (>=2 implies multilevel)
  # - "cluster": name(s) of clustering variable(s) if specified
  n_levels <- tryCatch(lavaan::lavInspect(fit, "nlevels"),
                       error = function(e) NA_integer_)
  cluster_var <- tryCatch(lavaan::lavInspect(fit, "cluster"),
                          error = function(e) NULL)
  n_clusters <- tryCatch(lavaan::lavInspect(fit, "nclusters"),
                         error = function(e) NA_integer_)
  avg_cluster_size <- tryCatch(lavaan::lavInspect(fit, "average.cluster.size"),
                               error = function(e) NA_real_)

  is_multilevel <- isTRUE(!is.na(n_levels) && n_levels >= 2L) ||
    (!is.null(cluster_var) && length(cluster_var) > 0L)

  list(
    is_single_group    = is_single_group,
    is_categorical     = is_categorical,
    observed_variables = observed_variables,
    latent_variables   = latent_variables,
    group_labels       = group_labels,
    n_obs              = n_obs,
    n_groups           = n_groups,
    # --- multilevel summary ---
    is_multilevel      = is_multilevel,
    n_levels           = n_levels,
    cluster_var        = cluster_var,
    n_clusters         = n_clusters,
    average_cluster_size = avg_cluster_size
  )
}
