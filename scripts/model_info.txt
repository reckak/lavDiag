#' Summarize basic model properties from a lavaan/blavaan fit
#'
#' @param fit A fitted lavaan or blavaan object.
#' @return A named list with model properties (groups, variables, parameterization, etc.).
#' @noRd
model_info <- function(fit) {
  # -- Defensive check ---------------------------------------------------------
  .assert_lavaan_fit(fit)  # requires lavaan-compatible object

  # -- Groups ------------------------------------------------------------------
  n_groups        <- tryCatch(lavaan::lavInspect(fit, "ngroups"),
                              error = function(e) NA_integer_)
  is_single_group <- isTRUE(n_groups == 1L)

  # -- Basic meta --------------------------------------------------------------
  converged        <- tryCatch(lavaan::lavInspect(fit, "converged"),
                               error = function(e) NA)
  has_meanstructure<- tryCatch(lavaan::lavInspect(fit, "meanstructure"),
                               error = function(e) NA)

  # Pull lavaan options once; extract estimator & parameterization from there
  opts <- tryCatch(lavaan::lavInspect(fit, "options"), error = function(e) NULL)

  # Estimator (e.g., "ML", "WLSMV", "Bayes"); NA if missing
  estimator <- if (!is.null(opts) && !is.null(opts$estimator)) opts$estimator else NA_character_

  # Parameterization ("delta" default if missing)
  parameterization <- if (!is.null(opts) && !is.null(opts$parameterization))
    opts$parameterization else NA_character_

  # -- Categorical flag (lavaan-level) ----------------------------------------
  is_categorical <- tryCatch(lavaan::lavInspect(fit, "categorical"),
                             error = function(e) NA)

  # -- Variables ---------------------------------------------------------------
  observed_variables <- tryCatch(lavaan::lavNames(fit, type = "ov"),
                                 error = function(e) character())
  latent_variables   <- tryCatch(lavaan::lavNames(fit, type = "lv"),
                                 error = function(e) character())

  # Separate ordinal vs. continuous observed variables
  ov_ordinal    <- tryCatch(lavaan::lavNames(fit, type = "ov.ord"),
                            error = function(e) character())
  ov_continuous <- setdiff(observed_variables, ov_ordinal)

  # -- Group meta --------------------------------------------------------------
  group_labels <- tryCatch(lavaan::lavInspect(fit, "group.label"),
                           error = function(e) NULL)
  n_obs        <- tryCatch(lavaan::lavInspect(fit, "nobs"),
                           error = function(e) NULL)
  group_var    <- tryCatch(lavaan::lavInspect(fit, "group"),
                           error = function(e) NULL)

  # -- Multilevel detection ----------------------------------------------------
  n_levels         <- tryCatch(lavaan::lavInspect(fit, "nlevels"),
                               error = function(e) NA_integer_)
  cluster_var      <- tryCatch(lavaan::lavInspect(fit, "cluster"),
                               error = function(e) NULL)
  n_clusters       <- tryCatch(lavaan::lavInspect(fit, "nclusters"),
                               error = function(e) NA_integer_)
  avg_cluster_size <- tryCatch(lavaan::lavInspect(fit, "average.cluster.size"),
                               error = function(e) NA_real_)

  is_multilevel <- isTRUE(!is.na(n_levels) && n_levels >= 2L) ||
    (!is.null(cluster_var) && length(cluster_var) > 0L)

  list(
    # --- high-level status ---
    converged          = converged,
    has_meanstructure  = has_meanstructure,
    estimator          = estimator,
    parameterization   = parameterization,

    # --- grouping ---
    is_single_group    = is_single_group,
    n_groups           = n_groups,
    group_var          = group_var,
    group_labels       = group_labels,
    n_obs              = n_obs,

    # --- variables ---
    observed_variables = observed_variables,
    latent_variables   = latent_variables,
    ov_ordinal         = ov_ordinal,
    ov_continuous      = ov_continuous,

    # --- categorical flag ---
    is_categorical     = is_categorical,

    # --- multilevel summary ---
    is_multilevel        = is_multilevel,
    n_levels             = n_levels,
    cluster_var          = cluster_var,
    n_clusters           = n_clusters,
    average_cluster_size = avg_cluster_size
  )
}
