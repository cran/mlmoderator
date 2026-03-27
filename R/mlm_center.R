#' Center variables for multilevel modeling
#'
#' Performs grand-mean centering, group-mean centering, or both (within-between
#' decomposition) on one or more variables in a data frame. Group-mean centering
#' is the standard preparation for cross-level interaction models.
#'
#' @param data A data frame.
#' @param vars Character vector of variable names to center.
#' @param cluster Character scalar: name of the clustering variable (required
#'   when `type` is `"group"` or `"both"`).
#' @param type One of `"grand"`, `"group"`, or `"both"`.
#'   * `"grand"`: subtract the overall mean.
#'   * `"group"`: subtract the cluster mean (within-person / within-school
#'     centering).
#'   * `"both"`: return both the within-cluster-centered value *and* the
#'     cluster mean (between component), appended as new columns.
#' @param suffix_within Suffix appended to within-centered variable names when
#'   `type = "both"`. Default is `"_within"`.
#' @param suffix_between Suffix appended to between (cluster mean) variable
#'   names when `type = "both"`. Default is `"_between"`.
#'
#' @return The input data frame with new centered columns appended. Original
#'   columns are not modified.
#'
#' @examples
#' data(school_data)
#'
#' # Grand-mean center SES
#' d1 <- mlm_center(school_data, vars = "ses", type = "grand")
#' head(d1[, c("ses", "ses_c")])
#'
#' # Group-mean center SES within schools
#' d2 <- mlm_center(school_data, vars = "ses", cluster = "school", type = "group")
#' head(d2[, c("ses", "ses_c")])
#'
#' # Within-between decomposition
#' d3 <- mlm_center(school_data, vars = "ses", cluster = "school", type = "both")
#' head(d3[, c("ses", "ses_within", "ses_between")])
#'
#' @export
mlm_center <- function(data,
                       vars,
                       cluster = NULL,
                       type = c("grand", "group", "both"),
                       suffix_within  = "_within",
                       suffix_between = "_between") {
  type <- match.arg(type)

  if (!is.data.frame(data)) {
    rlang::abort("`data` must be a data frame.")
  }
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0) {
    rlang::abort(paste0("Variable(s) not found in data: ",
                        paste(missing_vars, collapse = ", ")))
  }
  if (type %in% c("group", "both") && is.null(cluster)) {
    rlang::abort("`cluster` must be specified for group-mean or within-between centering.")
  }
  if (!is.null(cluster) && !cluster %in% names(data)) {
    rlang::abort(paste0("Cluster variable '", cluster, "' not found in data."))
  }

  for (v in vars) {
    if (!is.numeric(data[[v]])) {
      rlang::abort(paste0("Variable '", v, "' must be numeric."))
    }

    if (type == "grand") {
      gm <- mean(data[[v]], na.rm = TRUE)
      data[[paste0(v, "_c")]] <- data[[v]] - gm

    } else if (type == "group") {
      grp_means <- tapply(data[[v]], data[[cluster]], mean, na.rm = TRUE)
      matched   <- as.numeric(grp_means[as.character(data[[cluster]])])
      data[[paste0(v, "_c")]] <- data[[v]] - matched

    } else if (type == "both") {
      grp_means <- tapply(data[[v]], data[[cluster]], mean, na.rm = TRUE)
      between   <- as.numeric(grp_means[as.character(data[[cluster]])])
      within    <- data[[v]] - between
      data[[paste0(v, suffix_within)]]  <- as.numeric(within)
      data[[paste0(v, suffix_between)]] <- as.numeric(between)
    }
  }

  data
}
