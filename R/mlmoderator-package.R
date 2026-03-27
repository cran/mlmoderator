#' mlmoderator: Probing and Visualizing Multilevel Interaction Effects
#'
#' @description
#' `mlmoderator` provides a unified workflow for estimating, probing, and
#' visualizing multilevel moderation effects from mixed-effects models fitted
#' with `lme4::lmer()`. It supports simple slopes analysis, Johnson---Neyman
#' intervals, publication-ready interaction plots, and grand- or group-mean
#' centering.
#'
#' ## Core functions
#'
#' | Function | Description |
#' |---|---|
#' | `mlm_center()` | Grand- or group-mean center variables |
#' | `mlm_probe()` | Compute simple slopes at selected moderator values |
#' | `mlm_jn()` | Johnson-Neyman significance regions |
#' | `mlm_plot()` | Publication-ready interaction plot |
#' | `mlm_summary()` | Consolidated moderation report |
#' | `mlm_variance_decomp()` | Decompose slope uncertainty into fixed + random components |
#' | `mlm_surface()` | Slope surface heatmap over the full predictor x moderator space |
#' | `mlm_sensitivity()` | Robustness of interaction to ICC shift and omitted confounders |
#'
#' ## Typical workflow
#'
#' ```r
#' library(mlmoderator)
#' library(lme4)
#'
#' data(school_data)
#'
#' # 1. Center variables
#' dat <- mlm_center(school_data, vars = "ses", cluster = "school", type = "group")
#'
#' # 2. Fit model
#' mod <- lmer(math ~ ses * climate + gender + (1 + ses | school), data = dat)
#'
#' # 3. Probe interaction
#' mlm_probe(mod, pred = "ses", modx = "climate")
#'
#' # 4. Johnson-Neyman interval
#' mlm_jn(mod, pred = "ses", modx = "climate")
#'
#' # 5. Plot
#' mlm_plot(mod, pred = "ses", modx = "climate")
#'
#' # 6. Summary report
#' mlm_summary(mod, pred = "ses", modx = "climate")
#' ```
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom lme4 lmer fixef getME VarCorr lmerControl
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_point geom_hline
#'   geom_vline geom_contour geom_contour_filled geom_errorbar geom_text
#'   annotate scale_color_manual scale_fill_manual scale_x_discrete
#'   scale_fill_gradient2 scale_fill_viridis_c scale_fill_viridis_d
#'   coord_cartesian labs theme_classic theme
#'   element_line element_rect element_text unit
#' @importFrom grDevices colorRampPalette
#' @importFrom stats vcov sd quantile predict setNames qt pt formula
#' @importFrom rlang abort .data
#' @importFrom utils globalVariables
## usethis namespace: end
NULL
