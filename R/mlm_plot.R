# Suppress R CMD check NOTEs for ggplot2 data-masking variables
utils::globalVariables(c(
  ".data", "ci_lo", "ci_hi", "fitted", ".modx_label",
  "modx_value", "slope", "ci_lower", "ci_upper", "region",
  "pred_val", "modx_val", "yhat", "y",
  "icc", "se_int_adj", "sig", "jn_bound_adj",
  "pi_lower", "pi_upper", "pct_random",
  "b_int", "influential", "label", "cluster"
))

#' Publication-ready interaction plot for multilevel models
#'
#' Creates a `ggplot2`-based interaction plot showing predicted values of the
#' outcome across levels of `pred`, with separate lines for each selected value
#' of `modx`. Confidence bands and raw data overlay are optional.
#'
#' @param model An `lmerMod` object with a two-way interaction between `pred`
#'   and `modx` in the fixed-effects.
#' @param pred Character scalar. Focal predictor (x-axis).
#' @param modx Character scalar. Moderator (separate lines).
#' @param modx.values Strategy for moderator values. Same options as
#'   `mlm_probe()`: `"mean-sd"`, `"quartiles"`, `"tertiles"`, `"custom"`.
#' @param at Numeric vector of custom moderator values (used when
#'   `modx.values = "custom"`).
#' @param interval Logical. Draw confidence bands? Default `TRUE`.
#' @param conf.level Confidence level for bands. Default `0.95`.
#' @param points Logical. Overlay raw data points? Default `FALSE`.
#' @param point_alpha Transparency for raw data points. Default `0.3`.
#' @param colors Character vector of colours for moderator lines. If `NULL`,
#'   uses a accessible default palette.
#' @param line_size Line width for predicted lines. Default `1`.
#' @param x_label Label for x-axis. Defaults to `pred`.
#' @param y_label Label for y-axis. Defaults to response variable name.
#' @param legend_title Label for the legend. Defaults to `modx`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(
#'   y   = rnorm(200),
#'   x   = rnorm(200),
#'   m   = rep(rnorm(20), each = 10),
#'   grp = factor(rep(1:20, each = 10))
#' )
#' dat$y <- dat$y + dat$x * dat$m
#' mod <- lme4::lmer(y ~ x * m + (1 | grp), data = dat,
#'                   control = lme4::lmerControl(optimizer = "bobyqa"))
#' mlm_plot(mod, pred = "x", modx = "m")
#' mlm_plot(mod, pred = "x", modx = "m", modx.values = "quartiles")
#'
#' @export
mlm_plot <- function(model,
                     pred,
                     modx,
                     modx.values  = c("mean-sd", "quartiles", "tertiles", "custom"),
                     at           = NULL,
                     interval     = TRUE,
                     conf.level   = 0.95,
                     points       = FALSE,
                     point_alpha  = 0.3,
                     colors       = NULL,
                     line_size    = 1,
                     x_label      = NULL,
                     y_label      = NULL,
                     legend_title = NULL) {

  .check_lmer(model)
  .validate_terms(model, pred, modx)

  modx.values <- match.arg(modx.values)
  mf          <- model@frame
  outcome     <- names(mf)[1]
  modx_x      <- mf[[modx]]

  modx_vals <- .pick_modx_values(modx_x, modx.values = modx.values, at = at)
  grid      <- .make_prediction_grid(model, pred, modx, modx_vals)

  # Predict using fixed effects only (re.form = NA)
  grid$fitted <- stats::predict(model, newdata = grid, re.form = NA)

  # Confidence bands via delta method over the prediction
  if (interval) {
    ci_list <- lapply(seq_len(nrow(grid)), function(i) {
      w <- grid[[modx]][i]
      ss <- .simple_slope_linear(model, pred, modx, w, conf.level = conf.level)
      # SE of predicted value at this pred value
      # y_hat = b0 + b_pred*x + b_modx*w + b_int*x*w + ...
      # For plotting, we compute SE of the simple slope * pred value + intercept uncertainty
      # Use a simpler approximation: SE of fit from vcov
      se_fit <- .se_of_fit(model, grid[i, , drop = FALSE], pred, modx)
      data.frame(se_fit = se_fit)
    })
    ci_df <- do.call(rbind, ci_list)
    df_resid <- get_residual_df(model)
    t_crit <- stats::qt(1 - (1 - conf.level) / 2, df = df_resid)
    grid$ci_lo <- grid$fitted - t_crit * ci_df$se_fit
    grid$ci_hi <- grid$fitted + t_crit * ci_df$se_fit
  }

  # Build factor labels for moderator
  modx_labels <- .build_modx_labels(modx_vals, modx, modx.values)
  grid$.modx_label <- factor(
    modx_labels[as.character(grid$.modx_val)],
    levels = modx_labels
  )

  # Colors
  if (is.null(colors)) {
    colors <- c("#2166AC", "#D6604D", "#1A9641", "#7B2D8B", "#FF7F00")
    colors <- colors[seq_along(modx_vals)]
  }

  # x / y / legend labels
  if (is.null(x_label))      x_label      <- pred
  if (is.null(y_label))      y_label      <- outcome
  if (is.null(legend_title)) legend_title <- modx

  p <- ggplot2::ggplot(
    grid,
    ggplot2::aes(
      x     = .data[[pred]],
      y     = .data[["fitted"]],
      color = .data[[".modx_label"]],
      fill  = .data[[".modx_label"]]
    )
  )

  if (points) {
    p <- p + ggplot2::geom_point(
      data  = mf,
      ggplot2::aes(x = .data[[pred]], y = .data[[outcome]]),
      color = "grey60",
      alpha = point_alpha,
      size  = 1,
      inherit.aes = FALSE
    )
  }

  if (interval) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_lo, ymax = ci_hi),
      alpha = 0.15, color = NA
    )
  }

  p <- p +
    ggplot2::geom_line(linewidth = line_size) +
    ggplot2::scale_color_manual(values = colors, name = legend_title) +
    ggplot2::scale_fill_manual(values  = colors, name = legend_title) +
    ggplot2::labs(x = x_label, y = y_label) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      legend.position   = "right",
      panel.grid.major  = ggplot2::element_line(color = "grey92"),
      strip.background  = ggplot2::element_rect(fill = "grey95", color = NA)
    )

  p
}

# ---- Internal helpers --------------------------------------------------------

#' SE of fitted value using the delta method over fixed-effects
#' @noRd
.se_of_fit <- function(model, newrow, pred, modx) {
  fe    <- lme4::fixef(model)
  vcv   <- .extract_vcov(model)
  fe_nm <- names(fe)

  x_val <- newrow[[pred]]
  w_val <- newrow[[modx]]
  int_term <- .get_interaction_term(model, pred, modx)

  # Build gradient vector (partial derivatives wrt each beta)
  grad <- stats::setNames(rep(0, length(fe)), fe_nm)

  for (nm in fe_nm) {
    if (nm == "(Intercept)") {
      grad[nm] <- 1
    } else if (nm == pred) {
      grad[nm] <- x_val
    } else if (nm == modx) {
      grad[nm] <- w_val
    } else if (nm == int_term) {
      grad[nm] <- x_val * w_val
    }
    # Other covariates held at mean => value is mean => gradient * mean
    # (already absorbed into fitted value; SE contribution modest)
  }

  sqrt(as.numeric(t(grad) %*% vcv %*% grad))
}
