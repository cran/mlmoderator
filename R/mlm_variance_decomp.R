#' Decompose uncertainty in the simple slope of a multilevel interaction
#'
#' In a random-slope model, the uncertainty around a simple slope has two
#' distinct sources that standard `mlm_probe()` collapses into one SE:
#'
#' 1. **Fixed-effect uncertainty** -- imprecision in the estimated average
#'    slope (\eqn{\beta_1 + \beta_3 \cdot w}), captured by the fixed-effect
#'    variance-covariance matrix.
#' 2. **Random-slope variance** -- genuine between-cluster heterogeneity in
#'    the slope of `pred` (\eqn{\tau_{11}}), which is *not* estimation error
#'    but real variation in effects across clusters.
#'
#' These answer different questions:
#' * Fixed-effect uncertainty: "How precisely do we know the *average* slope
#'   at this moderator value?"
#' * Random-slope variance: "How much does the slope *actually vary* across
#'   clusters, regardless of what the moderator does?"
#'
#' The function also reports a **prediction interval** for the slope in a
#' new (unobserved) cluster, which combines both sources and is the
#' appropriate uncertainty interval for making cluster-level predictions.
#'
#' @param model An `lmerMod` object with a random slope for `pred` and a
#'   two-way interaction between `pred` and `modx`.
#' @param pred Character scalar. Focal predictor name.
#' @param modx Character scalar. Moderator name.
#' @param modx.values Strategy for moderator values. See `mlm_probe()`.
#' @param at Optional numeric vector of custom moderator values.
#' @param conf.level Confidence level for fixed-effect CIs. Default `0.95`.
#'
#' @return An object of class `mlm_variance_decomp` (a list) with:
#'   * `decomp`: data frame with columns `modx_value`, `slope`,
#'     `se_fixed` (SE from fixed-effect vcov only),
#'     `tau11` (random-slope SD, if estimable),
#'     `se_total` (combined SE for prediction in new cluster),
#'     `ci_lower`, `ci_upper` (fixed-effect CI),
#'     `pi_lower`, `pi_upper` (prediction interval for new cluster),
#'     `pct_random` (% of total slope variance from random effects).
#'   * `tau11`: the random-slope variance (\eqn{\tau_{11}}).
#'   * `has_random_slope`: logical -- does the model include a random slope
#'     for `pred`?
#'   * Metadata: `pred`, `modx`, `conf.level`.
#'
#' @details
#' **Random-slope variance (\eqn{\tau_{11}}) interpretation:**
#' If `tau11` is large relative to the fixed slope, the effect of `pred`
#' varies substantially across clusters even *after* accounting for the
#' moderator. This is important: a significant interaction does not mean
#' the moderator fully explains between-cluster slope heterogeneity.
#'
#' **Prediction interval interpretation:**
#' The prediction interval answers: "For a randomly sampled new cluster at
#' this moderator value, what range of slopes should we expect?" It will
#' always be wider than the confidence interval because it incorporates
#' \eqn{\tau_{11}}.
#'
#' **Percentage of variance from random effects:**
#' \eqn{\% \text{random} = \frac{\tau_{11}}{\tau_{11} + \text{Var}(\hat{\beta}_{\text{simple slope}})} \times 100}
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(
#'   y   = rnorm(200), x = rnorm(200),
#'   m   = rep(rnorm(20), each = 10),
#'   grp = factor(rep(1:20, each = 10))
#' )
#' dat$y <- dat$y + dat$x * dat$m
#' mod <- lme4::lmer(y ~ x * m + (1 + x | grp), data = dat,
#'                   control = lme4::lmerControl(optimizer = "bobyqa"))
#' vd <- mlm_variance_decomp(mod, pred = "x", modx = "m")
#' print(vd)
#'
#' @export
mlm_variance_decomp <- function(model,
                                pred,
                                modx,
                                modx.values = c("mean-sd", "quartiles",
                                                "tertiles", "custom"),
                                at         = NULL,
                                conf.level = 0.95) {

  .check_lmer(model)
  .validate_terms(model, pred, modx)
  modx.values <- match.arg(modx.values)

  mf        <- model@frame
  modx_x    <- mf[[modx]]
  modx_vals <- .pick_modx_values(modx_x, modx.values = modx.values, at = at)

  # --- Extract random-slope variance (tau11) ---------------------------------------------------
  vc          <- lme4::VarCorr(model)
  grp_name    <- names(lme4::getME(model, "flist"))[1]
  re_mat      <- as.matrix(vc[[grp_name]])
  re_nms      <- rownames(re_mat)

  has_random_slope <- pred %in% re_nms
  tau11 <- if (has_random_slope) re_mat[pred, pred] else 0

  # --- Fixed-effect simple slopes ---------------------------------------------------------------------------
  df_resid <- get_residual_df(model)
  alpha    <- 1 - conf.level
  t_crit   <- stats::qt(1 - alpha / 2, df = df_resid)

  rows <- lapply(modx_vals, function(w) {

    ss        <- .simple_slope_linear(model, pred, modx, w,
                                      conf.level = conf.level)
    var_fixed <- ss$se^2

    # Total variance = fixed uncertainty + random-slope variance
    var_total <- var_fixed + tau11
    se_total  <- sqrt(var_total)
    pct_rand  <- if (var_total > 0) 100 * tau11 / var_total else 0

    # Prediction interval for slope in a NEW cluster
    pi_lo <- ss$slope - t_crit * se_total
    pi_hi <- ss$slope + t_crit * se_total

    data.frame(
      modx_value  = w,
      slope       = ss$slope,
      se_fixed    = ss$se,
      tau11_sd    = sqrt(tau11),
      se_total    = se_total,
      ci_lower    = ss$ci_lower,
      ci_upper    = ss$ci_upper,
      pi_lower    = pi_lo,
      pi_upper    = pi_hi,
      pct_random  = pct_rand,
      stringsAsFactors = FALSE
    )
  })

  decomp_df <- do.call(rbind, rows)
  rownames(decomp_df) <- NULL

  structure(
    list(
      decomp            = decomp_df,
      tau11             = tau11,
      tau11_sd          = sqrt(tau11),
      has_random_slope  = has_random_slope,
      pred              = pred,
      modx              = modx,
      modx.values       = modx.values,
      conf.level        = conf.level,
      grp_name          = grp_name
    ),
    class = "mlm_variance_decomp"
  )
}

#' @export
print.mlm_variance_decomp <- function(x, digits = 3, ...) {

  cat("\n========================================\n")
  cat("  Slope Variance Decomposition\n")
  cat("========================================\n")
  cat("Focal predictor :", x$pred, "\n")
  cat("Moderator       :", x$modx, "\n")
  cat("Cluster grouping:", x$grp_name, "\n")
  cat("Confidence level:", x$conf.level, "\n\n")

  if (!x$has_random_slope) {
    cat("NOTE: No random slope for '", x$pred, "' found in this model.\n",
        "Prediction intervals equal confidence intervals.\n",
        "Consider fitting (1 + ", x$pred, " | cluster) for full decomposition.\n\n",
        sep = "")
  } else {
    cat(sprintf("Random-slope SD (tau11):  %.3f\n", x$tau11_sd))
    cat(sprintf("Random-slope Var (tau11): %.3f\n\n", x$tau11))
  }

  cat("--- Per-moderator-value breakdown ---\n\n")

  df <- x$decomp
  cat(sprintf("  %-10s  %8s  %10s  %10s  %10s  %9s  %9s  %9s\n",
              x$modx, "slope", "SE(fixed)", "SE(total)",
              "% random", "CI lower", "CI upper", "PI lower"))
  cat(sprintf("  %-10s  %8s  %10s  %10s  %10s  %9s  %9s  %9s\n",
              "", "", "", "", "(new clust)", "", "", "PI upper"))
  cat(strrep("-", 88), "\n")

  for (i in seq_len(nrow(df))) {
    cat(sprintf(
      "  %-10.3f  %8.3f  %10.3f  %10.3f  %10.1f  %9.3f  %9.3f\n",
      df$modx_value[i], df$slope[i], df$se_fixed[i],
      df$se_total[i], df$pct_random[i],
      df$ci_lower[i], df$ci_upper[i]
    ))
    cat(sprintf(
      "  %-10s  %8s  %10s  %10s  %10s  %9.3f  %9.3f\n",
      "", "", "", "", "",
      df$pi_lower[i], df$pi_upper[i]
    ))
  }

  if (x$has_random_slope && mean(x$decomp$pct_random) > 30) {
    cat("\nINTERPRETATION: Random-slope variance accounts for >30% of total\n")
    cat("slope uncertainty. The moderator does not fully explain between-\n")
    cat("cluster heterogeneity in the '", x$pred, "' slope.\n", sep = "")
  }

  cat("\n")
  invisible(x)
}

#' Plot the variance decomposition of simple slopes
#'
#' Shows the simple slope at each moderator value as a point with two
#' interval layers: an inner confidence interval (fixed-effect uncertainty
#' only) and an outer prediction interval (fixed + random-slope variance).
#' The gap between the two intervals is the contribution of random-slope
#' heterogeneity.
#'
#' @param x An `mlm_variance_decomp` object.
#' @param x_label x-axis label. Defaults to moderator name.
#' @param y_label y-axis label. Defaults to "Simple slope of pred".
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#' @export
plot.mlm_variance_decomp <- function(x,
                                     x_label = NULL,
                                     y_label = NULL,
                                     ...) {

  if (is.null(x_label)) x_label <- x$modx
  if (is.null(y_label)) y_label <- paste0("Simple slope of '", x$pred, "'")

  df <- x$decomp
  df[[x$modx]] <- df$modx_value   # for aes mapping

  p <- ggplot2::ggplot(df,
         ggplot2::aes(x = modx_value, y = slope)) +

    # Outer: prediction interval (fixed + random)
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pi_lower, ymax = pi_upper),
      width = 0.06, color = "#7B2D8B", linewidth = 0.8,
      linetype = "dashed"
    ) +

    # Inner: confidence interval (fixed only)
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.04, color = "#2166AC", linewidth = 1.2
    ) +

    # Slope estimate
    ggplot2::geom_point(size = 3.5, color = "#2166AC") +

    # Zero reference
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "grey50", linewidth = 0.6) +

    ggplot2::labs(
      x        = x_label,
      y        = y_label,
      title    = "Simple Slope Variance Decomposition",
      subtitle = paste0(
        "Inner bars: fixed-effect CI  |  ",
        "Outer bars: prediction interval (adds \u03C411 = ",
        round(x$tau11_sd, 3), ")"
      )
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "grey92"),
      plot.subtitle    = ggplot2::element_text(color = "grey40", size = 10)
    )

  p
}
