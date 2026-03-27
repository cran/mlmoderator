#' Johnson---Neyman interval for multilevel two-way interactions
#'
#' Computes the Johnson---Neyman (JN) interval: the region(s) of the moderator
#' (`modx`) where the simple slope of `pred` transitions between statistical
#' significance and non-significance. Useful for identifying *exactly* at which
#' moderator values an effect becomes significant.
#'
#' @param model An `lmerMod` object with a two-way interaction between `pred`
#'   and `modx` in the fixed-effects structure.
#' @param pred Character scalar. Focal predictor name.
#' @param modx Character scalar. Moderator name.
#' @param alpha Significance level. Default `0.05`.
#' @param modx.range Numeric vector of length 2 giving the range over which to
#'   evaluate the slope. Defaults to the observed range of `modx`.
#' @param grid Integer. Number of points at which to evaluate the slope across
#'   the moderator range. Default `200`.
#'
#' @return An object of class `mlm_jn` with components:
#'   * `jn_bounds`: numeric vector of moderator values where the slope crosses
#'     the significance threshold. `NA` if no finite region exists.
#'   * `slopes_df`: data frame of slope estimates and significance across the
#'     grid.
#'   * `pred`, `modx`, `alpha`, `modx.range`.
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(
#'   y   = rnorm(200), x = rnorm(200),
#'   m   = rep(rnorm(20), each = 10),
#'   grp = factor(rep(1:20, each = 10))
#' )
#' dat$y <- dat$y + dat$x * dat$m
#' mod <- lme4::lmer(y ~ x * m + (1 | grp), data = dat,
#'                   control = lme4::lmerControl(optimizer = "bobyqa"))
#' jn <- mlm_jn(mod, pred = "x", modx = "m")
#' print(jn)
#'
#' @export
mlm_jn <- function(model,
                   pred,
                   modx,
                   alpha       = 0.05,
                   modx.range  = NULL,
                   grid        = 200L) {

  .check_lmer(model)
  .validate_terms(model, pred, modx)

  mf     <- model@frame
  modx_x <- mf[[modx]]

  if (is.null(modx.range)) {
    modx.range <- range(modx_x, na.rm = TRUE)
  }

  modx_grid <- seq(modx.range[1], modx.range[2], length.out = grid)

  slopes_list <- lapply(modx_grid, function(w) {
    ss <- .simple_slope_linear(model, pred, modx, w, conf.level = 1 - alpha)
    data.frame(
      modx_value = w,
      slope      = ss$slope,
      se         = ss$se,
      t          = ss$t,
      df         = ss$df,
      p          = ss$p,
      sig        = ss$p < alpha,
      ci_lower   = ss$ci_lower,
      ci_upper   = ss$ci_upper,
      stringsAsFactors = FALSE
    )
  })
  slopes_df <- do.call(rbind, slopes_list)
  rownames(slopes_df) <- NULL

  # Detect sign changes in (p - alpha) to locate JN boundaries
  p_centered <- slopes_df$p - alpha
  sign_changes <- which(diff(sign(p_centered)) != 0)

  jn_bounds <- if (length(sign_changes) == 0) {
    NA_real_
  } else {
    sapply(sign_changes, function(i) {
      # Linear interpolation between grid points
      x0 <- slopes_df$modx_value[i]
      x1 <- slopes_df$modx_value[i + 1]
      y0 <- p_centered[i]
      y1 <- p_centered[i + 1]
      x0 + (0 - y0) * (x1 - x0) / (y1 - y0)
    })
  }

  structure(
    list(
      jn_bounds  = jn_bounds,
      slopes_df  = slopes_df,
      pred       = pred,
      modx       = modx,
      alpha      = alpha,
      modx.range = modx.range
    ),
    class = "mlm_jn"
  )
}

#' Plot a Johnson-Neyman interval from a multilevel model
#'
#' Creates a ggplot2 figure showing the simple slope of `pred` across the
#' full range of `modx`, with shading indicating regions of significance.
#' A vertical dashed line marks each Johnson-Neyman boundary.
#'
#' @param x An `mlm_jn` object from `mlm_jn()`.
#' @param x_label Label for the x-axis. Defaults to the moderator name.
#' @param y_label Label for the y-axis. Defaults to "Simple slope of pred".
#' @param sig_color Fill colour for the significant region. Default `"#2166AC"`.
#' @param nonsig_color Fill colour for the non-significant region. Default `"#D6604D"`.
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(
#'   y   = rnorm(200), x = rnorm(200),
#'   m   = rep(rnorm(20), each = 10),
#'   grp = factor(rep(1:20, each = 10))
#' )
#' dat$y <- dat$y + dat$x * dat$m
#' mod <- lme4::lmer(y ~ x * m + (1 | grp), data = dat,
#'                   control = lme4::lmerControl(optimizer = "bobyqa"))
#' jn <- mlm_jn(mod, pred = "x", modx = "m")
#' plot(jn)
#'
#' @export
plot.mlm_jn <- function(x,
                        x_label    = NULL,
                        y_label    = NULL,
                        sig_color  = "#2166AC",
                        nonsig_color = "#D6604D",
                        ...) {
  df <- x$slopes_df

  if (is.null(x_label)) x_label <- x$modx
  if (is.null(y_label)) y_label <- paste0("Simple slope of '", x$pred, "'")

  # Significance as factor for fill
  df$region <- ifelse(df$sig, "Significant", "Non-significant")
  df$region <- factor(df$region, levels = c("Significant", "Non-significant"))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = modx_value, y = slope)) +
    # Shaded ribbon for CI
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper, fill = region),
      alpha = 0.25
    ) +
    # Slope line coloured by significance
    ggplot2::geom_line(ggplot2::aes(color = region), linewidth = 1.1) +
    # Zero reference line
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "grey40", linewidth = 0.7) +
    # JN boundary lines
    {
      if (!all(is.na(x$jn_bounds))) {
        lapply(x$jn_bounds, function(b) {
          ggplot2::geom_vline(xintercept = b, linetype = "dashed",
                              color = "grey20", linewidth = 0.8)
        })
      }
    } +
    # JN boundary labels
    {
      if (!all(is.na(x$jn_bounds))) {
        lapply(x$jn_bounds, function(b) {
          ggplot2::annotate("text",
            x = b, y = max(df$ci_upper, na.rm = TRUE),
            label = paste0(x$modx, " = ", round(b, 2)),
            hjust = -0.05, vjust = 1, size = 3.5, color = "grey20"
          )
        })
      }
    } +
    ggplot2::scale_color_manual(
      values = c("Significant" = sig_color, "Non-significant" = nonsig_color),
      name   = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = c("Significant" = sig_color, "Non-significant" = nonsig_color),
      name   = NULL
    ) +
    ggplot2::labs(x = x_label, y = y_label,
                  title = "Johnson-Neyman Plot",
                  subtitle = paste0("Regions where slope of '", x$pred,
                                    "' is significant (p < ", x$alpha, ")")) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      legend.position  = "bottom",
      panel.grid.major = ggplot2::element_line(color = "grey92"),
      plot.subtitle    = ggplot2::element_text(color = "grey40", size = 11)
    )

  p
}

#' @export
print.mlm_jn <- function(x, digits = 3, ...) {
  cat("\n--- Johnson-Neyman Interval: mlm_jn ---\n")
  cat("Focal predictor :", x$pred, "\n")
  cat("Moderator       :", x$modx, "\n")
  cat("Alpha           :", x$alpha, "\n")
  cat("Moderator range :", round(x$modx.range[1], digits), "to",
      round(x$modx.range[2], digits), "\n\n")

  if (all(is.na(x$jn_bounds))) {
    cat("No Johnson-Neyman boundary found in the observed moderator range.\n")
    # Indicate whether slope is universally sig or non-sig
    pct_sig <- mean(x$slopes_df$sig, na.rm = TRUE)
    if (pct_sig > 0.99) {
      cat("The simple slope of", x$pred, "is significant across the entire",
          "moderator range.\n")
    } else if (pct_sig < 0.01) {
      cat("The simple slope of", x$pred, "is non-significant across the entire",
          "moderator range.\n")
    }
  } else {
    cat("Johnson-Neyman boundary/boundaries (", x$modx, "):\n")
    for (b in x$jn_bounds) {
      cat("  ", round(b, digits), "\n")
    }
    cat("\nInterpretation: The simple slope of '", x$pred, "' is statistically\n",
        "significant (p < ", x$alpha, ") for values of '", x$modx, "'\n", sep = "")

    df <- x$slopes_df
    sig_range  <- range(df$modx_value[df$sig], na.rm = TRUE)
    nsig_range <- range(df$modx_value[!df$sig], na.rm = TRUE)

    if (!any(df$sig)) {
      cat("  outside the observed range (effect is non-significant throughout).\n")
    } else if (!any(!df$sig)) {
      cat("  across the entire observed range.\n")
    } else {
      cat("  between ", round(sig_range[1], digits),
          " and ", round(sig_range[2], digits), ".\n", sep = "")
    }
  }
  cat("\n")
  invisible(x)
}
