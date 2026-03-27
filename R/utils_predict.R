# Internal prediction helpers

#' Choose moderator values based on strategy
#'
#' Returns a sorted numeric vector of moderator values to probe.
#' If `at` is supplied it always takes priority over `modx.values`.
#'
#' @param x         Numeric vector of the moderator (from model data).
#' @param modx.values Strategy: `"mean-sd"`, `"quartiles"`, `"tertiles"`,
#'   or `"custom"`.
#' @param at        Custom numeric values; required when
#'   `modx.values = "custom"`.
#' @noRd
.pick_modx_values <- function(x, modx.values = "mean-sd", at = NULL) {

  # Custom values via `at` always win, regardless of modx.values
  if (!is.null(at)) return(sort(as.numeric(at)))

  modx.values <- match.arg(
    modx.values,
    choices = c("mean-sd", "quartiles", "tertiles", "custom")
  )

  switch(modx.values,
    "mean-sd" = {
      m <- mean(x, na.rm = TRUE)
      s <- stats::sd(x, na.rm = TRUE)
      c(m - s, m, m + s)
    },
    "quartiles" = {
      as.numeric(stats::quantile(x, probs = c(0.25, 0.50, 0.75), na.rm = TRUE))
    },
    "tertiles" = {
      as.numeric(stats::quantile(x, probs = c(1/3, 2/3), na.rm = TRUE))
    },
    "custom" = {
      rlang::abort(
        "Supply moderator values via the `at` argument when modx.values = 'custom'."
      )
    }
  )
}

#' Make a prediction grid for a two-way interaction
#'
#' Builds a data frame with `n_pred` evenly-spaced values of `pred` crossed
#' with each value in `modx_vals`. All other covariates are held at their
#' means (numeric) or reference level (factor). The grouping variable is set
#' to its first observed level so `predict(..., re.form = NA)` works cleanly.
#'
#' @param model     `lmerMod` object.
#' @param pred      Focal predictor name (character).
#' @param modx      Moderator name (character).
#' @param modx_vals Numeric vector of moderator values to use.
#' @param n_pred    Number of points along the predictor range. Default 100.
#' @noRd
.make_prediction_grid <- function(model, pred, modx, modx_vals, n_pred = 100) {

  mf <- model@frame
  cluster_vars <- names(lme4::getME(model, "flist"))

  pred_range <- seq(
    min(mf[[pred]], na.rm = TRUE),
    max(mf[[pred]], na.rm = TRUE),
    length.out = n_pred
  )

  grids <- lapply(modx_vals, function(w) {

    g <- mf[rep(1L, n_pred), , drop = FALSE]
    rownames(g) <- NULL

    g[[pred]] <- pred_range
    g[[modx]] <- w

    # Hold all other variables at sensible defaults
    other_vars <- setdiff(names(g), c(pred, modx, cluster_vars))
    for (v in other_vars) {
      if (is.numeric(g[[v]])) {
        g[[v]] <- mean(mf[[v]], na.rm = TRUE)
      } else if (is.factor(g[[v]])) {
        g[[v]] <- factor(levels(mf[[v]])[1], levels = levels(mf[[v]]))
      } else if (is.character(g[[v]])) {
        g[[v]] <- mf[[v]][1]
      }
    }

    g$.modx_val <- w
    g
  })

  do.call(rbind, grids)
}

#' Build clean legend labels for selected moderator values
#'
#' Labels are short and do NOT repeat the moderator name --- the legend title
#' already carries that. SD-based: "-1 SD", "Mean", "+1 SD". Quartile-based:
#' "25th pct" etc. Custom/fallback: the rounded numeric value.
#'
#' @param vals     Numeric vector of moderator values.
#' @param modx     Moderator name (passed through but used as legend title only).
#' @param strategy The `modx.values` strategy that produced `vals`.
#' @noRd
.build_modx_labels <- function(vals, modx, strategy) {

  lbls <- if (strategy == "mean-sd" && length(vals) == 3L) {
    c("-1 SD", "Mean", "+1 SD")
  } else if (strategy == "quartiles") {
    c("25th pct", "50th pct", "75th pct")
  } else if (strategy == "tertiles") {
    c("1st tertile", "2nd tertile")
  } else {
    as.character(round(vals, 2))
  }

  stats::setNames(lbls, as.character(vals))
}
