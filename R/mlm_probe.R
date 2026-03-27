#' Probe simple slopes from a multilevel interaction
#'
#' Computes simple slopes of a focal predictor (`pred`) at selected values of a
#' moderator (`modx`) from a two-level mixed-effects model fitted with
#' [lme4::lmer()]. Returns estimates, standard errors, *t*-values, *p*-values,
#' and confidence intervals in a tidy data frame.
#'
#' @param model An `lmerMod` object containing a two-way interaction between
#'   `pred` and `modx` in the fixed-effects structure.
#' @param pred Character scalar. Name of the focal predictor variable.
#' @param modx Character scalar. Name of the moderator variable.
#' @param modx.values Strategy for selecting moderator values. One of:
#'   * `"mean-sd"` (default): mean --- 1 SD, mean, mean + 1 SD.
#'   * `"quartiles"`: 25th, 50th, 75th percentiles.
#'   * `"tertiles"`: 33rd and 67th percentiles.
#'   * `"custom"`: use values supplied via `at`.
#' @param at Numeric vector of custom moderator values. Used when
#'   `modx.values = "custom"`, or to override any strategy.
#' @param conf.level Confidence level for intervals. Default `0.95`.
#'
#' @return An object of class `mlm_probe` (a list) with components:
#'   * `slopes`: a data frame with columns `modx_value`, `slope`, `se`, `t`,
#'     `df`, `p`, `ci_lower`, `ci_upper`.
#'   * `pred`, `modx`: names of the predictor and moderator.
#'   * `modx.values`: the strategy used.
#'   * `conf.level`: the confidence level.
#'   * `model`: the original model (stored for downstream use).
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
#' mlm_probe(mod, pred = "x", modx = "m")
#' mlm_probe(mod, pred = "x", modx = "m", modx.values = "quartiles")
#' mlm_probe(mod, pred = "x", modx = "m", modx.values = "custom", at = c(-1, 0, 1))
#'
#' @export
mlm_probe <- function(model,
                      pred,
                      modx,
                      modx.values = c("mean-sd", "quartiles", "tertiles", "custom"),
                      at          = NULL,
                      conf.level  = 0.95) {

  .check_lmer(model)
  .validate_terms(model, pred, modx)

  modx.values <- match.arg(modx.values)
  mf     <- model@frame
  modx_x <- mf[[modx]]

  modx_vals <- .pick_modx_values(modx_x, modx.values = modx.values, at = at)

  rows <- lapply(modx_vals, function(w) {
    ss <- .simple_slope_linear(model, pred, modx, w, conf.level = conf.level)
    data.frame(
      modx_value = w,
      slope      = ss$slope,
      se         = ss$se,
      t          = ss$t,
      df         = ss$df,
      p          = ss$p,
      ci_lower   = ss$ci_lower,
      ci_upper   = ss$ci_upper,
      stringsAsFactors = FALSE
    )
  })

  slopes_df <- do.call(rbind, rows)
  rownames(slopes_df) <- NULL

  structure(
    list(
      slopes      = slopes_df,
      pred        = pred,
      modx        = modx,
      modx.values = modx.values,
      conf.level  = conf.level,
      model       = model
    ),
    class = "mlm_probe"
  )
}

#' @export
print.mlm_probe <- function(x, digits = 3, ...) {
  cat("\n--- Simple Slopes: mlm_probe ---\n")
  cat("Focal predictor :", x$pred, "\n")
  cat("Moderator       :", x$modx, "\n")
  cat("Confidence level:", x$conf.level, "\n\n")

  df <- x$slopes
  df_print <- data.frame(
    modx_value = round(df$modx_value, digits),
    slope      = round(df$slope,      digits),
    se         = round(df$se,         digits),
    t          = round(df$t,          digits),
    df         = round(df$df,         1),
    p          = format_pval(df$p),
    ci_lower   = round(df$ci_lower,   digits),
    ci_upper   = round(df$ci_upper,   digits)
  )
  names(df_print)[1] <- x$modx

  print(df_print, row.names = FALSE)
  cat("\n")
  invisible(x)
}

#' Format p-values for printing
#' @noRd
format_pval <- function(p) {
  ifelse(p < 0.001, "< .001", formatC(p, digits = 3, format = "f"))
}
