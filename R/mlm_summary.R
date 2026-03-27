#' Summary table for a multilevel moderation effect
#'
#' Returns a consolidated summary of the moderation effect: the focal
#' interaction coefficient, simple slopes at selected moderator values, and
#' (optionally) the Johnson---Neyman interval. Designed for quick reporting and
#' results sections.
#'
#' @param model An `lmerMod` object with a two-way interaction between `pred`
#'   and `modx`.
#' @param pred Character scalar. Focal predictor name.
#' @param modx Character scalar. Moderator name.
#' @param modx.values Moderator value strategy. See `mlm_probe()`.
#' @param at Optional numeric vector of custom moderator values.
#' @param conf.level Confidence level. Default `0.95`.
#' @param jn Logical. Include Johnson-Neyman region? Default `TRUE`.
#' @param alpha Alpha for JN interval. Default `0.05`.
#'
#' @return An object of class `mlm_summary` (a list) with components:
#'   * `interaction`: one-row data frame for the interaction term.
#'   * `simple_slopes`: data frame from `mlm_probe()`.
#'   * `jn`: output of `mlm_jn()` (or `NULL` if `jn = FALSE`).
#'   * Other metadata.
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
#' mlm_summary(mod, pred = "x", modx = "m")
#'
#' @export
mlm_summary <- function(model,
                        pred,
                        modx,
                        modx.values = c("mean-sd", "quartiles", "tertiles", "custom"),
                        at          = NULL,
                        conf.level  = 0.95,
                        jn          = TRUE,
                        alpha       = 0.05) {

  .check_lmer(model)
  .validate_terms(model, pred, modx)
  modx.values <- match.arg(modx.values)

  int_term <- .get_interaction_term(model, pred, modx)
  fe       <- lme4::fixef(model)
  vcv      <- .extract_vcov(model)

  b_int  <- fe[int_term]
  se_int <- sqrt(vcv[int_term, int_term])
  df_r   <- get_residual_df(model)
  t_int  <- b_int / se_int
  p_int  <- 2 * stats::pt(abs(t_int), df = df_r, lower.tail = FALSE)
  t_crit <- stats::qt(1 - (1 - conf.level) / 2, df = df_r)

  interaction_row <- data.frame(
    term     = int_term,
    estimate = b_int,
    se       = se_int,
    t        = t_int,
    df       = df_r,
    p        = p_int,
    ci_lower = b_int - t_crit * se_int,
    ci_upper = b_int + t_crit * se_int,
    stringsAsFactors = FALSE
  )

  probe_out <- mlm_probe(model, pred = pred, modx = modx,
                         modx.values = modx.values, at = at,
                         conf.level = conf.level)

  jn_out <- NULL
  if (jn) {
    jn_out <- mlm_jn(model, pred = pred, modx = modx, alpha = alpha)
  }

  structure(
    list(
      interaction   = interaction_row,
      simple_slopes = probe_out$slopes,
      jn            = jn_out,
      pred          = pred,
      modx          = modx,
      modx.values   = modx.values,
      conf.level    = conf.level,
      alpha         = alpha,
      model         = model
    ),
    class = "mlm_summary"
  )
}

#' @export
print.mlm_summary <- function(x, digits = 3, ...) {
  cat("\n========================================\n")
  cat("  Multilevel Moderation Summary\n")
  cat("========================================\n")
  cat("Focal predictor :", x$pred, "\n")
  cat("Moderator       :", x$modx, "\n")
  cat("Confidence level:", x$conf.level, "\n")

  cat("\n--- Interaction Term ---\n")
  ir <- x$interaction
  cat(sprintf(
    "  %-28s  b = %7.3f  SE = %6.3f  t(%g) = %6.3f  p = %s  [%s, %s]\n",
    ir$term,
    ir$estimate, ir$se, round(ir$df, 0), ir$t,
    format_pval(ir$p),
    round(ir$ci_lower, digits),
    round(ir$ci_upper, digits)
  ))

  cat("\n--- Simple Slopes ---\n")
  df <- x$simple_slopes
  header <- sprintf("  %-10s  %8s  %8s  %8s  %8s  %8s  %8s  %8s",
                    x$modx, "slope", "se", "t", "df", "p",
                    paste0(round(x$conf.level * 100), "% CI lo"),
                    paste0(round(x$conf.level * 100), "% CI hi"))
  cat(header, "\n")
  cat(strrep("-", nchar(header)), "\n")
  for (i in seq_len(nrow(df))) {
    cat(sprintf("  %-10.3f  %8.3f  %8.3f  %8.3f  %8.1f  %8s  %8.3f  %8.3f\n",
                df$modx_value[i], df$slope[i], df$se[i], df$t[i],
                df$df[i], format_pval(df$p[i]),
                df$ci_lower[i], df$ci_upper[i]))
  }

  if (!is.null(x$jn)) {
    cat("\n--- Johnson-Neyman Region ---\n")
    jb <- x$jn$jn_bounds
    if (all(is.na(jb))) {
      pct_sig <- mean(x$jn$slopes_df$sig, na.rm = TRUE)
      if (pct_sig > 0.99) {
        cat("  Slope of '", x$pred, "' is significant across the entire range of '",
            x$modx, "'.\n", sep = "")
      } else if (pct_sig < 0.01) {
        cat("  Slope of '", x$pred, "' is non-significant across the entire range of '",
            x$modx, "'.\n", sep = "")
      } else {
        cat("  No clear JN boundary found.\n")
      }
    } else {
      cat("  JN boundary/boundaries for '", x$modx, "': ",
          paste(round(jb, digits), collapse = ", "), "\n", sep = "")
    }
  }

  cat("\n========================================\n\n")
  invisible(x)
}
