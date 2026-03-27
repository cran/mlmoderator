# Internal variance-covariance helpers

#' Extract fixed-effect variance-covariance matrix
#' @noRd
.extract_vcov <- function(model) {
  as.matrix(stats::vcov(model))
}

#' Compute simple slope and SE using the delta method
#'
#' For a model: y ~ b1*pred + b2*modx + b3*(pred:modx) + ...
#' The simple slope of pred at modx = w is:
#'   slope = b1 + b3 * w
#'   var(slope) = var(b1) + w^2 * var(b3) + 2 * w * cov(b1, b3)
#'
#' @noRd
.simple_slope_linear <- function(model, pred, modx, modx_val, conf.level = 0.95) {
  fe      <- lme4::fixef(model)
  vcv     <- .extract_vcov(model)
  int_term <- .get_interaction_term(model, pred, modx)

  b_pred <- fe[pred]
  b_int  <- fe[int_term]

  slope  <- b_pred + b_int * modx_val

  v_pred <- vcv[pred, pred]
  v_int  <- vcv[int_term, int_term]
  cov_pi <- vcv[pred, int_term]

  var_slope <- v_pred + modx_val^2 * v_int + 2 * modx_val * cov_pi
  se_slope  <- sqrt(var_slope)

  df_resid <- get_residual_df(model)
  t_val    <- slope / se_slope
  p_val    <- 2 * stats::pt(abs(t_val), df = df_resid, lower.tail = FALSE)

  alpha <- 1 - conf.level
  t_crit <- stats::qt(1 - alpha / 2, df = df_resid)
  ci_lo  <- slope - t_crit * se_slope
  ci_hi  <- slope + t_crit * se_slope

  list(
    slope    = slope,
    se       = se_slope,
    t        = t_val,
    df       = df_resid,
    p        = p_val,
    ci_lower = ci_lo,
    ci_upper = ci_hi
  )
}

#' Approximate residual degrees of freedom (Satterthwaite-like fallback)
#' @noRd
get_residual_df <- function(model) {
  n <- nrow(model@frame)
  p <- length(lme4::fixef(model))
  max(n - p, 1)
}
