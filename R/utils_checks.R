# Internal validation helpers

#' Validate that pred, modx, and their interaction are in the model
#' @noRd
.validate_terms <- function(model, pred, modx) {
  fe <- names(lme4::fixef(model))
  if (!pred %in% fe) {
    rlang::abort(paste0("Predictor '", pred, "' not found in fixed effects."))
  }
  if (!modx %in% fe) {
    rlang::abort(paste0("Moderator '", modx, "' not found in fixed effects."))
  }
  interaction_terms <- c(
    paste0(pred, ":", modx),
    paste0(modx, ":", pred)
  )
  if (!any(interaction_terms %in% fe)) {
    rlang::abort(paste0(
      "No interaction term found between '", pred, "' and '", modx, "'. ",
      "Make sure your model includes ", pred, " * ", modx, " or ",
      pred, ":", modx, "."
    ))
  }
  invisible(TRUE)
}

#' Check that model is an lmerMod object
#' @noRd
.check_lmer <- function(model) {
  if (!inherits(model, "lmerMod")) {
    rlang::abort("'model' must be an lmerMod object fitted with lme4::lmer().")
  }
  invisible(TRUE)
}

#' Get interaction term name from fixed effects
#' @noRd
.get_interaction_term <- function(model, pred, modx) {
  fe <- names(lme4::fixef(model))
  term1 <- paste0(pred, ":", modx)
  term2 <- paste0(modx, ":", pred)
  if (term1 %in% fe) return(term1)
  if (term2 %in% fe) return(term2)
  rlang::abort("Interaction term not found.")
}
