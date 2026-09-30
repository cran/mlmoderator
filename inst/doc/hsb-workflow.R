## ----setup, include = FALSE---------------------------------------------------
has_mlmrev <- requireNamespace("mlmRev", quietly = TRUE)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6.5,
                      fig.height = 4, eval = has_mlmrev)

## ----eval = !has_mlmrev, echo = FALSE, results = "asis"-----------------------
# cat("*This vignette uses the `Hsb82` data from the 'mlmRev' package;",
#     "install it to run the code.*")

## ----model--------------------------------------------------------------------
library(mlmoderator)
library(lme4)
data("Hsb82", package = "mlmRev")

fit <- lmer(mAch ~ cses * meanses + cses * sector + (1 + cses | school),
            data = Hsb82)

## ----summary------------------------------------------------------------------
mlm_summary(fit, pred = "cses", modx = "meanses")

## ----df-compare---------------------------------------------------------------
sapply(c("satterthwaite", "kenward-roger", "between", "residual"),
       function(m) {
         r <- mlm_summary(fit, "cses", "meanses", jn = FALSE,
                          df_method = m)$interaction
         round(c(SE = r$se, df = r$df, p = r$p), 5)
       })

## ----plot---------------------------------------------------------------------
mlm_plot(fit, pred = "cses", modx = "meanses",
         x_label = "Student SES (school-centred)",
         y_label = "Mathematics achievement",
         legend_title = "School mean SES")

## ----jn-----------------------------------------------------------------------
plot(mlm_jn(fit, pred = "cses", modx = "meanses"))

## ----decomp-------------------------------------------------------------------
vd <- mlm_variance_decomp(fit, pred = "cses", modx = "meanses")
vd
plot(vd)

## ----loco---------------------------------------------------------------------
sens <- mlm_sensitivity(fit, pred = "cses", modx = "meanses",
                        df_method = "between")
sens
plot(sens)

