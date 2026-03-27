## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.width  = 7,
  fig.height = 4.5,
  fig.align  = "center"
)

## ----setup-pkg----------------------------------------------------------------
library(mlmoderator)
library(lme4)

data(school_data)

## ----centering----------------------------------------------------------------
dat <- mlm_center(school_data,
                  vars    = "ses",
                  cluster = "school",
                  type    = "group")   # within-school SES

dat <- mlm_center(dat,
                  vars = "climate",
                  type = "grand")      # grand-mean-centred climate

# Check
cat("Within-school SES mean (should be ~0):", round(mean(dat$ses_c), 4), "\n")
cat("Grand-mean climate (should be ~0):",     round(mean(dat$climate_c), 4), "\n")

## ----fit----------------------------------------------------------------------
mod <- lmer(
  math ~ ses_c * climate_c + gender + (1 + ses_c | school),
  data = dat,
  REML = TRUE
)
summary(mod)

## ----probe--------------------------------------------------------------------
probe <- mlm_probe(mod,
                   pred        = "ses_c",
                   modx        = "climate_c",
                   modx.values = "mean-sd",
                   conf.level  = 0.95)
probe

## ----jn-----------------------------------------------------------------------
jn <- mlm_jn(mod, pred = "ses_c", modx = "climate_c")
jn

## ----plot---------------------------------------------------------------------
mlm_plot(
  mod,
  pred         = "ses_c",
  modx         = "climate_c",
  modx.values  = "mean-sd",
  interval     = TRUE,
  x_label      = "Student SES (group-mean centred)",
  y_label      = "Mathematics Achievement",
  legend_title = "School Climate"
)

## ----summary------------------------------------------------------------------
mlm_summary(mod,
            pred        = "ses_c",
            modx        = "climate_c",
            modx.values = "mean-sd",
            jn          = TRUE)

