## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 7,
  fig.height = 4.5,
  fig.align  = "center"
)

## ----load---------------------------------------------------------------------
library(mlmoderator)
library(lme4)

data(school_data)
head(school_data)
str(school_data)

## ----center-------------------------------------------------------------------
# Group-mean center SES within school
dat <- mlm_center(school_data,
                  vars    = "ses",
                  cluster = "school",
                  type    = "group")

# Grand-mean center school climate (level-2 variable)
dat <- mlm_center(dat,
                  vars = "climate",
                  type = "grand")

head(dat[, c("school", "ses", "ses_c", "climate", "climate_c")])

## ----model--------------------------------------------------------------------
mod <- lmer(math ~ ses_c * climate_c + gender + (1 + ses_c | school),
            data = dat, REML = TRUE)
summary(mod)

## ----probe--------------------------------------------------------------------
probe_out <- mlm_probe(mod, pred = "ses_c", modx = "climate_c")
probe_out

## ----probe-quartiles----------------------------------------------------------
mlm_probe(mod, pred = "ses_c", modx = "climate_c",
          modx.values = "quartiles")

## ----jn-----------------------------------------------------------------------
jn_out <- mlm_jn(mod, pred = "ses_c", modx = "climate_c")
jn_out

## ----plot-basic---------------------------------------------------------------
mlm_plot(mod, pred = "ses_c", modx = "climate_c")

## ----plot-custom--------------------------------------------------------------
mlm_plot(mod,
         pred         = "ses_c",
         modx         = "climate_c",
         points       = TRUE,
         point_alpha  = 0.2,
         x_label      = "Student SES (group-mean centred)",
         y_label      = "Mathematics Achievement",
         legend_title = "School Climate")

## ----summary------------------------------------------------------------------
mlm_summary(mod, pred = "ses_c", modx = "climate_c")

