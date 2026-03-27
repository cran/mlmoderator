# mlmoderator

<!-- badges: start -->
[![R-CMD-check](https://github.com/causalfragility-lab/mlmoderator/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/causalfragility-lab/mlmoderator/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

**Robustness diagnostics and visualization tools for cross-level interaction
effects in multilevel models.**

`mlmoderator` provides a unified workflow for probing, plotting, and assessing
the robustness of two-way interaction effects in mixed-effects models fitted
with `lme4::lmer()`. It is designed for researchers in education, psychology,
biostatistics, epidemiology, organizational science, and any field where
outcomes are clustered within higher-level units.

---

## Installation

```r
# Development version from GitHub
remotes::install_github("causalfragility-lab/mlmoderator")
```

---

## Core functions

| Function | Description |
|---|---|
| `mlm_center()` | Grand- or group-mean center variables |
| `mlm_probe()` | Simple slopes at selected moderator values |
| `mlm_jn()` | Johnson-Neyman significance regions |
| `mlm_plot()` | Publication-ready interaction plot |
| `mlm_surface()` | Contour plot of predicted outcomes over pred x modx space |
| `mlm_summary()` | Consolidated moderation report |
| `mlm_variance_decomp()` | Decompose slope uncertainty into fixed vs random components |
| `mlm_sensitivity()` | ICC-shift robustness + leave-one-cluster-out (LOCO) diagnostics |

---

## Usage examples

### Education: SES x School Climate

```r
library(mlmoderator)
library(lme4)

data(school_data)

# Center variables
dat <- mlm_center(school_data, vars = "ses",     cluster = "school", type = "group")
dat <- mlm_center(dat,         vars = "climate",                      type = "grand")

# Fit model
mod <- lmer(math ~ ses_c * climate_c + gender + (1 + ses_c | school),
            data = dat, REML = TRUE)

# Simple slopes at -1 SD, Mean, +1 SD of school climate
mlm_probe(mod, pred = "ses_c", modx = "climate_c")
#>  climate_c  slope    se      t   df      p ci_lower ci_upper
#>     -1.036  0.944 0.145  6.504 2995 < .001    0.659    1.228
#>      0.000  1.426 0.101 14.075 2995 < .001    1.228    1.625
#>      1.036  1.909 0.144 13.277 2995 < .001    1.627    2.191

# Johnson-Neyman interval
jn <- mlm_jn(mod, pred = "ses_c", modx = "climate_c")
print(jn)
plot(jn)

# Interaction plot
mlm_plot(mod,
         pred         = "ses_c",
         modx         = "climate_c",
         x_label      = "Student SES (group-mean centred)",
         y_label      = "Mathematics Achievement",
         legend_title = "School Climate")

# Contour surface - predicted math over the full ses x climate space
mlm_surface(mod, pred = "ses_c", modx = "climate_c")

# Full summary report
mlm_summary(mod, pred = "ses_c", modx = "climate_c")

# Variance decomposition: fixed-effect CI vs prediction interval
vd <- mlm_variance_decomp(mod, pred = "ses_c", modx = "climate_c")
print(vd)
plot(vd)

# Robustness diagnostics: ICC shift + LOCO
sens <- mlm_sensitivity(mod, pred = "ses_c", modx = "climate_c")
print(sens)
plot(sens)
```

---

### Biostatistics / Clinical Trials: Treatment x Baseline Severity

Does the effect of a treatment vary depending on how severe the patient's
condition was at baseline? Patients are nested within hospitals.

```r
# mod <- lmer(
#   recovery ~ treatment * baseline_severity + age + comorbidity +
#     (1 + treatment | hospital),
#   data = trial_data
# )

# Where does the treatment effect become significant?
# mlm_jn(mod, pred = "treatment", modx = "baseline_severity")

# How stable is this across hospitals?
# sens <- mlm_sensitivity(mod, pred = "treatment", modx = "baseline_severity")
# plot(sens)  # LOCO panel shows hospital-level influence

# How much of the treatment slope variation is unexplained by severity?
# vd <- mlm_variance_decomp(mod, pred = "treatment", modx = "baseline_severity")
# print(vd)   # % random tells you how much tau11 remains after conditioning on severity
```

---

### Epidemiology: Individual SES x Neighborhood Deprivation

Does individual socioeconomic status moderate the effect of neighborhood
deprivation on BMI? Individuals nested within neighborhoods.

```r
# mod <- lmer(
#   bmi ~ neighborhood_deprivation * individual_ses + age + smoking +
#     (1 + individual_ses | neighborhood),
#   data = health_data
# )

# Full predicted BMI surface over the deprivation x SES space
# mlm_surface(mod,
#             pred    = "individual_ses",
#             modx    = "neighborhood_deprivation",
#             x_label = "Individual SES",
#             y_label = "Neighborhood Deprivation Index")

# Contour rotation reveals where the cross-level interaction is strongest
```

---

### Organizational Psychology: Workload x Leadership Quality

Does leadership quality buffer the effect of workload on burnout?
Employees nested within teams.

```r
# mod <- lmer(
#   burnout ~ workload * leadership + tenure + role +
#     (1 + workload | team),
#   data = org_data
# )

# At what level of leadership quality does workload stop predicting burnout?
# mlm_jn(mod, pred = "workload", modx = "leadership")

# Are results driven by specific teams?
# sens <- mlm_sensitivity(mod, pred = "workload", modx = "leadership",
#                         verbose = TRUE)
# plot(sens)  # LOCO panel flags influential teams
```

---

### Economics / Policy: Education x Regional GDP

Does regional economic development moderate the return to education in wages?
Workers nested within regions.

```r
# mod <- lmer(
#   log_wage ~ education * regional_gdp + experience + sector +
#     (1 + education | region),
#   data = labor_data
# )

# Contour surface: predicted log wage over education x regional GDP
# mlm_surface(mod,
#             pred         = "education",
#             modx         = "regional_gdp",
#             legend_title = "log(Wage)",
#             x_label      = "Years of Education",
#             y_label      = "Regional GDP per capita (standardised)")
```

---

### Neuroscience / fMRI: Age x Task Difficulty

Does task difficulty moderate the age-related decline in neural activation?
Observations nested within subjects.

```r
# mod <- lmer(
#   bold_signal ~ age * task_difficulty + motion + session +
#     (1 + age | subject),
#   data = fmri_data
# )

# Simple slopes: age effect at low, medium, high task difficulty
# mlm_probe(mod, pred = "age", modx = "task_difficulty",
#           modx.values = "tertiles")

# Variance decomp: is the age slope variable across subjects
# beyond what task difficulty explains?
# vd <- mlm_variance_decomp(mod, pred = "age", modx = "task_difficulty")
# print(vd)
```

---

## Model requirements

| Requirement | Detail |
|---|---|
| Model class | `lmerMod` from `lme4::lmer()` |
| Interaction | Two-way: `pred * modx` or `pred : modx` in fixed effects |
| Outcome | Continuous (linear mixed model) |
| Clustering | Any single grouping variable |
| Predictors | Continuous x continuous, or continuous x binary |

---

## Scope of sensitivity analysis

`mlm_sensitivity()` provides **robustness diagnostics**, not a full causal
sensitivity analysis:

- **ICC-shift**: how does the interaction SE and JN boundary change if the
  intraclass correlation were different?
- **LOCO**: does the interaction hold when each cluster is dropped in turn?
  Which clusters are most influential?

A level-2 omitted variable bound specifically designed for cross-level
interactions is under development as a separate methodological contribution.
Standard OLS-based sensitivity tools (E-value, ITCV) are not appropriate
for multilevel interaction effects and are intentionally excluded.

---

## Citation

If you use `mlmoderator` in your research, please cite:

```
Hait, S. (2026). mlmoderator: Robustness Diagnostics and Visualization
for Cross-Level Interaction Effects in Multilevel Models.
R package version 0.2.0.
https://github.com/causalfragility-lab/mlmoderator
```

Or in BibTeX:

```bibtex
@Manual{mlmoderator,
  title  = {mlmoderator: Robustness Diagnostics and Visualization
            for Cross-Level Interaction Effects in Multilevel Models},
  author = {Subir Hait},
  year   = {2026},
  note   = {R package version 0.2.0},
  url    = {https://github.com/causalfragility-lab/mlmoderator}
}
```

---

## License

MIT
