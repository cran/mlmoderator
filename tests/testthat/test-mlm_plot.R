## Tests for mlm_plot()

make_test_model <- function() {
  set.seed(99)
  n_grp <- 20
  n_obs <- 15
  N     <- n_grp * n_obs
  grp   <- rep(seq_len(n_grp), each = n_obs)
  x     <- rnorm(N)
  m     <- rep(rnorm(n_grp), each = n_obs)
  u0    <- rep(rnorm(n_grp, sd = 1), each = n_obs)
  y     <- 2 + 1.5 * x + 0.8 * m + 0.6 * x * m + u0 + rnorm(N)
  dat   <- data.frame(y = y, x = x, m = m, grp = factor(grp))
  suppressMessages(
    lme4::lmer(y ~ x * m + (1 | grp), data = dat, REML = FALSE)
  )
}

test_that("mlm_plot returns a ggplot object", {
  mod <- make_test_model()
  p   <- mlm_plot(mod, pred = "x", modx = "m")
  expect_s3_class(p, "gg")
})

test_that("mlm_plot works without confidence bands", {
  mod <- make_test_model()
  p   <- mlm_plot(mod, pred = "x", modx = "m", interval = FALSE)
  expect_s3_class(p, "gg")
})

test_that("mlm_plot works with raw data points", {
  mod <- make_test_model()
  p   <- mlm_plot(mod, pred = "x", modx = "m", points = TRUE)
  expect_s3_class(p, "gg")
})

test_that("mlm_plot works with quartile values", {
  mod <- make_test_model()
  p   <- mlm_plot(mod, pred = "x", modx = "m", modx.values = "quartiles")
  expect_s3_class(p, "gg")
})

test_that("mlm_plot works with custom at values", {
  mod <- make_test_model()
  p   <- mlm_plot(mod, pred = "x", modx = "m",
                  modx.values = "custom", at = c(-1, 0, 1))
  expect_s3_class(p, "gg")
})

test_that("mlm_plot custom labels work", {
  mod <- make_test_model()
  p   <- mlm_plot(mod, pred = "x", modx = "m",
                  x_label = "SES", y_label = "Math", legend_title = "Climate")
  expect_equal(p$labels$x, "SES")
  expect_equal(p$labels$y, "Math")
})

test_that("mlm_plot errors on non-lmer model", {
  expect_error(mlm_plot(list(), pred = "x", modx = "m"), "lmerMod")
})
