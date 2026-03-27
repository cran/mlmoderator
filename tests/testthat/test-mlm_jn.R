## Tests for mlm_jn()

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

test_that("mlm_jn returns an mlm_jn object", {
  mod <- make_test_model()
  out <- mlm_jn(mod, pred = "x", modx = "m")
  expect_s3_class(out, "mlm_jn")
})

test_that("mlm_jn slopes_df has correct columns", {
  mod <- make_test_model()
  out <- mlm_jn(mod, pred = "x", modx = "m")
  expect_named(out$slopes_df,
    c("modx_value", "slope", "se", "t", "df", "p", "sig", "ci_lower", "ci_upper"))
})

test_that("mlm_jn slopes_df has grid rows", {
  mod <- make_test_model()
  out <- mlm_jn(mod, pred = "x", modx = "m", grid = 100L)
  expect_equal(nrow(out$slopes_df), 100L)
})

test_that("mlm_jn respects custom modx.range", {
  mod <- make_test_model()
  out <- mlm_jn(mod, pred = "x", modx = "m", modx.range = c(-1, 1))
  expect_equal(out$modx.range, c(-1, 1))
  expect_equal(range(out$slopes_df$modx_value), c(-1, 1))
})

test_that("mlm_jn jn_bounds is numeric", {
  mod <- make_test_model()
  out <- mlm_jn(mod, pred = "x", modx = "m")
  expect_true(is.numeric(out$jn_bounds))
})

test_that("mlm_jn sig column is logical", {
  mod <- make_test_model()
  out <- mlm_jn(mod, pred = "x", modx = "m")
  expect_true(is.logical(out$slopes_df$sig))
})

test_that("print.mlm_jn does not error", {
  mod <- make_test_model()
  out <- mlm_jn(mod, pred = "x", modx = "m")
  expect_output(print(out), "Johnson-Neyman")
})

test_that("mlm_jn errors on non-lmer model", {
  expect_error(mlm_jn(list(), pred = "x", modx = "m"), "lmerMod")
})
