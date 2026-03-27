## Tests for mlm_probe()
## These tests use a small simulated dataset so they run quickly.

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

test_that("mlm_probe returns an mlm_probe object", {
  mod <- make_test_model()
  out <- mlm_probe(mod, pred = "x", modx = "m")
  expect_s3_class(out, "mlm_probe")
})

test_that("mlm_probe returns 3 rows for mean-sd", {
  mod <- make_test_model()
  out <- mlm_probe(mod, pred = "x", modx = "m")
  expect_equal(nrow(out$slopes), 3)
})

test_that("mlm_probe returns 2 rows for tertiles", {
  mod <- make_test_model()
  out <- mlm_probe(mod, pred = "x", modx = "m", modx.values = "tertiles")
  expect_equal(nrow(out$slopes), 2)
})

test_that("mlm_probe returns 3 rows for quartiles", {
  mod <- make_test_model()
  out <- mlm_probe(mod, pred = "x", modx = "m", modx.values = "quartiles")
  expect_equal(nrow(out$slopes), 3)
})

test_that("mlm_probe works with custom at values", {
  mod <- make_test_model()
  out <- mlm_probe(mod, pred = "x", modx = "m",
                   modx.values = "custom", at = c(-1, 0, 1))
  expect_equal(nrow(out$slopes), 3)
  expect_equal(out$slopes$modx_value, c(-1, 0, 1))
})

test_that("mlm_probe slopes have correct structure", {
  mod <- make_test_model()
  out <- mlm_probe(mod, pred = "x", modx = "m")
  expect_named(out$slopes,
    c("modx_value", "slope", "se", "t", "df", "p", "ci_lower", "ci_upper"))
  expect_true(all(out$slopes$se > 0))
  expect_true(all(out$slopes$p >= 0 & out$slopes$p <= 1))
  expect_true(all(out$slopes$ci_lower < out$slopes$ci_upper))
})

test_that("mlm_probe errors without interaction term", {
  set.seed(1)
  dat <- data.frame(y = rnorm(60), x = rnorm(60), m = rnorm(60),
                    grp = factor(rep(1:10, 6)))
  mod <- suppressMessages(lme4::lmer(y ~ x + m + (1 | grp), data = dat))
  expect_error(mlm_probe(mod, pred = "x", modx = "m"), "interaction")
})

test_that("mlm_probe errors on non-lmer model", {
  expect_error(mlm_probe(list(), pred = "x", modx = "m"), "lmerMod")
})

test_that("print.mlm_probe does not error", {
  mod <- make_test_model()
  out <- mlm_probe(mod, pred = "x", modx = "m")
  expect_output(print(out), "Simple Slopes")
})
