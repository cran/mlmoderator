test_that("mlm_center: grand-mean centering works", {
  df <- data.frame(x = c(1, 2, 3, 4, 5), g = c(1, 1, 2, 2, 2))

  out <- mlm_center(df, vars = "x", type = "grand")
  expect_true("x_c" %in% names(out))
  expect_equal(mean(out$x_c), 0, tolerance = 1e-10)
  expect_equal(out$x_c, out$x - mean(out$x))
})

test_that("mlm_center: group-mean centering works", {
  df <- data.frame(x = c(1, 3, 2, 4, 6), g = factor(c(1, 1, 2, 2, 2)))

  out <- mlm_center(df, vars = "x", cluster = "g", type = "group")
  expect_true("x_c" %in% names(out))

  g1_mean <- mean(c(1, 3))
  g2_mean <- mean(c(2, 4, 6))
  expected <- c(1 - g1_mean, 3 - g1_mean, 2 - g2_mean, 4 - g2_mean, 6 - g2_mean)
  expect_equal(out$x_c, expected, tolerance = 1e-10)
})

test_that("mlm_center: both decomposition works", {
  df <- data.frame(x = c(1, 3, 2, 4, 6), g = factor(c(1, 1, 2, 2, 2)))

  out <- mlm_center(df, vars = "x", cluster = "g", type = "both")
  expect_true("x_within" %in% names(out))
  expect_true("x_between" %in% names(out))

  # within + between should reconstruct x
  expect_equal(out$x_within + out$x_between, out$x, tolerance = 1e-10)
})

test_that("mlm_center: errors on missing cluster for group centering", {
  df <- data.frame(x = 1:5, g = factor(1:5))
  expect_error(mlm_center(df, vars = "x", type = "group"), "cluster")
})

test_that("mlm_center: errors on non-existent variable", {
  df <- data.frame(x = 1:5)
  expect_error(mlm_center(df, vars = "z", type = "grand"), "not found")
})

test_that("mlm_center: errors on non-numeric variable", {
  df <- data.frame(x = letters[1:5])
  expect_error(mlm_center(df, vars = "x", type = "grand"), "numeric")
})

test_that("mlm_center: custom suffixes work", {
  df <- data.frame(x = c(1, 2, 3), g = factor(c(1, 1, 2)))
  out <- mlm_center(df, vars = "x", cluster = "g", type = "both",
                    suffix_within = "_w", suffix_between = "_b")
  expect_true("x_w" %in% names(out))
  expect_true("x_b" %in% names(out))
})
