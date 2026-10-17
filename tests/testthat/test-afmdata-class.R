# Pure unit tests: independent of any golden file.

test_that("afmdata() builds a valid afmdata object", {
  z <- seq(1e-6, 2e-6, length.out = 10)
  d <- afmdata(data.frame(Z = z, Force = z * 2, Segment = rep("approach", 10)))
  expect_s3_class(d, "afmdata")
  expect_true(is.afmdata(d))
  expect_false(is.afmdata(list(1, 2)))
  expect_false(is.afmdata(1:10))
})

test_that("afmdata() rejects invalid inputs", {
  expect_error(afmdata(data.frame(x = 1)), "valid AFM data frame")
  expect_error(afmdata(list(a = 1)), "does not contain")
  expect_error(afmdata(42), "Not valid data")
})

test_that("afmexperiment() builds a valid afmexperiment", {
  z <- seq(1e-6, 2e-6, length.out = 10)
  df <- data.frame(Z = z, Force = z * 2, Segment = rep("approach", 10))
  e <- afmexperiment(list(afmdata(df), afmdata(df)))
  expect_true(is.afmexperiment(e))
  expect_equal(length(e), 2L)
  expect_error(afmexperiment(list(afmdata(df), 2)), "afmdata class")
})

test_that("append.afmdata appends and replaces fields", {
  z <- seq(1e-6, 2e-6, length.out = 10)
  d <- afmdata(data.frame(Z = z, Force = z * 2, Segment = rep("approach", 10)))
  d2 <- append.afmdata(d, list(a = 1), name = "extra")
  expect_equal(d2$extra, list(a = 1))
  expect_true(is.afmdata(d2))
  expect_setequal(names(d2), c(names(d), "extra"))
  # Appending with the same name replaces the field
  d3 <- append.afmdata(d2, list(b = 2), name = "extra")
  expect_equal(d3$extra, list(b = 2))
  expect_equal(length(d3), length(d2))
})

test_that("trapz integrates correctly", {
  # Trapezoidal rule of y = x^2 over [0, 10] with unit steps
  expect_equal(trapz(0:10, (0:10)^2), 335)
})

test_that("windowedFit matches a direct rolling lm.fit", {
  set.seed(42)
  n <- 200
  width <- 5
  x <- seq(0, 1, length.out = n)
  y <- sin(6 * x) + rnorm(n, sd = 0.05)
  X <- matrix(c(rep(1, n), x, y), nrow = n, ncol = 3)
  b <- as.numeric(windowedFit(X, width))
  expect_equal(length(b), n - 2 * width)
  ref <- vapply((width + 1):(n - width), function(i) {
    idx <- (i - width):(i + width)
    coef(lm.fit(X[idx, 1:2], X[idx, 3]))[2]
  }, numeric(1))
  expect_equal(b, as.numeric(ref))
})
