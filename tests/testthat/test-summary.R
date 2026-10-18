context("summary.afmdata")

test_that("afmSummaryTable includes only the fields actually computed", {
  data <- read_jpk_3h()
  tbl <- capture.output(res <- afmSummaryTable(data))
  expect_equal(names(res$output), c("# of segments", "Spring Constant"))
  expect_null(res$typeoffit)

  data <- afmContactPoint(data, width = 20, mul1 = 1, mul2 = 10)
  tbl <- capture.output(res <- afmSummaryTable(data))
  expect_equal(names(res$output), c("# of segments", "Spring Constant", "Contact Point"))
})

test_that("summary() with plt = TRUE works on a curve with no analysis fields", {
  # Regression test: the parameters table in the approach plot used to
  # break with an "row names supplied are of the wrong length" error when
  # CP, Slopes or YoungModulus were missing.
  data <- read_jpk_2h()
  expect_error(capture.output(summary(data, plt = TRUE)), NA)
})

test_that("summary() with plt = TRUE works with a partial analysis", {
  data <- read_jpk_3h()
  data <- afmContactPoint(data, width = 20, mul1 = 1, mul2 = 10)
  data <- afmBaselineCorrection(data)
  expect_error(capture.output(summary(data, plt = TRUE)), NA)
})
