context("summary.afmdata")

test_that("afmSummaryTable includes only the fields actually computed", {
  data <- afmReadJPK("force-save-JPK-3h.txt.gz", path = path.package("afmToolkit"))
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
  data <- afmReadJPK("force-save-JPK-2h.txt.gz", path = path.package("afmToolkit"))
  expect_error(capture.output(summary(data, plt = TRUE)), NA)
})

test_that("summary() with plt = TRUE works with a partial analysis", {
  data <- afmReadJPK("force-save-JPK-3h.txt.gz", path = path.package("afmToolkit"))
  data <- afmContactPoint(data, width = 20, mul1 = 1, mul2 = 10)
  data <- afmBaselineCorrection(data)
  expect_error(capture.output(summary(data, plt = TRUE)), NA)
})
