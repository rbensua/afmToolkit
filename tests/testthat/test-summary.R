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

test_that("afmSummaryTable shows the relaxation fit formula and coefficients", {
  data <- std_pipeline(read_jpk_3h())
  data <- afmRelax(data, model = "exp", nexp = 2, type = "CH", plt = FALSE,
                   start = c(a0 = 8.2e-7, a1 = 1.675e-7, a2 = 1.675e-7,
                             tau1 = 5, tau2 = 0.1))
  out <- capture.output(res <- afmSummaryTable(data))
  # The general table comes first...
  expect_equal(names(res$output),
               c("# of segments", "Spring Constant", "Contact Point",
                 "Zero force pt.", "Type of Experiment"))
  expect_equal(res$typeoffit, "Constant Height")
  # The table always has a single row, so no row numbers are printed.
  expect_false(any(grepl("^1\\s", out)))
  # ...and the relaxation fit is printed as its own section, with the
  # fitted formula and the coefficients table.
  expect_true(any(grepl("Relaxation fit \\(model: exp, Constant Height\\):", out)))
  expect_true(any(grepl("ForceCorrected ~ a0 + a1", out)))
  expect_true(any(grepl("Estimate", out)))
  expect_true(any(grepl("tau2", out)))
  # The general table is printed before the relaxation section
  expect_true(which(grepl("Spring Constant", out)) <
              which(grepl("Relaxation fit", out)))
})

test_that("afmSummaryPausePlot draws the two-panel pause figure without warnings", {
  # Regression test: the pause plot used to combine its two panels by
  # hand via gtable internals (g$widths with unit.pmax() and g$layout),
  # which the gtable documentation discourages. The relaxed curve is
  # used so the fit overlay and the coefficients table are exercised.
  data <- std_pipeline(read_jpk_3h())
  data <- afmRelax(data, model = "exp", nexp = 2, type = "CH", plt = FALSE,
                   start = c(a0 = 8.2e-7, a1 = 1.675e-7, a2 = 1.675e-7,
                             tau1 = 5, tau2 = 0.1))
  out <- capture.output(res <- afmSummaryTable(data))
  expect_warning(capture.output(
    afmSummaryPausePlot(data, res$typeoffit)), NA)
})
