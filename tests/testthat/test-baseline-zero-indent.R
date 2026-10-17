# Golden-master tests: baseline correction, zero-force point, indentation.

d2 <- std_pipeline(read_jpk_2h())
d3 <- std_pipeline(read_jpk_3h())

test_that("afmBaselineCorrection adds a ForceCorrected column", {
  expect_true("ForceCorrected" %in% names(d2$data))
  expect_equal(nrow(d2$data), nrow(read_jpk_2h()$data))
  # After correction the off-contact region of the approach must be ~ 0
  app <- d2$data$ForceCorrected[d2$data$Segment == "approach" &
                                  d2$data$Z > d2$CP$CP]
  expect_lt(abs(mean(app)), abs(mean(d2$data$Force[d2$data$Segment ==
                                                    "approach"])))
  expect_golden("baseline-2h", probe(d2$data$ForceCorrected))
  expect_golden("baseline-3h", probe(d3$data$ForceCorrected))
})

test_that("afmBaselineCorrection needs the contact point", {
  expect_error(afmBaselineCorrection(read_jpk_2h()), "ZPointApp")
})

test_that("afmZeroPointSlope finds the zero-force point and the slope", {
  expect_true("Slopes" %in% names(d2))
  expect_true(is.numeric(d2$Slopes$Z0Point))
  expect_true(is.numeric(d2$Slopes$Slope))
  expect_golden("slopes-2h", d2$Slopes)
  expect_golden("slopes-3h", d3$Slopes)
})

test_that("afmZeroPointSlope requires a baseline correction", {
  raw <- afmContactPoint(read_jpk_2h(), width = 20, mul1 = 1, mul2 = 10,
                         silent = TRUE)
  expect_error(afmZeroPointSlope(raw, segment = "approach"), "Baseline")
})

test_that("afmIndentation adds the Indentation column", {
  expect_true("Indentation" %in% names(d2$data))
  expect_equal(nrow(d2$data), nrow(read_jpk_2h()$data))
  expect_golden("indentation-2h", probe(d2$data$Indentation))
})

test_that("afmIndentation requires the previous steps", {
  expect_error(afmIndentation(read_jpk_2h()), "Baseline")
})
