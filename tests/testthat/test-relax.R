# Golden-master tests: stress relaxation / creep fits.

d3 <- std_pipeline(read_jpk_3h())

test_that("afmRelax fits a power-law relaxation", {
  d <- afmRelax(d3, model = "power", type = "CH", tmax = 7.5, plt = FALSE)
  expect_true("Relax" %in% names(d))
  expect_equal(d$Relax$model, "power")
  expect_s3_class(d$Relax$relaxModel, "lm")
  expect_equal(length(d$Relax$relaxFit), sum(d3$data$Segment == "pause" &
                                               d3$data$Time <= 7.5))
  expect_golden("relax-3h-power-ch-tmax7.5", list(
    coefs = unname(coef(d$Relax$relaxModel)),
    fit = probe(d$Relax$relaxFit)
  ))
})

test_that("afmRelax fits a two-exponential Prony series", {
  start <- c(a0 = 8.2e-7, a1 = 1.675e-7, a2 = 1.675e-7, tau1 = 5, tau2 = 0.1)
  d <- afmRelax(d3, model = "exp", nexp = 2, type = "CH", tmax = 7.5,
                plt = FALSE, start = start)
  expect_s3_class(d$Relax$relaxModel, "nls")
  expect_golden("relax-3h-exp2-ch-tmax7.5", list(
    coefs = unname(coef(d$Relax$relaxModel)),
    fit = probe(d$Relax$relaxFit)
  ), tolerance = 1e-5)
})

test_that("afmRelax fits a single-exponential decay", {
  start <- c(a0 = 8.2e-7, a1 = 3.35e-7, tau1 = 5)
  d <- afmRelax(d3, model = "exp", nexp = 1, type = "CH", tmax = 7.5,
                plt = FALSE, start = start)
  expect_s3_class(d$Relax$relaxModel, "nls")
  expect_golden("relax-3h-exp1-ch-tmax7.5", list(
    coefs = unname(coef(d$Relax$relaxModel)),
    fit = probe(d$Relax$relaxFit)
  ), tolerance = 1e-5)
})

test_that("afmRelax requires a pause segment", {
  d2 <- std_pipeline(read_jpk_2h())
  expect_error(afmRelax(d2, model = "power", type = "CH", plt = FALSE),
               "pause")
})

test_that("afmRelax handles multi-indentation experiments", {
  skip_slow()
  d <- afmReadJPKMultiIndent("force-save-JPK-multiIndent.txt.gz",
                             path = afm_example_path(), silent = TRUE)
  d <- afmContactPoint(d, width = 100, mul1 = 1, mul2 = 10, silent = TRUE)
  d <- afmDetachPoint(d, width = 100, mul1 = 1, mul2 = 10, silent = TRUE)
  d <- afmBaselineCorrection(d)
  d <- afmZeroPointSlope(d, fstar = 1, silent = TRUE)
  d <- afmIndentation(d)
  d <- afmRelax(d, model = "power", type = "CH", tmax = 5, plt = FALSE)
  expect_true(is.list(d$Relax$relaxModel))
  expect_gte(length(d$Relax$relaxModel), 1)
  expect_golden("relax-multi-power-ch", list(
    model = d$Relax$model,
    coefs = lapply(d$Relax$relaxModel, function(m) unname(coef(m)))
  ))
})
