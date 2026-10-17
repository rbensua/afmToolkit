# Golden-master test for the full batch workflow over the 14-curve
# batchExperiment data set (parameters follow the R Journal paper).
# Slowest test of the suite: set AFMTOOLKIT_SKIP_SLOW=true to skip locally.

test_that("full batch pipeline over batchExperiment", {
  skip_slow()
  e <- afmContactPoint(batchExperiment, width = 50, mul1 = 1, mul2 = 10,
                       silent = TRUE)
  e <- afmDetachPoint(e, width = 50, mul1 = 1, mul2 = 10, silent = TRUE)
  e <- afmBaselineCorrection(e)
  e <- afmZeroPointSlope(e, segment = "approach", silent = TRUE)
  e <- afmIndentation(e)
  e <- afmYoungModulus(e, thickness = 2.5e-7, geometry = "paraboloid",
                       params = list(R = 1e-8), silent = TRUE)
  e <- afmRelax(e, model = "power", type = "CH", tmax = 5, plt = FALSE)
  e <- afmAdhesionEnergy(e, width = 5, mul = 15)

  expect_true(is.afmexperiment(e))
  expect_equal(length(e), length(batchExperiment))
  expect_true(all(vapply(e, function(x) {
    all(c("CP", "DP", "YoungModulus", "Relax", "AdhEner") %in% names(x))
  }, logical(1))))

  expect_golden("batch-cp-w50",
                setNames(vapply(e, function(x) x$CP$CP, numeric(1)), names(e)))
  expect_golden("batch-ym-paraboloid",
                vapply(e, function(x) x$YoungModulus$YoungModulus, numeric(1)))
  expect_golden("batch-adhesion",
                lapply(e, function(x) as.numeric(x$AdhEner$Energies[1, ])))

  out <- afmExtract(e, params = list("YM", "AE", "RE"), opt.param = "type")
  expect_type(out, "list")
  expect_s3_class(out$General, "data.frame")
  expect_s3_class(out$Relax, "data.frame")
  expect_equal(nrow(out$General), length(e))
  expect_true("type" %in% names(out$General))
  expect_equal(nrow(out$Relax), 2 * length(e))
  expect_golden("extract-general", probe_df(out$General))
  expect_golden("extract-relax", probe_df(out$Relax))
})

test_that("afmExtract warns for parameters that were not computed", {
  expect_warning(afmExtract(batchExperiment, params = list("YM")),
                 "Young Modulus")
})
