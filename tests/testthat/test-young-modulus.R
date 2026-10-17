# Golden-master tests: Young's modulus estimation (Hertz contact models).

d2 <- std_pipeline(read_jpk_2h())
d3 <- std_pipeline(read_jpk_3h())

test_that("afmYoungModulus works with a pyramidal tip", {
  d <- afmYoungModulus(d2, thickness = 1e-8, params = list(alpha = 22),
                       silent = TRUE)
  expect_true("YoungModulus" %in% names(d))
  expect_s3_class(d$YoungModulus$fitYM, "lm")
  expect_gte(d$YoungModulus$YoungModulus, 0)
  expect_gt(nrow(d$YoungModulus$fitdata), 3)
  expect_golden("ym-2h-pyramid-th1e-8", list(
    YM = d$YoungModulus$YoungModulus,
    r2 = summary(d$YoungModulus$fitYM)$r.squared,
    nfit = nrow(d$YoungModulus$fitdata),
    fitdata = probe_df(d$YoungModulus$fitdata)
  ))
})

test_that("afmYoungModulus works with a paraboloid tip", {
  d <- afmYoungModulus(d3, thickness = 5e-8, geometry = "paraboloid",
                       params = list(R = 1e-8), silent = TRUE)
  expect_true("YoungModulus" %in% names(d))
  expect_golden("ym-3h-paraboloid-th5e-8", list(
    YM = d$YoungModulus$YoungModulus,
    r2 = summary(d$YoungModulus$fitYM)$r.squared,
    nfit = nrow(d$YoungModulus$fitdata)
  ))
})

test_that("afmYoungModulus works with a conical tip", {
  d <- afmYoungModulus(d2, thickness = 1e-8, geometry = "cone",
                       params = list(alpha = 22), silent = TRUE)
  expect_golden("ym-2h-cone-th1e-8", d$YoungModulus$YoungModulus)
})

test_that("afmYoungModulus requires the previous steps", {
  expect_error(afmYoungModulus(read_jpk_2h(), thickness = 1e-8,
                               params = list(alpha = 22)),
               "Baseline")
})
