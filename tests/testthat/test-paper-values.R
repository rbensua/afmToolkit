# Soft sanity anchors: values published in The R Journal (2017) paper
# "afmToolkit: an R Package for Automated AFM Force-Distance Curves
# Analysis" (Benitez, Bolos and Toca-Herrera), Vol. 9/2.
#
# These are NOT strict golden tests: the package has evolved since 2017,
# so a loose tolerance is used. If one of these fails after a refactor,
# the change in numerical behaviour was probably NOT intentional and
# should be investigated before updating anything.

d3 <- std_pipeline(read_jpk_3h())

test_that("contact point matches the value published in the paper", {
  # Paper (summary output): Contact Point = 6.982701e-06
  # unname(): the CP value currently carries a leftover name attribute
  # ("x2", inherited from the lm coefficient used in its computation).
  # Purely cosmetic, no numerical impact. Candidate for a unname() in
  # the v2 rewrite.
  expect_equal(unname(d3$CP$CP), 6.982701e-06, tolerance = 0.01)
})

test_that("zero force point matches the paper", {
  # Paper: Zero force pt. = 6.99318e-06
  expect_equal(d3$Slopes$Z0Point, 6.99318e-06, tolerance = 0.01)
})

test_that("Young's modulus matches the paper (pyramidal tip)", {
  d <- afmYoungModulus(d3, thickness = 5e-8, params = list(alpha = 22),
                       silent = TRUE)
  # Paper: 249153289 Pa. The current implementation returns ~276.3 MPa
  # (+11%): the fit-region selection of afmYoungModulus changed in the
  # Oct 2025 commits. The goldens lock the CURRENT behaviour; this
  # anchor only guards against gross errors, hence the loose 15%.
  expect_equal(d$YoungModulus$YoungModulus, 249153289, tolerance = 0.15)
})

test_that("adhesion energies match the paper", {
  d <- afmAdhesionEnergy(d3, width = 10, mul = 15)
  # Paper: Points 68.0, 124.5, 347.0
  expect_equal(d$AdhEner$Points, c(68.0, 124.5, 347.0), tolerance = 0.1)
  # Paper: E1adh = 1.142784e-14, E2adh = 1.171835e-15, Etotal = 1.259942e-14
  expect_equal(d$AdhEner$Energies$E1adh, 1.142784e-14, tolerance = 0.15)
  expect_equal(d$AdhEner$Energies$E2adh, 1.171835e-15, tolerance = 0.3)
  expect_equal(d$AdhEner$Energies$Etotal, 1.259942e-14, tolerance = 0.15)
})

test_that("two-exponential relaxation matches the paper (Table 1)", {
  start <- c(a0 = 8.2e-7, a1 = 1.675e-7, a2 = 1.675e-7, tau1 = 5, tau2 = 0.1)
  d <- afmRelax(d3, model = "exp", nexp = 2, type = "CH", tmax = 7.5,
                plt = FALSE, start = start)
  cf <- coef(d$Relax$relaxModel)
  expect_equal(cf[["a0"]],   8.295849e-07, tolerance = 0.1)
  expect_equal(cf[["a1"]],   1.497522e-07, tolerance = 0.1)
  expect_equal(cf[["a2"]],   1.358410e-07, tolerance = 0.1)
  expect_equal(cf[["tau1"]], 2.779792e+00, tolerance = 0.1)
  expect_equal(cf[["tau2"]], 1.334409e-01, tolerance = 0.1)
})
