# Golden-master tests: adhesion energies and hysteresis.

d3 <- std_pipeline(read_jpk_3h())

test_that("afmAdhesionEnergy computes the adhesion energies", {
  d <- afmAdhesionEnergy(d3, width = 10, mul = 15)
  expect_true("AdhEner" %in% names(d))
  expect_length(d$AdhEner$Points, 3L)
  expect_equal(nrow(d$AdhEner$Energies), 1L)
  expect_equal(nrow(d$AdhEner$minForce), 1L)
  expect_golden("adhesion-3h-w10-mul15", list(
    points = d$AdhEner$Points,
    energies = as.numeric(d$AdhEner$Energies[1, ]),
    min_force = as.numeric(d$AdhEner$minForce[1, ])
  ))
})

test_that("afmHysteresis quantifies the hysteresis loop", {
  d <- afmHysteresis(d3, silent = TRUE)
  expect_true("Hysteresis" %in% names(d))
  expect_setequal(names(d$Hysteresis$info), c("t1", "t2", "t3", "t4"))
  expect_true(is.numeric(d$Hysteresis$Hysteresis))
  expect_true(is.numeric(d$Hysteresis$Hyst_ratio))
  expect_golden("hysteresis-3h", d$Hysteresis)
})
