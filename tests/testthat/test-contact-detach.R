# Golden-master tests: contact point and detach point detection
# (Benitez et al. 2013 algorithm).

d2 <- read_jpk_2h()
d3 <- read_jpk_3h()

test_that("afmContactPoint finds the contact point (2-segment curve)", {
  d <- afmContactPoint(d2, width = 20, mul1 = 1, mul2 = 10, silent = TRUE)
  expect_true("CP" %in% names(d))
  n_approach <- sum(d2$data$Segment == "approach")
  expect_equal(length(d$CP$delta), n_approach)
  expect_gte(d$CP$CP, min(d2$data$Z))
  expect_lte(d$CP$CP, max(d2$data$Z))
  expect_golden("cp-2h-w20", list(
    CP = d$CP$CP,
    iCP = d$CP$iCP,
    noise = d$CP$noise,
    delta = probe(d$CP$delta)
  ))
})

test_that("afmDetachPoint finds the detach point (2-segment curve)", {
  d <- afmDetachPoint(d2, width = 20, mul1 = 1, mul2 = 10, silent = TRUE)
  expect_true("DP" %in% names(d))
  expect_golden("dp-2h-w20", list(
    DP = d$DP$DP,
    iDP = d$DP$iDP,
    noise = d$DP$noise,
    delta = probe(d$DP$delta)
  ))
})

test_that("contact and detach points on the 3-segment curve", {
  cp <- afmContactPoint(d3, width = 20, mul1 = 1, mul2 = 10, silent = TRUE)
  dp <- afmDetachPoint(d3, width = 20, mul1 = 1, mul2 = 10, silent = TRUE)
  expect_golden("cp-3h-w20",
                list(CP = cp$CP$CP, iCP = cp$CP$iCP, noise = cp$CP$noise))
  expect_golden("dp-3h-w20",
                list(DP = dp$DP$DP, iDP = dp$DP$iDP, noise = dp$DP$noise))
})

test_that("afmContactPoint rejects non-afmdata input", {
  expect_error(afmContactPoint(1:10, width = 5, mul1 = 1, mul2 = 5), "afmdata")
  expect_error(afmDetachPoint(1:10, width = 5, mul1 = 1, mul2 = 5), "afmdata")
})

test_that("afmContactPoint dispatches over an afmexperiment", {
  e <- afmReadJPKFolder(afm_example_path("afmexperiment"), silent = TRUE)
  e <- afmContactPoint(e, width = 20, mul1 = 1, mul2 = 10, silent = TRUE)
  expect_true(is.afmexperiment(e))
  expect_true(all(vapply(e, function(x) "CP" %in% names(x), logical(1))))
  expect_golden("cp-folder-w20",
                setNames(vapply(e, function(x) x$CP$CP, numeric(1)), names(e)))
})
