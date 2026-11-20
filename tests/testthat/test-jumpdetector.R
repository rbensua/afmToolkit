context("jumpdetector")

test_that("jumpdetector detects a single clean step", {
  x <- c(rep(0, 300), rep(1, 300))
  jd <- jumpdetector(x)
  expect_false(is.null(jd$jump_idx))
  expect_length(jd$jump_idx, 1)
  # The step is between indices 300 and 301.
  expect_true(jd$jump_idx[1] >= 295 && jd$jump_idx[1] <= 306)
  # The amplitude window contains both signal levels.
  expect_equal(jd$Fmin, 0)
  expect_equal(jd$Fmax, 1)
})

test_that("jumpdetector validates its inputs", {
  expect_error(jumpdetector("not a signal"), "numeric")
  expect_error(jumpdetector(c(rep(0, 20), NA)), "missing")
  expect_error(jumpdetector(rep(0, 600), jumpsep = 0), "jumpsep")
  expect_error(jumpdetector(rep(0, 600), jumpdur = 0), "jumpdur")
  expect_error(jumpdetector(rep(0, 10)), "too short")
})

test_that("the figure is drawn without errors", {
  x <- c(rep(0, 300), rep(1, 300))
  expect_error(capture.output(jumpdetector(x, makefigure = TRUE)), NA)
})
