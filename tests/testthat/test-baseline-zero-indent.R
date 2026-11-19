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

# Synthetic curve with an extra modulation segment and piecewise linear
# forces, so that the baseline fits are exact.
mk_modulated <- function() {
  n <- 40
  zApp <- seq(0, 1e-5, length.out = n)
  zMid <- rep(1.05e-5, n)
  zRet <- seq(1.2e-5, 2e-5, length.out = n)
  afmdata(data.frame(
    Z = c(zApp, zMid, zMid, zRet),
    Time = rep(seq_len(n) * 1e-3, times = 4),
    Force = c(3 + 0.5 * zApp,
              rep(10, n),
              20 - 0.1 * zMid,
              -1 + 0.2 * zRet),
    Segment = factor(rep(c("approach", "pause", "modulation", "retract"),
                         each = n))
  ))
}

test_that("afmBaselineCorrection corrects every segment of the curve", {
  d <- mk_modulated()
  dc <- afmBaselineCorrection(d, ZPointApp = 4e-6, ZPointRet = 1.6e-5)
  seg <- d$data$Segment
  fc <- dc$data$ForceCorrected
  # approach and retract are corrected with their own fits: their
  # exactly linear off-contact regions are left at zero
  expect_lt(max(abs(fc[seg == "approach"])), 1e-12)
  expect_lt(max(abs(fc[seg == "retract"])), 1e-12)
  # pause and modulation are corrected with the approach fit (default)
  expect_equal(fc[seg == "pause"],
               d$data$Force[seg == "pause"] - (3 + 0.5 * 1.05e-5))
  expect_equal(fc[seg == "modulation"],
               d$data$Force[seg == "modulation"] - (3 + 0.5 * 1.05e-5))
})

test_that("fitpause = 'retract' corrects the other segments with the retract fit", {
  d <- mk_modulated()
  dc <- afmBaselineCorrection(d, ZPointApp = 4e-6, ZPointRet = 1.6e-5,
                              fitpause = "retract")
  seg <- d$data$Segment
  expect_equal(dc$data$ForceCorrected[seg == "pause"],
               d$data$Force[seg == "pause"] - (-1 + 0.2 * 1.05e-5))
  expect_equal(dc$data$ForceCorrected[seg == "modulation"],
               d$data$Force[seg == "modulation"] - (-1 + 0.2 * 1.05e-5))
})

test_that("fitpause = 'none' keeps the raw force of the other segments", {
  # Regression: the raw forces of the approach and retract segments used
  # to be assigned to the intermediate segments.
  d <- mk_modulated()
  dc <- afmBaselineCorrection(d, ZPointApp = 4e-6, ZPointRet = 1.6e-5,
                              fitpause = "none")
  seg <- d$data$Segment
  expect_equal(dc$data$ForceCorrected[seg == "pause"],
               d$data$Force[seg == "pause"])
  expect_equal(dc$data$ForceCorrected[seg == "modulation"],
               d$data$Force[seg == "modulation"])
  expect_lt(max(abs(dc$data$ForceCorrected[seg == "approach"])), 1e-12)
  expect_lt(max(abs(dc$data$ForceCorrected[seg == "retract"])), 1e-12)
})

test_that("afmBaselineCorrection corrects curves without a retract segment", {
  d <- mk_modulated()
  noRetract <- d$data[d$data$Segment != "retract", ]
  noRetract$Segment <- droplevels(noRetract$Segment)
  dNoRet <- afmdata(noRetract)
  # Previously this failed while trying to fit the absent retract segment
  dc <- afmBaselineCorrection(dNoRet, ZPointApp = 4e-6)
  seg <- dNoRet$data$Segment
  expect_equal(dc$data$ForceCorrected[seg == "pause"],
               dNoRet$data$Force[seg == "pause"] - (3 + 0.5 * 1.05e-5))
  expect_equal(dc$data$ForceCorrected[seg == "modulation"],
               dNoRet$data$Force[seg == "modulation"] - (3 + 0.5 * 1.05e-5))
  expect_error(
    afmBaselineCorrection(dNoRet, ZPointApp = 4e-6, fitpause = "retract"),
    "retract segment"
  )
})

test_that("vsTime corrects every segment with the approach time fit", {
  d <- mk_modulated()
  d$data$Force <- 0.5 + 2 * d$data$Time
  dc <- afmBaselineCorrection(d, ZPointApp = 4e-6, vsTime = TRUE)
  expect_lt(max(abs(dc$data$ForceCorrected)), 1e-12)
})

test_that("vsTime with a sinusoidal fit corrects every segment", {
  d <- mk_modulated()
  w <- 0.2 / (2 * pi)
  # Time must span a good part of the harmonic period for the fit to be
  # identifiable; a sin(Time) term plays the role of non-modelled noise.
  d$data$Time <- rep(seq(0, 150, length.out = 40), times = 4)
  d$data$Force <- with(d$data, -0.5 - 0.001 * Time +
                         0.001 * cos(w * Time + pi / 4) +
                         0.001 * sin(Time))
  dc <- afmBaselineCorrection(d, ZPointApp = 4e-6, vsTime = TRUE,
                              sinusoidal = TRUE)
  expect_true(all(is.finite(dc$data$ForceCorrected)))
  expect_lt(max(abs(dc$data$ForceCorrected)), 0.05)
})
