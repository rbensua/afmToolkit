# Golden-master tests for the file readers.

test_that("afmReadJPK reads a 2-header JPK file", {
  d <- read_jpk_2h()
  expect_true(is.afmdata(d))
  expect_setequal(levels(d$data$Segment), c("approach", "retract"))
  expect_setequal(names(d$data), c("Z", "Force", "Time", "Segment"))
  expect_gt(d$params$SpringConstant, 0)
  expect_equal(d$params$curvename, "force-save-JPK-2h.txt.gz")
  expect_golden("reader-jpk2h", list(
    nrow = nrow(d$data),
    rows_per_segment = as.list(table(d$data$Segment)),
    data = probe_df(d$data, c("Z", "Force", "Time")),
    spring = d$params$SpringConstant,
    speeds = d$params$speeds,
    index = d$params$Index,
    xpos = d$params$xpos,
    ypos = d$params$ypos
  ))
})

test_that("afmReadJPK reads a 3-header JPK file", {
  d <- read_jpk_3h()
  expect_true(is.afmdata(d))
  expect_setequal(levels(d$data$Segment), c("approach", "pause", "retract"))
  expect_golden("reader-jpk3h", list(
    nrow = nrow(d$data),
    rows_per_segment = as.list(table(d$data$Segment)),
    data = probe_df(d$data, c("Z", "Force", "Time")),
    spring = d$params$SpringConstant,
    speeds = d$params$speeds,
    index = d$params$Index,
    xpos = d$params$xpos,
    ypos = d$params$ypos
  ))
})

test_that("afmReadJPK keeps every numeric data row of each segment", {
  # Regression test: the end of each non-final segment used to be cut
  # one row short, silently dropping the last data row of the segment.
  for (file in c("force-save-JPK-2h.txt.gz", "force-save-JPK-3h.txt.gz")) {
    lines <- readLines(gzfile(afm_example_path(file)))
    n_numeric <- sum(!grepl("^#", lines) & trimws(lines) != "")
    d <- afmReadJPK(file, path = afm_example_path(), silent = TRUE)
    expect_equal(nrow(d$data), n_numeric, info = file)
  }
})

test_that("afmReadVeeco reads a Veeco ascii file", {
  d <- afmReadVeeco("veeco_file.txt.gz", path = afm_example_path(), silent = TRUE)
  expect_true(is.afmdata(d))
  expect_setequal(levels(d$data$Segment), c("approach", "retract"))
  expect_golden("reader-veeco", list(
    nrow = nrow(d$data),
    rows_per_segment = as.list(table(d$data$Segment)),
    data = probe_df(d$data, c("Z", "Force", "Time")),
    spring = d$params$SpringConstant
  ))
})

test_that("afmReadJPKFolder reads all the JPK files of a folder", {
  e <- afmReadJPKFolder(afm_example_path("afmexperiment"), silent = TRUE)
  expect_true(is.afmexperiment(e))
  expect_equal(length(e), 4L)
  expect_true(all(vapply(e, is.afmdata, logical(1))))
  expect_golden("reader-jpkfolder", list(
    names = names(e),
    springs = vapply(e, function(x) x$params$SpringConstant, numeric(1)),
    nrows = vapply(e, function(x) nrow(x$data), integer(1))
  ))
})

test_that("afmReadVeecoFolder reads all the Veeco files of a folder", {
  e <- afmReadVeecoFolder(afm_example_path("veecoFolder"), silent = TRUE)
  expect_true(is.afmexperiment(e))
  expect_equal(length(e), 4L)
  expect_golden("reader-veecofolder", list(
    names = names(e),
    springs = vapply(e, function(x) x$params$SpringConstant, numeric(1)),
    nrows = vapply(e, function(x) nrow(x$data), integer(1))
  ))
})

test_that("afmReadJPKMultiIndent reads a multi-indentation experiment", {
  d <- afmReadJPKMultiIndent("force-save-JPK-multiIndent.txt.gz",
                             path = afm_example_path(), silent = TRUE)
  expect_true(is.afmdata(d))
  expect_true(is.afmmulti(d))
  segs <- levels(d$data$Segment)
  expect_true(all(grepl("^(approach|pause)[0-9]+$|^retract$", segs)))
  expect_golden("reader-jpkmulti", list(
    nrow = nrow(d$data),
    segments = segs,
    rows_per_segment = as.list(table(d$data$Segment)),
    spring = d$params$SpringConstant
  ))
})

test_that("afmReadJPK reads a single-segment retract-only JPK file", {
  d <- afmReadJPK(test_path("jpk-retract-only.txt"), silent = TRUE)
  expect_true(is.afmdata(d))
  expect_setequal(names(d$data), c("Z", "Force", "Time", "Segment"))
  expect_false("SegmentID" %in% names(d$data))
  expect_setequal(levels(d$data$Segment), "retract")
  expect_equal(nrow(d$data), 100L)
  # Z comes from "Height (measured & smoothed)" and Force from
  # "Vertical Deflection", not from "Vertical Tip Position"
  expect_equal(d$data$Z[1], -2.4674994E-8)
  expect_equal(d$data$Force[1], 1.9822173E-9)
  expect_equal(d$data$Time[1], 4.0E-4)
  expect_gt(d$params$speeds$retract, 0)
  # This export carries neither a spring constant nor index/position keys
  expect_true(is.na(d$params$SpringConstant))
  expect_true(is.na(d$params$Index))
  expect_true(is.na(d$params$xpos))
  expect_true(is.na(d$params$ypos))
})

test_that("afmReadJPK warns when the file has no spring constant", {
  expect_warning(
    capture.output(afmReadJPK(test_path("jpk-retract-only.txt"), silent = FALSE)),
    "No spring constant"
  )
})

test_that("afmReadJPK keeps 'approach' when a single-segment file declares no segment", {
  lines <- c(
    "# segmentIndex: 1",
    "# columns: height vDeflection time",
    "# fancyNames: \"Height (measured)\" \"Vertical Deflection\" \"Segment Time\"",
    "# units: m N s",
    "#",
    "0 0 0",
    "1e-9 1e-9 0.1",
    "2e-9 2e-9 0.2"
  )
  f <- tempfile(fileext = ".txt")
  writeLines(lines, f)
  d <- afmReadJPK(f, silent = TRUE)
  expect_setequal(levels(d$data$Segment), "approach")
  expect_equal(d$data$Z, c(0, 1e-9, 2e-9))
  unlink(f)
})

test_that("afmReadJPK fails informatively when no column matches the pattern", {
  lines <- c(
    "# spring constant: 0.5",
    "# segment: retract",
    "# columns: zSensor vDeflection time",
    "# fancyNames: \"Z Sensor\" \"Vertical Deflection\" \"Segment Time\"",
    "# units: m N s",
    "#",
    "0 0 0",
    "1e-9 1e-9 0.1",
    "2e-9 2e-9 0.2"
  )
  f <- tempfile(fileext = ".txt")
  writeLines(lines, f)
  expect_error(afmReadJPK(f, silent = TRUE), "could not uniquely identify")
  unlink(f)
})

test_that("afmReadJPK maps 'extend' to 'approach' and keeps other declared names", {
  lines <- c(
    "# segmentIndex: 1",
    "# segment: extend",
    "# columns: height vDeflection time",
    "# fancyNames: \"Height (measured)\" \"Vertical Deflection\" \"Segment Time\"",
    "# units: m N s",
    "#",
    "0 0 0",
    "1e-9 1e-9 0.1",
    "2e-9 2e-9 0.2",
    " ",
    "# segmentIndex: 2",
    "# segment: retract",
    "# columns: height vDeflection time",
    "# fancyNames: \"Height (measured)\" \"Vertical Deflection\" \"Segment Time\"",
    "# units: m N s",
    "#",
    "3e-9 3e-9 0.3",
    "4e-9 4e-9 0.4",
    "5e-9 5e-9 0.5"
  )
  f <- tempfile(fileext = ".txt")
  writeLines(lines, f)
  d <- afmReadJPK(f, silent = TRUE)
  expect_setequal(levels(d$data$Segment), c("approach", "retract"))
  expect_equal(as.integer(table(d$data$Segment)), c(3L, 3L))
  expect_equal(d$data$Z[d$data$Segment == "approach"], c(0, 1e-9, 2e-9))
  expect_equal(d$data$Z[d$data$Segment == "retract"], c(3e-9, 4e-9, 5e-9))
  expect_equal(sort(names(d$params$speeds)), c("approach", "retract"))
  unlink(f)
})

test_that("afmReadJPK reads every header, including modulation and repeated segments", {
  lines <- c(
    "# segmentIndex: 1",
    "# segment: extend",
    "# columns: height vDeflection time",
    "# fancyNames: \"Height (measured)\" \"Vertical Deflection\" \"Segment Time\"",
    "# units: m N s",
    "#",
    "0 0 0",
    "1e-9 1e-9 0.1",
    " ",
    "# segmentIndex: 2",
    "# segment: modulation",
    "# columns: height vDeflection time",
    "# fancyNames: \"Height (measured)\" \"Vertical Deflection\" \"Segment Time\"",
    "# units: m N s",
    "#",
    "2e-9 2e-9 0.2",
    "3e-9 3e-9 0.3",
    " ",
    "# segmentIndex: 3",
    "# segment: extend",
    "# columns: height vDeflection time",
    "# fancyNames: \"Height (measured)\" \"Vertical Deflection\" \"Segment Time\"",
    "# units: m N s",
    "#",
    "4e-9 4e-9 0.4",
    "5e-9 5e-9 0.5",
    " ",
    "# segmentIndex: 4",
    "# segment: retract",
    "# columns: height vDeflection time",
    "# fancyNames: \"Height (measured)\" \"Vertical Deflection\" \"Segment Time\"",
    "# units: m N s",
    "#",
    "6e-9 6e-9 0.6",
    "7e-9 7e-9 0.7"
  )
  f <- tempfile(fileext = ".txt")
  writeLines(lines, f)
  d <- afmReadJPK(f, silent = TRUE)
  expect_setequal(levels(d$data$Segment), c("approach", "modulation", "retract"))
  expect_equal(as.integer(table(d$data$Segment)), c(4L, 2L, 2L))
  expect_equal(nrow(d$data), 8L)
  expect_equal(sort(names(d$params$speeds)), c("approach", "modulation", "retract"))
  # Two "extend" (approach) segments: the SegmentID column disambiguates
  # them; unique names keep their plain segment name
  expect_true("SegmentID" %in% names(d$data))
  expect_equal(as.character(d$data$SegmentID),
               c("approach_1", "approach_1", "modulation", "modulation",
                 "approach_2", "approach_2", "retract", "retract"))
  unlink(f)
})

test_that("afmReadJPK falls back to positional names when no segment key is declared", {
  lines <- c(
    "# segmentIndex: 1",
    "# columns: height vDeflection time",
    "# fancyNames: \"Height (measured)\" \"Vertical Deflection\" \"Segment Time\"",
    "# units: m N s",
    "#",
    "0 0 0",
    "1e-9 1e-9 0.1",
    " ",
    "# segmentIndex: 2",
    "# columns: height vDeflection time",
    "# fancyNames: \"Height (measured)\" \"Vertical Deflection\" \"Segment Time\"",
    "# units: m N s",
    "#",
    "2e-9 2e-9 0.2",
    "1e-9 1e-9 0.3"
  )
  f <- tempfile(fileext = ".txt")
  writeLines(lines, f)
  d <- afmReadJPK(f, silent = TRUE)
  expect_setequal(levels(d$data$Segment), c("approach", "retract"))
  expect_equal(as.integer(table(d$data$Segment)), c(2L, 2L))
  # The retract segment is normalised to increasing Z even without the key
  expect_equal(d$data$Z[d$data$Segment == "retract"], c(1e-9, 2e-9))
  unlink(f)
})

test_that("afmReadJPK reads JPK 5.0 force-save exports (fancyNames keys, consecutive headers)", {
  lines <- c(
    "# TEXT EXPORT",
    "# index: 0",
    "# xPosition: 1.0E-13",
    "# yPosition: -2.0E-13",
    "# force-settings.segment.0.identifier.prefix-fancy: ",
    "# force-settings.segment.0.identifier.suffix-fancy:  (1)",
    "#",
    "# segmentIndex: 0",
    "# segment: extend",
    "# columns: smoothedMeasuredHeight vDeflection measuredHeight seriesTime time",
    "# fancyNames: \"Height (measured & smoothed)\" \"Vertical Deflection\" \"Height (measured)\" \"Series Time\" \"Segment Time\"",
    "# springConstant: 0.12",
    "# units: m N m s s",
    "#",
    "1e-8 0 0 0 0",
    "1.1e-8 1e-9 1e-9 0.1 0.1",
    "1.2e-8 2e-9 2e-9 0.2 0.2",
    "#",
    "# segmentIndex: 1",
    "# segment: pause",
    "# columns: smoothedMeasuredHeight vDeflection measuredHeight seriesTime time",
    "# fancyNames: \"Height (measured & smoothed)\" \"Vertical Deflection\" \"Height (measured)\" \"Series Time\" \"Segment Time\"",
    "# units: m N m s s",
    "#",
    "1.3e-8 3e-9 3e-9 0.3 0.3",
    "1.4e-8 4e-9 4e-9 0.4 0.4",
    "#",
    "# segmentIndex: 2",
    "# segment: retract",
    "# columns: smoothedMeasuredHeight vDeflection measuredHeight seriesTime time",
    "# fancyNames: \"Height (measured & smoothed)\" \"Vertical Deflection\" \"Height (measured)\" \"Series Time\" \"Segment Time\"",
    "# units: m N m s s",
    "#",
    "1.6e-8 6e-9 6e-9 0.6 0.6",
    "1.5e-8 5e-9 5e-9 0.7 0.7"
  )
  f <- tempfile(fileext = ".txt")
  writeLines(lines, f)
  d <- afmReadJPK(f, silent = TRUE)
  expect_true(is.afmdata(d))
  expect_setequal(levels(d$data$Segment), c("approach", "pause", "retract"))
  expect_equal(as.integer(table(d$data$Segment)), c(3L, 2L, 2L))
  expect_false("SegmentID" %in% names(d$data))
  # Z comes from the exact "Height (measured)" channel, not from
  # "Height (measured & smoothed)", despite the "prefix-fancy"/
  # "suffix-fancy" keys and the consecutively written header blocks
  expect_equal(d$data$Z, c(0, 1e-9, 2e-9, 3e-9, 4e-9, 5e-9, 6e-9))
  expect_equal(d$data$Force, c(0, 1e-9, 2e-9, 3e-9, 4e-9, 5e-9, 6e-9))
  expect_equal(sort(names(d$params$speeds)), c("approach", "pause", "retract"))
  expect_equal(d$params$SpringConstant, 0.12)
  expect_equal(d$params$Index, 0)
  unlink(f)
})
