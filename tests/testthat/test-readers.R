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
