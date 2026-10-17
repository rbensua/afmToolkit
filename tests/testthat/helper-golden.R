# Golden-master (regression) test infrastructure --------------------------
#
# Workflow (see GOLDEN_TESTS.md at the repo root):
#
# 1. First run: every expect_golden() stores the value it receives as an
#    .rds file under tests/testthat/_goldens/ and the comparison is
#    skipped. Commit the generated _goldens/ directory afterwards.
# 2. Every later run compares against the stored goldens: any change in
#    the numerical behaviour of the package fails the tests.
# 3. If a numerical change is INTENTIONAL, regenerate the goldens,
#    review the .rds diff and commit it together with the change:
#      Sys.setenv(AFMTOOLKIT_REGEN_GOLDENS = "true"); devtools::test()
#
# Long vectors / data frames are reduced with probe()/probe_df() so the
# golden files stay small and readable.

afm_golden_dir <- function() {
  testthat::test_path("_goldens")
}

afm_golden_path <- function(name) {
  file.path(afm_golden_dir(), paste0(name, ".rds"))
}

afm_golden_regen <- function() {
  identical(Sys.getenv("AFMTOOLKIT_REGEN_GOLDENS"), "true")
}

expect_golden <- function(name, value, tolerance = 1e-6) {
  force(value)
  path <- afm_golden_path(name)
  if (afm_golden_regen() || !file.exists(path)) {
    dir.create(afm_golden_dir(), recursive = TRUE, showWarnings = FALSE)
    saveRDS(value, path)
    testthat::skip(sprintf(
      "Golden '%s' %s; run the tests again to compare.",
      name, if (afm_golden_regen()) "regenerated" else "created"
    ))
  }
  expected <- readRDS(path)
  testthat::expect_equal(
    value, expected,
    tolerance = tolerance,
    info = sprintf("golden '%s'", name)
  )
}

# Compact but highly sensitive digest of a long numeric vector
probe <- function(x, n = 25) {
  if (!is.numeric(x)) return(head(x, n))
  if (length(x) == 0L) return(c(n = 0))
  idx <- unique(round(seq.int(1L, length(x), length.out = min(n, length(x)))))
  c(
    n      = length(x),
    sum    = sum(x),
    mean   = mean(x),
    sd     = stats::sd(x),
    min    = min(x),
    q1     = unname(stats::quantile(x, 0.25)),
    median = stats::median(x),
    q3     = unname(stats::quantile(x, 0.75)),
    max    = max(x),
    strided = unname(x[idx])
  )
}

# Column-wise probe of a data frame (numeric columns are probed, the
# rest are just headed)
probe_df <- function(df, cols = names(df), n = 25) {
  stopifnot(is.data.frame(df))
  out <- lapply(df[cols], function(col) {
    if (is.numeric(col)) probe(col, n = n) else head(col, n)
  })
  names(out) <- cols
  out
}

# Shared fixtures -----------------------------------------------------------

afm_example_path <- function(...) {
  system.file(..., package = "afmToolkit")
}

read_jpk_2h <- function() {
  afmReadJPK("force-save-JPK-2h.txt.gz", path = afm_example_path(), silent = TRUE)
}

read_jpk_3h <- function() {
  afmReadJPK("force-save-JPK-3h.txt.gz", path = afm_example_path(), silent = TRUE)
}

# Standard preprocessing pipeline used by many tests
std_pipeline <- function(afmdata, width = 20, mul1 = 1, mul2 = 10) {
  afmdata <- afmContactPoint(afmdata, width = width, mul1 = mul1, mul2 = mul2,
                             silent = TRUE)
  afmdata <- afmDetachPoint(afmdata, width = width, mul1 = mul1, mul2 = mul2,
                            silent = TRUE)
  afmdata <- afmBaselineCorrection(afmdata)
  afmdata <- afmZeroPointSlope(afmdata, segment = "approach", silent = TRUE)
  afmIndentation(afmdata)
}

skip_slow <- function() {
  if (identical(Sys.getenv("AFMTOOLKIT_SKIP_SLOW"), "true")) {
    testthat::skip("Slow test (set AFMTOOLKIT_SKIP_SLOW=true to skip locally)")
  }
}
