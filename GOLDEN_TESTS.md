# Golden-master (regression) tests

This package uses **golden-master tests** to guarantee that the numerical
results of the analysis functions do not change while the code is being
refactored (the `v2` rewrite).

## How they work

- Every analysis pipeline is executed on the example data shipped with the
  package (`inst/*.txt.gz`, `data/batchExperiment`).
- `expect_golden()` (defined in `tests/testthat/helper-golden.R`) reduces
  each result to a compact digest with `probe()` / `probe_df()` and stores
  it as an `.rds` file under `tests/testthat/_goldens/`.
- On every subsequent run the freshly computed values are compared against
  the stored goldens: any change in the numerical behaviour fails the tests.

## Workflow

1. **Generate the goldens** (first run, no golden files committed yet):

   ```r
   devtools::test()
   ```

   The first run writes `tests/testthat/_goldens/*.rds` and skips the
   comparisons. **Commit the generated directory** afterwards.

2. **Refactor freely.** Run `devtools::test()` after every change; the tests
   now compare against the committed goldens.

3. **If a numerical change is intentional**, regenerate and commit the goldens
   together with the code change, explaining it in the commit message / PR:

   ```r
   Sys.setenv(AFMTOOLKIT_REGEN_GOLDENS = "true")
   devtools::test()
   ```

The batch test (`test-batch-extract.R`) processes the 14-curve
`batchExperiment` data set and is the slowest one; set
`AFMTOOLKIT_SKIP_SLOW=true` to skip it locally.

`test-paper-values.R` additionally checks (with a loose tolerance) the values
published in the 2017 R Journal article. Those are sanity anchors, not strict
goldens: a failure there deserves investigation, not an automatic update.

## CI

`.github/workflows/R-CMD-check.yaml` runs `R CMD check` (which includes these
tests) on several platforms at every push and pull request.
