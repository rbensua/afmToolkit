# afmToolkit 2.0.1 (unreleased)

## Bug fixes

* `afmReadJPK()` can now read single-segment JPK exports that contain
  only the retract curve (e.g. protein-unfolding experiments exported
  from the JPK SPM map viewer), which previously failed with the cryptic
  error "tipo (NULL) inválido para la variable 'Z'" (#13). The channel
  names of the "fancyNames" header line are now matched exactly first
  and then as substrings, so exports shipping channels such as "Height
  (measured & smoothed)" are read correctly, and when several channels
  match the force pattern the deflection channel is preferred. When no
  unique column can be identified, an informative error is thrown. The
  single segment of a one-header file now takes its name from the
  "# segment:" header line instead of always being labelled "approach",
  and a warning is issued when the file carries no spring constant
  (it stays `NA`).
* `afmReadJPK()` failed on JPK 5.0 force-save exports (e.g. CellHesion
  multi-segment experiments) with "could not uniquely identify the
  column matching 'Vertical'. Available columns: .": the "fancy"
  pattern used to locate the column-names line also matched the
  "prefix-fancy"/"suffix-fancy" keys of those exports, and their
  header blocks are written consecutively, without blank lines. The
  column-names line is now located through the "fancyNames" key, the
  data of each segment is delimited by its "# segmentIndex:" key, and
  only the numeric rows of each segment are parsed (#13).
* `afmReadJPK()` silently dropped the last data row of every segment
  except the last one: the end of each non-final segment was cut one
  row short. Every numeric data row of each segment is now kept, so
  files are read with one additional row per non-final segment and
  downstream results (speeds, contact point, slopes, hysteresis, ...)
  change accordingly. The golden master files of the test suite were
  regenerated to lock the corrected behaviour (#13).

## Improvements

* `afmReadJPK()` now reads every segment header of the file instead of
  at most three. The name of each segment is taken from its
  "# segment:" key ("extend" is renamed to "approach"; "pause",
  "retract", "modulation" or any other declared name is kept as is),
  so experiments with modulation segments or with several approach or
  pause segments are read completely. Header blocks without the key
  keep the positional names approach/pause/retract (#13).

# afmToolkit 2.0.0

## Major internal refactoring, no user-facing API changes

This is a full internal rewrite of the package, carried out on top of a
new golden-master regression test suite so that every intermediate step
was verified to produce identical (or explicitly documented) outputs.

### Testing infrastructure

* New testthat test suite (edition 3) with golden-master regression
  tests covering the readers, the whole analysis pipeline and the values
  published in the R Journal 2017 paper (#1).
* Continuous integration via a GitHub Actions workflow running
  R CMD check on Linux, macOS and Windows (#1).
* `GOLDEN_TESTS.md` documents the test workflow, including how to
  regenerate the goldens (`AFMTOOLKIT_REGEN_GOLDENS=true`).

### Performance

* `windowedFit()` now computes every rolling least-squares slope from
  prefix sums in O(n), instead of running one `lm.fit()` per window
  (O(n * w)). This speeds up `afmContactPoint()`, `afmDetachPoint()`,
  `afmAdhesionEnergy()` and the batch processing of experiments (#2).
* `afmAdhesionEnergy()` computes the cumulative adhesion energy via
  `cumsum()` in O(n), instead of an O(n^2) loop of repeated trapezoidal
  integrations (#2).

### Bug fixes

* `afmBaselineCorrection()` applied to an `afmexperiment` silently
  ignored the `fitpause`, `vsTime` and `sinusoidal` arguments: every
  curve was corrected with the default options. They are now honoured
  at experiment level, exactly as for single curves (#5).
* `afmHysteresis()` applied to an `afmexperiment` silently dropped the
  extra `...` arguments, so every curve was processed with the default
  options of `afmZeroPointSlope()`. They are now forwarded to every
  curve (#5).
* The contact point (CP) and detach point (DP) values returned by
  `afmContactPoint()` and `afmDetachPoint()` no longer carry a leftover
  `"x2"` name attribute, an artifact inherited from the `lm.fit()`
  coefficient names that propagated through the delta signal (#5).

### Internal improvements

* New internal `afmDispatch()` helper centralises the
  `afmdata`/ `afmexperiment` dispatch, the "Processing curve" message
  and the invalid-input error; all nine analysis functions are now thin
  wrappers over internal single-curve functions (#4).
* `afmReadJPK()` no longer repeats the same segment-parsing block three
  times (approach/pause/retract); it uses two internal helpers (#3).
* Dead commented-out code and unused locals removed throughout (#6).

### Note for users comparing with the published paper

The Young's modulus computed for the example data with a pyramidal tip
(thickness = 5e-8 m, alpha = 22 degrees) currently gives ~276 MPa,
versus ~249 MPa published in the R Journal 2017 article. This drift was
introduced by changes made after the paper's publication and is locked
and documented by the new test suite (with a 15% tolerance).
