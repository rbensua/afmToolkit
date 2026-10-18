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
