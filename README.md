# afmToolkit

[![CRAN status](https://www.r-pkg.org/badges/version/afmToolkit)](https://CRAN.R-project.org/package=afmToolkit)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/afmToolkit)](https://CRAN.R-project.org/package=afmToolkit)

**afmToolkit** is an R package for the automated analysis of Atomic Force Microscope (AFM)
force-distance and force-time curves, both individually and in batch mode.

It takes a raw curve from the microscope file to the physical quantities:
contact point, adhesion energy, Young's modulus, stress-relaxation fits,
or discrete events such as tethers and protein unfolding ruptures.

## Features

- **Readers**: JPK ASCII files (single curves, whole folders and multi-indent
  maps) and Veeco files, plus batch folder import.
- **Pre-processing**: automatic contact point (`afmContactPoint()`) and detach
  point (`afmDetachPoint()`) detection, and baseline correction
  (`afmBaselineCorrection()`).
- **Mechanical characterization**: indentation (`afmIndentation()`), Young's
  modulus for spherical and conical indenters (`afmYoungModulus()`), adhesion
  energy (`afmAdhesionEnergy()`) and hysteresis (`afmHysteresis()`).
- **Discrete events**: `jumpdetector()` locates jumps (tethers, protein
  unfolding events, ...) in approach/retract signals through the peaks of
  their Haar wavelet transform.
- **Rheology**: `afmRelax()` fits stress-relaxation / creep curves with up to
  two exponentials or a power law, at constant height or constant force.
- **Batch mode**: every function operates transparently on a single curve or
  on a whole experiment (`afmReadJPKFolder()`, `afmReadJPKMultiFolder()`, ...),
  and `afmExtract()` returns one results table for the full experiment.
- **Plots and summaries**: segment-wise plotting (`plot.afmdata`) and a
  `summary()` method that prints the fitted parameters and draws all the
  available segments with their fits.

## Installation

Install the released version from CRAN:

```r
install.packages("afmToolkit")
```

or the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("rbensua/afmToolkit")
```

## Quick start

```r
library(afmToolkit)

## One of the example curves shipped with the package
data <- afmReadJPK("force-save-JPK-3h.txt.gz", path = path.package("afmToolkit"))
data <- afmContactPoint(data, width = 20, mul1 = 1, mul2 = 10)
data <- afmDetachPoint(data, width = 20, mul1 = 2, mul2 = 30)
data <- afmBaselineCorrection(data)
data <- afmZeroPointSlope(data, segment = "approach")
data <- afmIndentation(data)
data <- afmYoungModulus(data, thickness = 1e-7, params = list(alpha = 22),
                        silent = TRUE)
data <- afmRelax(data, model = "exp", nexp = 2, type = "CH")

summary(data, plt = TRUE)
```

The same functions applied to the output of `afmReadJPKFolder()` (or
`afmReadJPKMultiFolder()`) process a whole experiment at once;
`afmExtract()` collects the results into a single table.

## Development (v2)

The `master` branch hosts the v2.0.0 refactoring:

- vectorized rolling regression and adhesion-energy algorithms (O(n)),
- modular internals (`summary.afmdata` split into independent internal modules),
- golden-master regression tests protecting all numeric outputs
  (see `GOLDEN_TESTS.md`), plus input validation and bug fixes throughout.

## Citation

If you use afmToolkit in your research, please cite:

> Benitez, R., Bolos, V. J., & Toca-Herrera, J. L. (2017).
> afmToolkit: an R Package for Automated AFM Force-Distance Curves Analysis.
> *The R Journal*, 9(2), 291-308.
> [doi:10.32614/RJ-2017-045](https://doi.org/10.32614/RJ-2017-045)

## License

GPL.
