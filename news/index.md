# Changelog

## precmed (development version)

- Added optional C integration for Cox restricted mean survival times.
  Set `rmst.control = list(solver = "C")` for survival ATE, fit, or CV
  workflows to avoid allocating patient-by-event-time matrices. C is the
  default; the R backend remains available.

- Added an optional C backend for count contrast regression, selected
  with `contrast.control = list(solver = "C")` in
  [`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)/[`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)
  or their count variants. C is the default; fallback optimizer settings
  can be supplied through `optim.control`. Installing from source now
  requires a C compiler and R’s standard BLAS/LAPACK libraries.

- Compute Breslow censoring probabilities with one cumulative-hazard
  lookup rather than repeatedly scanning and summing hazards for every
  patient.

- Reuse numeric model designs and GLM families during count-data
  balanced splitting, avoiding repeated formula and model-frame
  construction.

- Avoid fitting and discarding a duplicate treated-arm boosting model in
  count cross-fitting. Seeded boosting results can change because fewer
  random draws are consumed.

- Added Smart Data Analysis and Statistics B.V. as copyright holder for
  its contributions after version 1.1.0, alongside the existing Biogen
  attribution.

- Added the book chapter by Jiang et al. (2026) to the package
  description, references, and README.

- Fixed display-math formatting of the weighted estimating equation in
  the theoretical details vignette.

- Updated validation plots to use ggplot2’s `linewidth` argument and
  require ggplot2 \>= 3.4.0.

- Updated column selections in validation plots and box plots for
  current tidyselect versions.

- Declared vignette and test dependencies, removed package installation
  during rendering, and limited boosting examples to one core for
  reliable builds.

- Corrected outdated boosting argument names in the vignettes.

- Invalidated vignette caches when R or package dependencies change.

## precmed 1.1.0

CRAN release: 2024-10-05

- Removed dependency on the MESS package.
- Fixed an issue in the cross-validation procedure that lead to an
  increased total sample size
- Updated documentation using roxygen2 (version 7.3.2) and resolved
  bracketing issues.

## precmed 1.0.0.9000

- Added a `NEWS.md` file to track changes to the package.
