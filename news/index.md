# Changelog

## precmed (development version)

- Added the book chapter by Jiang et al. (2026) to the package
  references and README.
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
