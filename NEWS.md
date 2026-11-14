# precmed (development version)

## Major changes

* Updated the default `simulate_survival()` coefficients to illustrate
  confounding and heterogeneous treatment effects. Simulated results for a
  given seed will differ from earlier development versions.

* Added the population-level `log.cumulative.hazard.ratio` at `tau0` to
  survival `atefit()`, with bootstrap inference. Inference is unavailable if
  the estimate or any bootstrap replicate is undefined. This replaces the
  conditional Cox `log.hazard.ratio` in `atefit()` output, a breaking change.
  Subgroup and validation HR outputs are unchanged.

* Added subgroup cumulative hazard ratios to survival `catefit()` and `catecv()`
  as `chr.*` outputs; `hr.*` outputs remain conditional Cox hazard ratios.
  Survival plots and box plots now accept `measure = "rmtl"` (default),
  `"cumulative_hr"`, or `"conditional_hr"`; validation curves can display
  multiple measures. `plot.hr` and `hazard.type` are deprecated. Older saved
  objects must be refitted for cumulative hazard plots.

* Added `simulate_survival()` to generate reproducible synthetic survival data.
  The survival vignette now uses it; `survivalExample` remains available.

* Added C backends, enabled by default, for count and survival contrast
  regression, Cox restricted mean survival integration, and LASSO propensity-score
  interaction screening. Installing from source now requires a C compiler and
  R's standard BLAS/LAPACK libraries.

* Added `contrast.control`, `rmst.control`, and `ps.control` for numerical
  backend settings. Set `solver = "R"` to use the R implementation; small
  floating-point differences between backends are possible.

* Replaced `gbm` with `gbm3` for boosting models. Existing tuning argument
  names and the learning rate of 0.1 are preserved; stored fits now have class
  `GBMFit`. Results can differ even with identical settings and random seeds.

* Fixed continuous boosting predictions by retaining the selected tree count
  for each treatment arm, and forwarding additional boosting controls.

## Minor changes and bug fixes

* Fixed survival Cox hazard ratios to censor follow-up at `tau0` instead of
  treating observations beyond the horizon as events. This affects
  subgroup and cross-validation hazard ratio estimates; RMST and RMTL
  calculations are unchanged.

* Fixed survival random forest event coding at `tau0` and integrated predicted
  survival curves through the full horizon when computing RMST. Survival forest
  scores and scores using forest initial predictions may change.

* Warn when an out-of-range doubly robust RMST estimate is replaced by an
  outcome-model estimate, including when this occurs during bootstrap inference.

* Fixed intercept-only propensity models in survival cross-fitting, subgroup
  evaluation, and bootstrap ATE estimation by preserving matrix dimensions.

* Fixed LASSO propensity-score dimension handling with two predictors or a single
  prediction row.

* Return `NA` with an informative warning for count-outcome responder subgroups
  missing one or both treatment arms, instead of failing the analysis.
  ABC is `NA` when a required validation estimate is unavailable.

* Use GBM3 cross-validation for single-predictor count boosting fits and
  boosting initial predictions, removing the legacy out-of-bag workaround and
  fixing a cross-fitting error when selecting the number of trees.

* Removed a duplicate treated-arm boosting fit in count cross-fitting. Results
  obtained with the same seed may change because fewer random draws are consumed.

* Improved performance of balanced splitting, continuous contrast regression,
  censoring probability calculations, and propensity-score LASSO fitting.

* Updated validation plots and box plots for current ggplot2 and tidyselect
  versions; ggplot2 >= 3.4.0 is now required.

* Added the book chapter by Jiang et al. (2026) to the package description and
  references.

* Added Smart Data Analysis and Statistics B.V. as copyright holder for its
  contributions after version 1.1.0, alongside the existing Biogen attribution.

# precmed 1.1.0

* Removed dependency on the MESS package.
* Fixed an issue in the cross-validation procedure that lead to an increased total sample size
* Updated documentation using roxygen2 (version 7.3.2) and resolved bracketing issues.

# precmed 1.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
