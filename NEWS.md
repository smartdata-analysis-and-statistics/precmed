# precmed (development version)

## Major changes

* Added a population-level cumulative hazard ratio at `tau0` to survival
  `atefit()` output. It uses the same outcome, propensity-score and censoring
  models as the RMST analysis, using augmented inverse probability of treatment
  weighting without discretizing event times and with bootstrap inference. The new
  `log.cumulative.hazard.ratio` is included in printing, plotting, `coef()` and
  `confint()`. Inference is unavailable if the estimate or any bootstrap
  replicate is undefined. The conditional Cox `log.hazard.ratio` is no longer
  returned by survival `atefit()`; this is a breaking change for code using
  that output. Subgroup and validation HR outputs are unchanged.

* Added subgroup cumulative hazard ratios to survival `catefit()` and `catecv()`
  as `chr.*` outputs, including training and validation estimates. Existing
  `hr.*` outputs continue to report conditional Cox hazard ratios. For survival
  plots and box plots, select `measure = "rmtl"` (the default),
  `"cumulative_hr"` (population-level cumulative hazard ratio at `tau0`), or
  `"conditional_hr"` (covariate-adjusted Cox hazard ratio fitted using
  follow-up through `tau0`). The Cox ratio summarizes the modeled instantaneous
  hazards over that period, whereas the cumulative hazard ratio compares
  cumulative hazards at `tau0`. Validation curves
  accept multiple measures, for example
  `measure = c("rmtl", "cumulative_hr", "conditional_hr")`, and display them
  in separate panels. Legacy `plot.hr` and `hazard.type` arguments remain
  supported but are deprecated. Older saved objects must be refitted for
  cumulative hazard plots.

* Added `simulate_survival()` to generate reproducible synthetic survival
  data with configurable sample size, propensity-score and outcome-model
  coefficients, and a conditional hazard ratio. The defaults include
  treatment-effect heterogeneity by age and previous relapses, and the
  treatment split approximates the original survival example. The survival
  vignette generates its example with this function; the existing
  `survivalExample` dataset remains available.

* Added C backends, enabled by default, for count and survival contrast
  regression, Cox restricted mean survival integration, and LASSO propensity-score
  interaction screening. Installing from source now requires a C compiler and
  R's standard BLAS/LAPACK libraries.

* Added `contrast.control`, `rmst.control`, and `ps.control` to select numerical
  backends independently of model-fitting arguments. Set `solver = "R"` in the
  relevant list to use the R implementation. `contrast.control` also accepts
  Newton-Raphson settings and nested `optim.control` options for fallback
  optimization. Small floating-point differences between backends are possible.

* Replaced the `gbm` dependency with `gbm3` for boosting models across count,
  survival, and internal continuous-outcome methods. Existing tuning argument
  names and the learning rate of 0.1 are preserved, rather than adopting the
  GBM3 learning-rate default of 0.001. Parallel controls map to GBM3 threads.
  Stored boosting fits now have class `GBMFit`. Results can differ from `gbm`
  even with identical settings and random seeds because random-number usage
  and Poisson updates for zero-event leaves differ between the packages.

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

* Reused model designs and GLM families during count and continuous-outcome
  balanced splitting to avoid repeated formula processing. Continuous-outcome
  public APIs remain disabled as before.

* Used `lm.fit()` for complete-data continuous contrast regression and reused the
  slope generalized inverse in the sandwich covariance, preserving coefficient
  labels and the formula path for incomplete inputs.

* Accelerated Breslow censoring probabilities using cumulative hazards and
  interval lookup instead of repeated scans and sums.

* Reused the full-data fit returned by `cv.glmnet()` instead of fitting the same
  propensity-score LASSO path separately.

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
