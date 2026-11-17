# Changelog

## precmed 1.1.0.9000 (development version)

### Major changes

- Added support for binary outcomes on the risk-difference scale to
  [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md),
  [`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md),
  and
  [`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md),
  including doubly robust estimation, logistic and boosting scores,
  cross-validated subgroup estimates, and a new vignette. Binary
  cross-validation uses treatment-stratified splits rather than the
  balanced-split controls `error.max` and `max.iter`.

- Standardized bootstrap failure handling across binary, count and
  survival
  [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md):
  affected measures now have unavailable inference if any replicate
  fails or is undefined, instead of silently omitting count replicates.
  Results record replicate-level reasons in `bootstrap.diagnostics` and
  [`diagnostics()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/diagnostics.md).

- Added the population-level `log.cumulative.hazard.ratio` at `tau0` to
  survival
  [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md),
  with bootstrap inference. Inference is unavailable if the estimate or
  any bootstrap replicate is undefined. This replaces the conditional
  Cox `log.hazard.ratio` in
  [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
  output, a breaking change. Subgroup and validation conditional HR
  outputs remain available as `hr.*`; their values can change because of
  the horizon correction described below.

- Added subgroup cumulative hazard ratios to survival
  [`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
  and
  [`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)
  as `chr.*` outputs; `hr.*` outputs remain conditional Cox hazard
  ratios. Survival plots and box plots now accept `measure = "rmtl"`
  (default), `"cumulative_hr"`, or `"conditional_hr"`; validation curves
  can display multiple measures. The previously supported `plot.hr`
  argument is deprecated; `plot.hr = TRUE` now selects the cumulative
  hazard ratio by default rather than the conditional Cox HR. Use
  `measure = "conditional_hr"` to plot the latter. Older saved objects
  must be refitted for cumulative hazard plots.

- Added C backends, enabled by default, for count and survival contrast
  regression, Cox restricted mean survival integration, and LASSO
  propensity-score interaction screening. Installing from source now
  requires a C compiler and R’s standard BLAS/LAPACK libraries.

- Added `contrast.control`, `rmst.control`, and `ps.control` for
  numerical backend settings. Set `solver = "R"` to use the R
  implementation; small floating-point differences between backends are
  possible.

- Replaced `gbm` with `gbm3` for boosting models. Existing tuning
  argument names and the learning rate of 0.1 are preserved; stored fits
  now have class `GBMFit`. Results can differ even with identical
  settings and random seeds.

- Added
  [`simulate_count()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/simulate_count.md),
  [`simulate_binary()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/simulate_binary.md),
  and
  [`simulate_survival()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/simulate_survival.md)
  to generate reproducible synthetic data for illustrating and testing
  the package’s methods with confounding and heterogeneous treatment
  effects.

### Minor changes

- Exported
  [`catefitbinary()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitbinary.md)
  as the binary-outcome counterpart to
  [`catefitcount()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitcount.md)
  and
  [`catefitsurv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitsurv.md),
  with a dedicated help page.

- Expanded the
  [`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
  help to compare CATE score methods and explain which use doubly robust
  estimation and their qualifications.

- Added a package-level overview accessible with
  [`?precmed`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/precmed-package.md),
  including a short example, effect scales, and links to online guides
  and the handbook.

- Added an introductory vignette explaining doubly robust estimation,
  propensity score weighting, outcome-model standardization
  (g-computation), their assumptions, and their use in `precmed`.

- Added
  [`subgroup_estimates()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/subgroup_estimates.md)
  to extract nested binary, count and survival subgroup estimates by
  scoring method and measure. Newly fitted survival
  [`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
  objects also retain subgroup treatment-specific risks and RMSTs and
  risk differences.

- Added [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`confint()`](https://rdrr.io/r/stats/confint.html) methods for
  `catefit` score extraction. Predictions return stored analysis-sample
  scores; coefficient intervals require available contrast-regression
  standard errors and are not available for survival fits.

- Survival
  [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
  bootstrap plots group RMST, risks and the cumulative hazard ratio by
  row by default; `metrics` still controls selection and order and
  accepts underscore spellings of dotted metric names.

- Survival
  [`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)
  now stores subgroup risk differences at `tau0`; validation plots
  support `measure = "risk_difference"`; dotted measure names are
  available as aliases.

- Survival
  [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
  now reports treatment-specific event risks and their difference at
  `tau0`, with bootstrap inference and plot support.

- `diagnostics(fit, type = "propensity")` now reports full-sample
  propensity score overlap, truncation, inverse treatment weights and
  effective sample sizes for newly fitted binary, count and survival
  [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
  results. The default message table is unchanged. Overlap diagnostics
  now default to a compact formatted table; use `layout = "arms"` or
  `layout = "metrics"` for unrounded numeric results; `type = "overlap"`
  remains an alias.

- Binary
  [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
  bootstrap plots now display the treatment-specific risks and risk
  difference by default. `metrics` selects and orders `risk0`, `risk1`,
  and `risk.difference`. Older objects without risk bootstrap estimates
  retain the risk-difference plot.

- Added [`plot()`](https://rdrr.io/r/graphics/plot.default.html) support
  for `catefit` score distributions, with histograms, optional density
  curves, method selection, and a configurable reference line. The score
  axis uses the original treatment labels when available; newly fitted
  `catefit` objects retain these labels, while older fits use treatment
  0/1. Binary, count and survival vignettes now use this method.

- The binary vignette compares logistic and two-regression scores on the
  same internal validation splits, uses 10 iterations and larger
  subgroups to reduce rank-deficient fits, removes the redundant manual
  validation example, and explains the printed average treatment effect
  before showing extraction and plotting commands.

- Standardized documentation terminology for outcome, treatment and
  censoring models, clarifying propensity scores and inverse probability
  of censoring weights.

- Clarified binary, count and survival outcome support in the package
  description, README, quick start and function documentation. The
  README and introductory vignettes also describe survival risk
  differences and distinguish their scale from survival CATE scores. The
  quick start presents the estimated count rate ratio and confidence
  interval in a concise table.

- Added an ATE methods vignette covering binary, count and survival
  outcomes, assumptions and bootstrap inference; the CATE methods
  vignette separates score development from validation.

- Added `combine` to
  [`abc()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/abc.md)
  to aggregate cross-validation iterations per scoring method with a
  mean, median, or user-supplied function. The default still returns the
  full matrix; missing values are not removed automatically.

- Added a dashed null reference line to count and survival
  [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
  bootstrap log-ratio plots. `null.value` defaults to 0 on the log scale
  and can be set to another value, `NULL`, or `NA`; treatment-specific
  RMST panels are unchanged.

- Count
  [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
  now retains bootstrap rate estimates and treatment labels; its
  bootstrap plot displays rates for both treatments alongside the log
  rate ratio. Older saved fits continue to display the log-ratio panel
  only.

- Added `metrics` to
  [`plot.atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/plot.atefit.md)
  to select and order available bootstrap metrics by their names in the
  fitted object.

- Improved bootstrap histogram styling with a minimal default theme and
  panel headings, and added a caption explaining the estimate and null
  reference lines. Supply `theme` to use another ggplot2 theme.
  Histogram bars now default to `alpha = 1` so the background grid does
  not show through them.

- Applied the same default theme to validation curves and subgroup box
  plots; their `theme` argument still accepts a custom ggplot2 theme.
  Captions explain the dashed overall treatment effect reference lines.
  Multi-measure survival plots use horizontal panel headings and matched
  y scales within each measure. Set `show.overall = FALSE` on validation
  curves or box plots to hide the overall effect line. Validation curves
  accept `layout = "method"` to compare training and validation curves
  within method-specific panels. Subgroup box plots now default to two
  columns; `ncol` controls their layout.

- Validate outcome types and outcome-specific CATE scoring methods
  before fitting, including empty and invalid inputs.

- Added `catecvCountExample`, a pre-fitted validation result used in the
  ABC and plotting examples to avoid repeatedly fitting models.

- Label survival
  [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
  bootstrap plot panels with the original treatment names when
  available, and use readable names for the treatment effects.

- Report automatic treatment recoding as an informational message rather
  than a warning, while retaining the mapping in existing result
  objects.

- Added
  [`diagnostics()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/diagnostics.md)
  to summarize recorded survival
  [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
  treatment recoding information and subgroup errors and warnings from
  [`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
  and
  [`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)
  results in a data frame.

- Added [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`confint()`](https://rdrr.io/r/stats/confint.html) methods for
  extracting supported treatment effects and confidence intervals from
  count and survival
  [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
  results.

- Improved performance of balanced splitting, continuous contrast
  regression, censoring probability calculations, and propensity-score
  LASSO fitting.

- Updated validation plots and box plots for current ggplot2 and
  tidyselect versions; ggplot2 \>= 3.4.0 is now required.

- Added the book chapter by Jiang et al. (2026) to the package
  description and references.

- Added Smart Data Analysis and Statistics B.V. as copyright holder for
  its contributions after version 1.1.0, alongside the existing Biogen
  attribution.

### Bug fixes

- Fixed continuous boosting predictions by retaining the selected tree
  count for each treatment arm, and forwarding additional boosting
  controls.

- Fixed survival Cox hazard ratios to censor follow-up at `tau0` instead
  of treating observations beyond the horizon as events. This affects
  subgroup and cross-validation hazard ratio estimates; this correction
  does not change the RMST and RMTL calculations.

- Fixed survival random forest event coding at `tau0` and integrated
  predicted survival curves through the full horizon when computing
  RMST. Survival forest scores and scores using forest initial
  predictions may change.

- Warn when an out-of-range doubly robust RMST estimate is replaced by
  an outcome-model estimate, including when this occurs during bootstrap
  inference.

- Fixed intercept-only propensity models in survival cross-fitting,
  subgroup evaluation, and bootstrap ATE estimation by preserving matrix
  dimensions.

- Fixed LASSO propensity-score dimension handling with two predictors or
  a single prediction row.

- Return `NA` with an informative warning for count-outcome responder
  subgroups missing one or both treatment arms, instead of failing the
  analysis. ABC is `NA` when a required validation estimate is
  unavailable.

- Use GBM3 cross-validation for single-predictor count boosting fits and
  boosting initial predictions, removing the legacy out-of-bag
  workaround and fixing a cross-fitting error when selecting the number
  of trees.

- Removed a duplicate treated-arm boosting fit in count cross-fitting.
  Results obtained with the same seed may change because fewer random
  draws are consumed.

## precmed 1.1.0

CRAN release: 2024-10-05

- Removed dependency on the MESS package.
- Fixed an issue in the cross-validation procedure that lead to an
  increased total sample size
- Updated documentation using roxygen2 (version 7.3.2) and resolved
  bracketing issues.

## precmed 1.0.0.9000

- Added a `NEWS.md` file to track changes to the package.
