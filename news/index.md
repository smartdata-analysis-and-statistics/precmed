# Changelog

## precmed 1.1.0.9000 (development version)

### Major changes

- Survival T-learners in
  [`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
  and
  [`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)
  support `outcome_method = "coxph"`, `"weibull"`, and `"log_normal"` to
  compare arm-specific restricted mean time lost at the selected
  horizon.

- [`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)
  and its outcome-specific functions use the same model and method
  arguments as
  [`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md).
  Defaults are logistic, Poisson, and Cox T-learners for binary, count,
  and survival outcomes. Two-regression and contrast regression can use
  a separate `benefit_model` and different outcome methods within one
  call. Deprecated argument names and positional calls remain supported
  with warnings.
  [`catecvbinary()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecvbinary.md)
  is now exported.

- Count
  [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
  now always fits separate outcome models by treatment arm, without the
  previous fallback to a pooled model for small arms. The `interactions`
  argument has been removed from ATE fitting, including the option for a
  pooled model without interactions. This is a breaking change: calls
  using `interactions` must be updated, and positional arguments that
  followed it have shifted. Results and bootstrap inference can change.

- Added support for binary outcomes on the risk-difference scale to
  [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md),
  [`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md),
  and
  [`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md),
  including doubly robust estimation, logistic and boosting scores,
  cross-validated subgroup estimates, and a new vignette. Binary
  [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
  fits separate logistic outcome models by treatment arm. Binary
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

- [`catefitcount()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitcount.md)
  no longer measures or prints runtime when verbose output is enabled.

- `initial.predictor.method` is no longer shown in
  [`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
  or its outcome-specific signatures. Use `outcome_method` to select
  initial outcome models; existing calls using the dotted argument
  remain supported through `...` with a deprecation warning. Supplying
  it together with `outcome_method` now produces an error, including
  when the values match. Legacy calls can pair the dotted argument with
  `benefit_method`.

- [`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md),
  [`catefitbinary()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitbinary.md),
  [`catefitcount()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitcount.md),
  and
  [`catefitsurv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitsurv.md)
  now show `benefit_method = "t_learner"` in their signatures. The
  outcome-specific signatures show `outcome_method = "logistic"`,
  `"poisson"`, and `"coxph"`, respectively. Specify `outcome_method` for
  the doubly robust methods.

- [`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
  and its outcome-specific interfaces accept a one-sided `benefit_model`
  formula for two-regression and contrast-regression scores. It defaults
  to the terms of `outcome_model`, which defines the outcome prediction
  models. T-learners continue to use `outcome_model`. The previous
  `init.model` argument remains accepted through `...` under its
  original meaning and emits a deprecation warning. It cannot be
  combined with `benefit_model`. Coefficient tables for mixed models may
  contain missing entries for terms a score does not use.

- Use `fractional_logit` for survival initial outcome models in
  two-regression and contrast regression. The former value `logistic`
  remains accepted with a deprecation warning, including in
  `initial.predictor.method`. The underlying calculation is unchanged.

- Reorganized the binary, count and survival CATE vignettes to introduce
  model formulas, score-construction methods and fitted examples in
  sequence.

- [`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
  and its outcome-specific interfaces now expose `ps_model` and
  `ps_method` after `benefit_method`. The dotted names remain accepted
  through `...` with deprecation warnings; existing positional CATE
  calls remain supported.

- Outcome-specific
  [`catefitbinary()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitbinary.md),
  [`catefitcount()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitcount.md),
  [`catefitsurv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitsurv.md),
  and
  [`catefitmean()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitmean.md)
  accept `outcome_model`, `outcome_method`, and `benefit_method` in the
  same order as
  [`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md).
  Existing dotted arguments and positional calls remain supported.
  Continuous-outcome fitting remains unimplemented.

- [`predict()`](https://rdrr.io/r/stats/predict.html) for `catefit`
  objects accepts `"two_regression"`, `"contrast_regression"`,
  `"contrastreg"`, `"negative_binomial"`, and `"random_forest"`
  alongside the score names. `"t_learner"` selects a T-learner when
  exactly one is present; otherwise select it by its outcome method.
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html), and
  [`subgroup_estimates()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/subgroup_estimates.md)
  accept the same selectors, subject to coefficient and inference
  availability. Prediction still returns scores for the fitting sample.

- [`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
  accepts the established method names `"twoReg"` and `"contrastReg"` in
  `benefit_method` as aliases for `"two_regression"` and
  `"contrast_regression"`.

- Clarified the package description and overview to identify AIPW ATE
  estimators, T-learner CATE scores, and the associated methodological
  work.

- Reordered the vignette navigation to introduce ATE methods and
  examples before CATE methods and examples. The CATE methods vignette
  now begins with a cross-outcome introduction to CATE scores,
  T-learners, and validation. The binary and count CATE guides explain
  the outcome and treatment formulas alongside score fitting.

- [`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
  accepts paired `benefit_method` and `outcome_method` vectors, with
  `outcome_model` for the common outcome formula. The previous named
  `score.method` and `cate.model` arguments remain available through
  `...` with warnings, and earlier positional calls remain supported.
  The new arguments follow `response` and `data` in the order
  `outcome_model`, `outcome_method`, `benefit_model`, and
  `benefit_method`. The two doubly robust procedures can use different
  initial outcome methods; identical methods share cross-fitted
  predictions, and treatment-model predictions are shared across both.

- Binary
  [`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
  and
  [`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)
  now fit separate logistic outcome models by treatment arm for
  `score.method = "logistic"`. Stored fits contain one model per arm;
  fitted scores and validation results can change, including when an
  arm-specific design is rank deficient.

- Survival fitting with
  [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md),
  [`atefitsurv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefitsurv.md),
  [`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md),
  and
  [`catefitsurv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitsurv.md)
  now uses `censoring_model` and `censoring_method`. The former
  `ipcw.model` and `ipcw.method` names remain available through `...`
  with deprecation warnings; supplying an old and new name together is
  an error. Positional argument order is unchanged.

- Consolidated guidance on `atefit` result fields and methods in
  [`?atefit`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md),
  with outcome-specific help pages linking to it.

- [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
  and its outcome-specific functions now use `ps_model` and `ps_method`
  for the treatment model and its fitting method. The former `ps.model`
  and `ps.method` names remain available through `...` with deprecation
  warnings; supplying an old and new name together is an error.

- Split the binary, count and survival worked examples into separate ATE
  and CATE vignettes. Existing example URLs continue to lead to the ATE
  guides; the README and vignette navigation link to both guides for
  each outcome.

- Survival
  [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
  supports `outcome_method = "weibull"` and `"log_normal"` (alias
  `"lognormal"`) in addition to the default Cox outcome model.
  Arm-specific accelerated failure time models supply RMST and horizon
  survival predictions. The Weibull model allows different shapes by
  arm; the log-normal model does not impose proportional hazards.
  Censoring models remain controlled separately by `censoring_method`
  and `censoring_model`.

- [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
  now accepts `outcome_method`: count outcomes use `"poisson"` by
  default and also support `"negative_binomial"` (alias `"negbin"`).
  Binary and survival outcomes accept `"logistic"` and `"coxph"`,
  respectively. Negative binomial outcome models are refitted in each
  bootstrap sample; their estimates and uncertainty can differ from
  Poisson-based results. The new argument follows `outcome_model`; calls
  placing the treatment-model formula immediately after the outcome
  formula remain supported.

- [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
  and outcome-specific ATE functions now use `outcome_model` for the
  outcome formula. The former named argument `cate.model` remains a
  deprecated alias and emits a warning; supplying both arguments is an
  error. The alias is accepted through `...` rather than listed as a
  formal argument. The help describes how the outcome models handle
  treatment. See the major change above for the removal of
  `interactions` and its positional impact.

- `init.model` can specify separate initial outcome-model covariates for
  `twoReg` and `contrastReg` across supported outcome types in fitting
  and validation. It defaults to `cate.model`. Continuous-outcome
  fitting and validation remain unavailable.

- Boosting performance plots are now disabled by default in model
  fitting and cross-validation. Set `plot.gbmperf = TRUE` to display
  them.

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

- Count and survival preprocessing warnings no longer print internal
  calls that expose deprecated argument names in vignette output.

- Validation accepts `prop.cutoff = 1` with `abc = FALSE` when only the
  overall nested subgroup estimate is requested.

- Corrected and reordered the CATE fitting help arguments to match the
  function signatures.

- Survival
  [`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
  now reports a clear error when a treatment arm is missing.

- Fixed
  [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
  argument evaluation for count and survival data supplied from a local
  calling environment.

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
