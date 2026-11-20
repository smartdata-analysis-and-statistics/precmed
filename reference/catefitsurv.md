# Estimation of the conditional average treatment effect (CATE) score for survival data

Provides singly robust and doubly robust estimation of CATE score for
survival data with up to 5 scoring methods among the following: Random
forest, boosting, poisson regression, two regressions, and contrast
regression.

## Usage

``` r
catefitsurv(
  data,
  outcome_model = NULL,
  outcome_method = "coxph",
  benefit_model = NULL,
  benefit_method = "t_learner",
  ps_model,
  ps_method = "glm",
  censoring_model = NULL,
  censoring_method = "breslow",
  minPS = 0.01,
  maxPS = 0.99,
  followup.time = NULL,
  tau0 = NULL,
  higher.y = TRUE,
  prop.cutoff = seq(0.5, 1, length = 6),
  surv.min = 0.025,
  tree.depth = 2,
  n.trees.rf = 1000,
  n.trees.boosting = 200,
  B = 3,
  Kfold = 5,
  plot.gbmperf = FALSE,
  error.maxNR = 0.001,
  max.iterNR = 100,
  tune = c(0.5, 2),
  seed = NULL,
  verbose = 0,
  rmst.control = list(),
  contrast.control = list(),
  ps.control = list(),
  rf.control = list(),
  ...
)
```

## Arguments

- data:

  A data frame containing the variables in the outcome, treatment, and
  censoring models (if specified); a data frame with `n` rows (1 row per
  observation).

- outcome_model:

  A formula specifying the outcome and baseline covariates.

- outcome_method:

  Outcome-model method paired with `benefit_method`. For a T-learner,
  `"coxph"` (default), `"weibull"`, `"log_normal"` (alias
  `"lognormal"`), `"poisson"`, `"boosting"`, or `"random_forest"`. For
  two-regression or contrast regression, `"fractional_logit"`,
  `"boosting"`, or `"random_forest"`; the older value `"logistic"` is
  accepted for fractional logit with a deprecation warning. The Cox,
  Weibull and log-normal T-learners compare arm-specific predicted
  restricted mean time lost at `tau0`. Cox assumes proportional hazards
  and estimates a baseline cumulative hazard; Weibull and log-normal are
  parametric accelerated failure time models. Causal interpretation
  requires appropriate model, identification and censoring assumptions.
  These three T-learners are not yet available in
  [`catecv`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md).
  Two-regression and contrast regression may use different initial
  outcome methods; matching methods share cross-fitted outcome
  predictions.

- benefit_model:

  Optional one-sided formula for the two-regression or
  contrast-regression score covariates. Defaults to the covariates and
  terms in `outcome_model`; T-learner scores do not use it. See
  [`catefit`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md).

- benefit_method:

  Score-construction methods; defaults to `"t_learner"`. See
  [`catefit`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md).

- ps_model:

  A formula for the treatment model, which estimates propensity scores.
  Put the treatment variable on the left-hand side; for an
  intercept-only treatment model, use for example `trt ~ 1`. See
  [`catefit`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md).

- ps_method:

  A character value for the method to estimate the propensity score.
  Allowed values include one of: `'glm'` for logistic regression with
  main effects only (default), or `'lasso'` for a logistic regression
  with main effects and LASSO penalization on two-way interactions
  (added to the model if interactions are not specified in `ps_model`).
  Relevant only when `ps_model` has more than one variable.

- censoring_model:

  A one-sided formula for the model that estimates probabilities of
  remaining uncensored for inverse probability of censoring weighting
  (IPCW). By default (`NULL`), the model uses the outcome-model
  covariates and treatment.

- censoring_method:

  Method for fitting the censoring model. The default `"breslow"` fits
  Cox regression with the Breslow estimator of the baseline survivor
  function. Accelerated failure time options are `"aft (exponential)"`,
  `"aft (weibull)"`, `"aft (lognormal)"`, and `"aft (loglogistic)"`.

- minPS:

  A numerical value (in \`\[0, 1\]\`) below which estimated propensity
  scores should be truncated. Default is `0.01`.

- maxPS:

  A numerical value (in \`(0, 1\]\`) above which estimated propensity
  scores should be truncated. Must be strictly greater than `minPS`.
  Default is `0.99`.

- followup.time:

  A column name in `data` specifying the maximum follow-up time,
  interpreted as the potential censoring time. Default is
  `followup.time = NULL`, which corresponds to unknown potential
  censoring time.

- tau0:

  The truncation time for defining restricted mean time lost. Default is
  `NULL`, which corresponds to setting the truncation time as the
  maximum survival time in the data.

- higher.y:

  A logical value indicating whether higher (`TRUE`) or lower (`FALSE`)
  values of the outcome are more desirable. Default is `TRUE`.

- prop.cutoff:

  A vector of numerical values (in \`(0, 1\]\`) specifying percentiles
  of the estimated log CATE scores to define nested subgroups. Each
  element represents the cutoff to separate observations in nested
  subgroups (below vs above cutoff). The length of `prop.cutoff` is the
  number of nested subgroups. An equally-spaced sequence of proportions
  ending with 1 is recommended. Default is `seq(0.5, 1, length = 6)`.

- surv.min:

  Lower truncation limit for the probability of remaining uncensored,
  used in IPCW. It must be a positive value and should be chosen close
  to 0. Default is `0.025`.

- tree.depth:

  A positive integer specifying the maximum depth of boosting trees.
  Used when `outcome_method = 'boosting'`. Default is `2`.

- n.trees.rf:

  A positive integer specifying the number of trees in random forest.
  Used with `outcome_method = 'random_forest'` for survival outcomes.
  Default is `1000`.

- n.trees.boosting:

  The maximum number of boosting trees, used when
  `outcome_method = 'boosting'`. Default is `200`.

- B:

  Number of repetitions of cross-fitting for
  `benefit_method = 'two_regression'` or `'contrast_regression'`.
  Default is `3`.

- Kfold:

  A positive integer specifying the number of folds used in
  cross-fitting to partition the data in
  `benefit_method = 'two_regression'` and `'contrast_regression'`.
  Default is `5`.

- plot.gbmperf:

  Whether to plot boosting performance measures when
  `outcome_method = 'boosting'`. Default is `FALSE`.

- error.maxNR:

  Convergence tolerance for the Newton-Raphson algorithm used by
  `benefit_method = 'contrast_regression'`. Default is `0.001`.

- max.iterNR:

  A positive integer indicating the maximum number of iterations in the
  Newton Raphson algorithm. Used only if
  `benefit_method = 'contrast_regression'`. Default is `100`.

- tune:

  A vector of 2 numerical values \> 0 specifying tuning parameters for
  the Newton Raphson algorithm. `tune[1]` is the step size, `tune[2]`
  specifies a quantity to be added to diagonal of the slope matrix to
  prevent singularity. Used only if
  `benefit_method = 'contrast_regression'`. Default is `c(0.5, 2)`.

- seed:

  An optional integer specifying an initial randomization seed for
  reproducibility. Default is `NULL`, corresponding to no seed.

- verbose:

  An integer value indicating what kind of intermediate progress
  messages should be printed. `0` means no outputs. `1` means only
  progress and run time. `2` means progress, run time, and all errors
  and warnings. Default is `0`.

- rmst.control:

  Named options for Cox restricted mean survival integration. Supported
  entry: `solver = "C"` (default) or `solver = "R"`. The C backend
  accumulates each patient's integral without allocating a
  patient-by-event-time matrix. Cox model fitting is unchanged.
  Available only for survival outcomes; independent of
  `contrast.control`.

- contrast.control:

  A named list of contrast-regression options: `solver` (`"C"`, the
  default, or `"R"`), `error.maxNR`, `max.iterNR`, `tune`, and
  `optim.control` (a named list passed to
  [`stats::optim`](https://rdrr.io/r/stats/optim.html) on fallback).
  List entries override top-level NR arguments. This is independent of
  `rmst.control`, which selects RMST integration, and is not passed to
  boosting via `...`.

- ps.control:

  Named options for LASSO propensity-score interaction screening.
  Supported entry: `solver = "C"` (default) or `solver = "R"`. Only
  affects screening when `ps_method = "lasso"`; the LASSO fit and
  interaction selection rule are unchanged. Kept separate from
  boosting's `...` and other numerical controls.

- rf.control:

  Named list of randomForestSRC options: `mtry`, `nodesize`, `nsplit`,
  `nodedepth`, and `splitrule`. Applies to both treatment arms and
  survival initial-predictor forests. Use `n.trees.rf` for the number of
  trees.

- ...:

  Additional boosting options. See
  [`boosting-control`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/boosting-control.md).
  Deprecated argument names from precmed 1.1.0 remain supported for
  backward compatibility and emit warnings.

## Value

Returns an object of the class `catefit` containing the following
components: When the doubly robust scores use distinct outcome methods,
`outcome_methods` records the method selected for each score.

- `ate.randomForest`: A vector of numerical values of length
  `prop.cutoff` containing the estimated ATE by the RMTL ratio in nested
  subgroups (defined by `prop.cutoff`) constructed based on the
  estimated CATE scores with random forest method. Only provided if the
  fit includes a random forest.

- `ate.boosting`: Same as `ate.randomForest`, but with the nested
  subgroups based the estimated CATE scores with boosting. Only provided
  if the fit includes boosting.

- `ate.poisson`: Same as `ate.randomForest`, but with the nested
  subgroups based the estimated CATE scores with poisson regression.
  Only provided if the fit includes Poisson regression.

- `ate.twoReg`: Same as `ate.randomForest`, but with the nested
  subgroups based the estimated CATE scores with two regressions. Only
  provided if the fit includes two-regression.

- `ate.contrastReg`: Same as `ate.randomForest`, but with the nested
  subgroups based the estimated CATE scores with contrast regression.
  Only provided if the fit includes contrast regression.

- `chr.<method>`: A vector of subgroup marginal cumulative hazard ratios
  at `tau0`, for each selected scoring method. Undefined estimates are
  `NA`; existing `hr.*` results are conditional Cox HRs.

- `hr.randomForest`: A vector of numerical values of length
  `prop.cutoff` containing the adjusted hazard ratio in nested subgroups
  (defined by `prop.cutoff`) constructed based on the estimated CATE
  scores with random forest method. Only provided if the fit includes a
  random forest.

- `hr.boosting`: Same as `hr.randomForest`, but with the nested
  subgroups based the estimated CATE scores with boosting. Only provided
  if the fit includes boosting.

- `hr.poisson`: Same as `hr.randomForest`, but with the nested subgroups
  based the estimated CATE scores with poisson regression. Only provided
  if the fit includes Poisson regression.

- `hr.twoReg`: Same as `hr.randomForest`, but with the nested subgroups
  based the estimated CATE scores with two regressions. Only provided if
  the fit includes two-regression.

- `hr.contrastReg`: Same as `hr.randomForest`, but with the nested
  subgroups based the estimated CATE scores with contrast regression.
  Only provided if the fit includes contrast regression.

- `score.randomForest`: A vector of numerical values of length n (number
  of observations in `data`) containing the estimated log-CATE scores
  according to random forest. Only provided if the fit includes a random
  forest.

- `score.boosting`: Same as `score.randomForest`, but with estimated
  log-CATE score according to boosting. Only provided if the fit
  includes boosting.

- `score.poisson`: Same as `score.randomForest`, but with estimated
  log-CATE score according to the Poisson regression. Only provided if
  the fit includes Poisson regression.

- `score.twoReg`: Same as `score.randomForest`, but with estimated
  log-CATE score according to two regressions. Only provided if the fit
  includes two-regression.

- `score.contrastReg`: Same as `score.randomForest`, but with estimated
  log-CATE score according to contrast regression. Only provided if the
  fit includes contrast regression.

- `fit`: Additional model-fitting details for random forests, boosting,
  and contrast regression.

  - `result.randomForest`: Details on the random forest model fitted to
    observations with treatment = 0 `($fit0.rf)` and to observations
    with treatment = 1 `($fit1.rf)`. Only provided if the fit includes a
    random forest.

  - `result.boosting`: Details on the boosting model fitted to
    observations with treatment = 0, `($fit0.boosting)` and
    `($fit0.gam)` and to observations with treatment = 1,
    `($fit1.boosting)` and `($fit1.gam)`. Only provided if the fit
    includes boosting.

  - `result.contrastReg$converge.contrastReg`: Whether the contrast
    regression algorithm converged or not. Only provided if the fit
    includes contrast regression.

- `coefficients`: A data frame with the coefficients of the estimated
  log-CATE score for each applicable method. The data frame has number
  of rows equal to the number of design columns in `outcome_model`,
  including the intercept, and one column per available linear scoring
  method. Standard errors are not calculated for survival scores.
  `'randomForest'` and `'boosting'` do not have coefficient results
  because tree-based methods typically do not express the log-CATE as a
  linear combination of coefficients and covariates.

- `risk0.<method>`, `risk1.<method>`, `rd.<method>`, `rmst0.<method>`,
  `rmst1.<method>`: Subgroup event risks, risk differences (treatment 1
  minus 0), and treatment-specific RMST at `tau0`. See
  [`subgroup_estimates`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/subgroup_estimates.md)
  for extraction.

- `errors/warnings`: A nested list of errors and warnings that were
  wrapped during the calculation of ATE. Errors and warnings are
  organized by scoring method.

## Details

The CATE score estimates a covariate-conditional average treatment
effect for survival data, estimated with random forest, boosting,
Poisson regression, and the doubly robust estimator (two regressions,
Yadlowsky, 2020) applied separately by treatment group or with the other
doubly robust estimators (contrast regression, Yadlowsky, 2020) applied
to the entire data set.

`catefitsurv()` provides CATE score coefficients for supported linear
scoring methods. The coefficients describe the fitted log-RMTL score;
coefficient standard errors and confidence intervals are not currently
calculated for survival fits.

`catefitsurv()` also provides the predicted CATE score of each
observation in the data set, for each scoring method. The predictions
allow ranking the observations from potentially high responders to the
treatment to potentially low or standard responders.

The estimated ATE among nested subgroups of high responders are also
provided by scoring method. Note that the ATEs in `catefitsurv()` are
derived based on the CATE score which is estimated using the same data
sample. Therefore, overfitting may be an issue.
[`catecvsurv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecvsurv.md)
is more suitable to inspect the estimated ATEs across scoring methods as
it implements internal cross validation to reduce optimism.

## References

Yadlowsky, S., Pellegrini, F., Lionetto, F., Braune, S., & Tian, L.
(2020). *Estimation and validation of ratio-based conditional average
treatment effects using observational data. Journal of the American
Statistical Association, 1-18.* DOI: 10.1080/01621459.2020.1772080.

## See also

[`catecvsurv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecvsurv.md)

## Examples

``` r
# \donttest{
library(survival)

tau0 <- with(survivalExample, min(quantile(y[trt == "drug1"], 0.95),
                               quantile(y[trt == "drug0"], 0.95)))

fit <- catefitsurv(data = survivalExample,
                   outcome_model = Surv(y, d) ~ age + female + previous_cost +
                                             previous_number_relapses,
                   outcome_method = "random_forest",
                   benefit_method = "t_learner",
                   ps_model = trt ~ age + previous_treatment,
                   censoring_model = ~ age + previous_cost + previous_treatment,
                   tau0 = tau0,
                   seed = 999)
#> Variable trt was recoded to 0/1 with drug0->0 and drug1->1.

coef(fit)
#> NULL
# }
```
