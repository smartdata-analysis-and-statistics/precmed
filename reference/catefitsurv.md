# Estimation of the conditional average treatment effect (CATE) score for survival data

Provides singly robust and doubly robust estimation of CATE score for
survival data with up to 5 scoring methods among the following: Random
forest, boosting, poisson regression, two regressions, and contrast
regression.

## Usage

``` r
catefitsurv(
  data,
  score.method,
  cate.model,
  ps.model,
  ps.method = "glm",
  initial.predictor.method = "randomForest",
  ipcw.model = NULL,
  ipcw.method = "breslow",
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
  plot.gbmperf = TRUE,
  error.maxNR = 0.001,
  max.iterNR = 100,
  tune = c(0.5, 2),
  seed = NULL,
  verbose = 0,
  ...
)
```

## Arguments

- data:

  A data frame containing the variables in the outcome, propensity
  score, and inverse probability of censoring models (if specified); a
  data frame with `n` rows (1 row per observation).

- score.method:

  A vector of one or multiple methods to estimate the CATE score.
  Allowed values are: `'randomForest'`, `'boosting'`, `'poisson'`,
  `'twoReg'`, and `'contrastReg'`.

- cate.model:

  A standard `Surv` formula describing the outcome model to be fitted.
  The outcome must appear on the left-hand side.

- ps.model:

  A formula describing the propensity score (PS) model to be fitted. The
  treatment must appear on the left-hand side. The treatment must be a
  numeric vector coded as 0/1. If data are from a randomized controlled
  trial, specify `ps.model = ~1` as an intercept-only model.

- ps.method:

  A character value for the method to estimate the propensity score.
  Allowed values include one of: `'glm'` for logistic regression with
  main effects only (default), or `'lasso'` for a logistic regression
  with main effects and LASSO penalization on two-way interactions
  (added to the model if interactions are not specified in `ps.model`).
  Relevant only when `ps.model` has more than one variable.

- initial.predictor.method:

  A character vector for the method used to get initial outcome
  predictions conditional on the covariates specified in `cate.model`.
  Only applies when `score.method` includes `'twoReg'` or
  `'contrastReg'`. Allowed values include one of `'randomForest'`,
  `'boosting'` and `'logistic'` (fastest). Default is `'randomForest'`.

- ipcw.model:

  A formula describing the inverse probability of censoring weighting
  (IPCW) model to be fitted. The left-hand side must be empty. Default
  is `ipcw.model = NULL`, which corresponds to specifying the IPCW model
  with the same covariates as the outcome model `cate.model` plus the
  treatment.

- ipcw.method:

  A character value for the censoring model. Allowed values are:
  `'breslow'` (Cox regression with Breslow estimator of the baseline
  survivor function), `'aft (exponential)'`, `'aft (weibull)'`,
  `'aft (lognormal)'` or `'aft (loglogistic)'` (accelerated failure time
  model with different distributions for y variable). Default is
  `'breslow'`.

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

  Lower truncation limit for the probability of being censored. It must
  be a positive value and should be chosen close to 0. Default is
  `0.025`.

- tree.depth:

  A positive integer specifying the depth of individual trees in
  boosting (usually 2-3). Used only if `score.method = 'boosting'` or if
  `initial.predictor.method = 'boosting'` with `score.method = 'twoReg'`
  or `'contrastReg'`. Default is 2.

- n.trees.rf:

  A positive integer specifying the maximum number of trees in random
  forest. Used if `score.method = 'ranfomForest'` or if
  `initial.predictor.method = 'randomForest'` with
  `score.method = 'twoReg'` or `'contrastReg'`. Only applies for
  survival outcomes. Default is `1000`.

- n.trees.boosting:

  A positive integer specifying the maximum number of trees in boosting
  (usually 100-1000). Used if `score.method = 'boosting'` or if
  `initial.predictor.method = 'boosting'` with `score.method = 'twoReg'`
  or `'contrastReg'`. Default is `200`.

- B:

  A positive integer specifying the number of time cross-fitting is
  repeated in `score.method = 'twoReg'` and `'contrastReg'`. Default is
  `3`.

- Kfold:

  A positive integer specifying the number of folds used in
  cross-fitting to partition the data in `score.method = 'twoReg'` and
  `'contrastReg'`. Default is `5`.

- plot.gbmperf:

  A logical value indicating whether to plot the performance measures in
  boosting. Used only if `score.method = 'boosting'` or if
  `score.method = 'twoReg'` or `'contrastReg'` and
  `initial.predictor.method = 'boosting'`. Default is `TRUE`.

- error.maxNR:

  A numerical value \> 0 indicating the minimum value of the mean
  absolute error in Newton Raphson algorithm. Used only if
  `score.method = 'contrastReg'`. Default is `0.001`.

- max.iterNR:

  A positive integer indicating the maximum number of iterations in the
  Newton Raphson algorithm. Used only if `score.method = 'contrastReg'`.
  Default is `150`.

- tune:

  A vector of 2 numerical values \> 0 specifying tuning parameters for
  the Newton Raphson algorithm. `tune[1]` is the step size, `tune[2]`
  specifies a quantity to be added to diagonal of the slope matrix to
  prevent singularity. Used only if `score.method = 'contrastReg'`.
  Default is `c(0.5, 2)`.

- seed:

  An optional integer specifying an initial randomization seed for
  reproducibility. Default is `NULL`, corresponding to no seed.

- verbose:

  An integer value indicating what kind of intermediate progress
  messages should be printed. `0` means no outputs. `1` means only
  progress and run time. `2` means progress, run time, and all errors
  and warnings. Default is `0`.

- ...:

  Additional arguments for `gbm()`

## Value

Returns an object of the class `catefit` containing the following
components:

- `ate.randomForest`: A vector of numerical values of length
  `prop.cutoff` containing the estimated ATE by the RMTL ratio in nested
  subgroups (defined by `prop.cutoff`) constructed based on the
  estimated CATE scores with random forest method. Only provided if
  `score.method` includes `'randomForest'`.

- `ate.boosting`: Same as `ate.randomForest`, but with the nested
  subgroups based the estimated CATE scores with boosting. Only provided
  if `score.method` includes `'boosting'`.

- `ate.poisson`: Same as `ate.randomForest`, but with the nested
  subgroups based the estimated CATE scores with poisson regression.
  Only provided if `score.method` includes `'poisson'`.

- `ate.twoReg`: Same as `ate.randomForest`, but with the nested
  subgroups based the estimated CATE scores with two regressions. Only
  provided if `score.method` includes `'twoReg'`.

- `ate.contrastReg`: Same as `ate.randomForest`, but with the nested
  subgroups based the estimated CATE scores with contrast regression.
  Only provided if `score.method` includes `'contrastReg'`.

- `hr.randomForest`: A vector of numerical values of length
  `prop.cutoff` containing the adjusted hazard ratio in nested subgroups
  (defined by `prop.cutoff`) constructed based on the estimated CATE
  scores with random forest method. Only provided if `score.method`
  includes `'randomForest'`.

- `hr.boosting`: Same as `hr.randomForest`, but with the nested
  subgroups based the estimated CATE scores with boosting. Only provided
  if `score.method` includes `'boosting'`.

- `hr.poisson`: Same as `hr.randomForest`, but with the nested subgroups
  based the estimated CATE scores with poisson regression. Only provided
  if `score.method` includes `'poisson'`.

- `hr.twoReg`: Same as `hr.randomForest`, but with the nested subgroups
  based the estimated CATE scores with two regressions. Only provided if
  `score.method` includes `'twoReg'`.

- `hr.contrastReg`: Same as `hr.randomForest`, but with the nested
  subgroups based the estimated CATE scores with contrast regression.
  Only provided if `score.method` includes `'contrastReg'`.

- `score.randomForest`: A vector of numerical values of length n (number
  of observations in `data`) containing the estimated log-CATE scores
  according to random forest. Only provided if `score.method` includes
  `'randomForest'`.

- `score.boosting`: Same as `score.randomForest`, but with estimated
  log-CATE score according to boosting. Only provided if `score.method`
  includes `'boosting'`.

- `score.poisson`: Same as `score.randomForest`, but with estimated
  log-CATE score according to the Poisson regression. Only provided if
  `score.method` includes `'poisson'`.

- `score.twoReg`: Same as `score.randomForest`, but with estimated
  log-CATE score according to two regressions. Only provided if
  `score.method` includes `'twoReg'`.

- `score.contrastReg`: Same as `score.randomForest`, but with estimated
  log-CATE score according to contrast regression. Only provided if
  `score.method` includes `'contrastReg'`.

- `fit`: Additional details on model fitting if `score.method` includes
  'randomForest', 'boosting' or 'contrastReg':

  - `result.randomForest`: Details on the random forest model fitted to
    observations with treatment = 0 `($fit0.rf)` and to observations
    with treatment = 1 `($fit1.rf)`. Only provided if `score.method`
    includes `'randomForest'`.

  - `result.boosting`: Details on the boosting model fitted to
    observations with treatment = 0, `($fit0.boosting)` and
    `($fit0.gam)` and to observations with treatment = 1,
    `($fit1.boosting)` and `($fit1.gam)`. Only provided if
    `score.method` includes `'boosting'`.

  - `result.contrastReg$converge.contrastReg`: Whether the contrast
    regression algorithm converged or not. Only provided if
    `score.method` includes `'contrastReg'`.

- `coefficients`: A data frame with the coefficients of the estimated
  log-CATE score by `score.method`. The data frame has number of rows
  equal to the number of covariates in `cate.model` and number of
  columns equal to `length(score.method)`. If `score.method` includes
  `'contrastReg'`, the data frame has an additional column containing
  the standard errors of the coefficients estimated with contrast
  regression. `'randomForest'` and `'boosting'` do not have coefficient
  results because tree-based methods typically do not express the
  log-CATE as a linear combination of coefficients and covariates.

- `errors/warnings`: A nested list of errors and warnings that were
  wrapped during the calculation of ATE. Errors and warnings are
  organized by `score.method`.

## Details

The CATE score represents an individual-level treatment effect for
survival data, estimated with random forest, boosting, Poisson
regression, and the doubly robust estimator (two regressions, Yadlowsky,
2020) applied separately by treatment group or with the other doubly
robust estimators (contrast regression, Yadlowsky, 2020) applied to the
entire data set.

`catefitsurv()` provides the coefficients of the CATE score for each
scoring method requested through `score.method`. Currently, contrast
regression is the only method which allows for inference of the CATE
coefficients by providing standard errors of the coefficients. The
coefficients can be used to learn the effect size of each variable and
predict the CATE score for a new observation.

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
                   score.method = "randomForest",
                   cate.model = Surv(y, d) ~ age + female + previous_cost +
                                             previous_number_relapses,
                   ps.model = trt ~ age + previous_treatment,
                   ipcw.model = ~ age + previous_cost + previous_treatment,
                   tau0 = tau0,
                   seed = 999)
#> Warning: Variable trt was recoded to 0/1 with drug0->0 and drug1->1.

coef(fit)
#> NULL
# }
```
