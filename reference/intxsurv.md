# Estimate the CATE model using specified scoring methods for survival outcomes

Coefficients of the CATE estimated with random forest, boosting, naive
Poisson, two regression, and contrast regression

## Usage

``` r
intxsurv(
  y,
  d,
  trt,
  x.cate,
  x.ps,
  x.ipcw,
  yf = NULL,
  tau0,
  surv.min = 0.025,
  score.method = c("randomForest", "boosting", "poisson", "twoReg", "contrastReg"),
  ps.method = "glm",
  minPS = 0.01,
  maxPS = 0.99,
  ipcw.method = "breslow",
  initial.predictor.method = "randomForest",
  tree.depth = 3,
  n.trees.rf = 1000,
  n.trees.boosting = 150,
  B = 3,
  Kfold = 5,
  plot.gbmperf = TRUE,
  error.maxNR = 0.001,
  max.iterNR = 100,
  tune = c(0.5, 2),
  contrast.control = list(),
  ps.control = list(),
  ...
)
```

## Arguments

- y:

  Observed survival or censoring time; vector of size `n`.

- d:

  The event indicator, normally `1 = event, 0 = censored`; vector of
  size `n`.

- trt:

  Treatment received; vector of size `n` with treatment coded as 0/1.

- x.cate:

  Matrix of `p.cate` baseline covariates specified in the outcome model;
  dimension `n` by `p.cate`.

- x.ps:

  Matrix of `p.ps` baseline covariates specified in the propensity score
  model; dimension `n` by `p.ps`.

- x.ipcw:

  Matrix of `p.ipw` baseline covariate specified in inverse probability
  of censoring weighting; dimension `n` by `p.ipw`.

- yf:

  Follow-up time, interpreted as the potential censoring time; vector of
  size `n` if the potential censoring time is known.

- tau0:

  The truncation time for defining restricted mean time lost.

- surv.min:

  Lower truncation limit for probability of being censored (positive and
  very close to 0).

- score.method:

  A vector of one or multiple methods to estimate the CATE score.
  Allowed values are: `'randomForest'`, `'boosting'`, `'poisson'`,
  `'twoReg'`, `'contrastReg'`. Default specifies all 5 methods.

- ps.method:

  A character vector for the method to estimate the propensity score.
  Allowed values include one of: `'glm'` for logistic regression with
  main effects only (default), or `'lasso'` for a logistic regression
  with main effects and LASSO penalization on two-way interactions
  (added to the model if interactions are not specified in `ps.model`).
  Relevant only when `ps.model` has more than one variable.

- minPS:

  A numerical value (in \`\[0, 1\]\`) below which estimated propensity
  scores should be truncated. Default is `0.01`.

- maxPS:

  A number above which estimated propensity scores should be trimmed;
  scalar

- ipcw.method:

  The censoring model. Allowed values are: `'breslow'` (Cox regression
  with Breslow estimator of the baseline survivor function),
  `'aft (exponential)'`, `'aft (weibull)'`, `'aft (lognormal)'` or
  `'aft (loglogistic)'`. Default is `'breslow'`.

- initial.predictor.method:

  A character vector for the method used to get initial outcome
  predictions conditional on the covariates in `cate.model` in
  `score.method = 'twoReg'` and `'contrastReg'`. Allowed values include
  one of `'randomForest'`, `'boosting'` and `'logistic'` (fastest).
  Default is `'randomForest'`.

- tree.depth:

  A positive integer specifying the depth of individual trees in
  boosting (usually 2-3). Used only if `score.method = 'boosting'` or if
  `score.method = 'twoReg'` or `'contrastReg'` and
  `initial.predictor.method = 'boosting'`. Default is `3`.

- n.trees.rf:

  A positive integer specifying the maximum number of trees in random
  forest. Used if `score.method = 'ranfomForest'` or if
  `initial.predictor.method = 'randomForest'` with
  `score.method = 'twoReg'` or `'contrastReg'`. Default is `1000`.

- n.trees.boosting:

  A positive integer specifying the maximum number of trees in boosting
  (usually 100-1000). Used if `score.method = 'boosting'` or if
  `initial.predictor.method = 'boosting'` with `score.method = 'twoReg'`
  or `'contrastReg'`. Default is `150`.

- B:

  A positive integer specifying the number of time cross-fitting is
  repeated in `score.method = 'twoReg'` and `'contrastReg'`. Default is
  `3`.

- Kfold:

  A positive integer specifying the number of folds (parts) used in
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
  Default is `100`.

- tune:

  A vector of 2 numerical values \> 0 specifying tuning parameters for
  the Newton Raphson algorithm. `tune[1]` is the step size, `tune[2]`
  specifies a quantity to be added to diagonal of the slope matrix to
  prevent singularity. Used only if `score.method = 'contrastReg'`.
  Default is `c(0.5, 2)`.

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
  affects screening when `ps.method = "lasso"`; the LASSO fit and
  interaction selection rule are unchanged. Kept separate from
  boosting's `...` and other numerical controls.

- ...:

  Additional boosting arguments; see
  [`boosting-control`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/boosting-control.md).

## Value

Depending on what score.method is, the outputs is a combination of the
following: result.randomForest: Results of random forest fit, for trt =
0 and trt = 1 separately result.boosting: Results of boosting fit, for
trt = 0 and trt = 1 separately result.poisson: Naive Poisson estimator
(beta1 - beta0); vector of length `p.cate` + 1 result.twoReg: Two
regression estimator (beta1 - beta0); vector of length `p.cate` + 1
result.contrastReg: A list of the contrast regression results with 2
elements: \$delta.contrastReg: Contrast regression DR estimator; vector
of length `p.cate` + 1 \$converge.contrastReg: Indicator that the Newton
Raphson algorithm converged for `delta_0`; boolean
