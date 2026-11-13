# Estimate the CATE model using specified scoring methods

Coefficients of the CATE estimated with boosting, naive Poisson, two
regression, contrast regression, negative binomial

## Usage

``` r
intxcount(
  y,
  trt,
  x.cate,
  x.ps,
  time,
  score.method = c("boosting", "poisson", "twoReg", "contrastReg", "negBin"),
  ps.method = "glm",
  minPS = 0.01,
  maxPS = 0.99,
  initial.predictor.method = "boosting",
  xvar.smooth = NULL,
  tree.depth = 2,
  n.trees.boosting = 200,
  B = 3,
  Kfold = 6,
  plot.gbmperf = TRUE,
  error.maxNR = 0.001,
  max.iterNR = 150,
  tune = c(0.5, 2),
  contrast.control = list(),
  ...
)
```

## Arguments

- y:

  Observed outcome; vector of size `n` (observations)

- trt:

  Treatment received; vector of size `n` units with treatment coded as
  0/1

- x.cate:

  Matrix of `p.cate` baseline covariates; dimension `n` by `p.cate`
  (covariates in the outcome model)

- x.ps:

  Matrix of `p.ps` baseline covariates (plus a leading column of 1 for
  the intercept); dimension `n` by `p.ps + 1` (covariates in the
  propensity score model plus intercept)

- time:

  Log-transformed person-years of follow-up; vector of size `n`

- score.method:

  A vector of one or multiple methods to estimate the CATE score.
  Allowed values are: `'boosting'`, `'poisson'`, `'twoReg'`,
  `'contrastReg'`, `'negBin'`. Default specifies all 5 methods.

- ps.method:

  A character value for the method to estimate the propensity score.
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

- initial.predictor.method:

  A character vector for the method used to get initial outcome
  predictions conditional on the covariates in `cate.model` in
  `score.method = 'twoReg'` and `'contrastReg'`. Allowed values include
  one of `'poisson'` (fastest), `'boosting'` (default) and `'gam'`.

- xvar.smooth:

  A vector of characters indicating the name of the variables used as
  the smooth terms if `initial.predictor.method = 'gam'`. The variables
  must be selected from the variables listed in `cate.model`. Default is
  `NULL`, which uses all variables in `cate.model`.

- tree.depth:

  A positive integer specifying the depth of individual trees in
  boosting (usually 2-3). Used only if `score.method = 'boosting'` or if
  `score.method = 'twoReg'` or `'contrastReg'` and
  `initial.predictor.method = 'boosting'`. Default is `2`.

- n.trees.boosting:

  A positive integer specifying the maximum number of trees in boosting
  (usually 100-1000). Used only if `score.method = 'boosting'` or if
  `score.method = 'twoReg'` or `'contrastReg'` and
  `initial.predictor.method = 'boosting'`. Default is `200`.

- B:

  A positive integer specifying the number of time cross-fitting is
  repeated in `score.method = 'twoReg'` and `'contrastReg'`. Default is
  `3`.

- Kfold:

  A positive integer specifying the number of folds (parts) used in
  cross-fitting to partition the data in `score.method = 'twoReg'` and
  `'contrastReg'`. Default is `6`.

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

- contrast.control:

  A named list of options for count-outcome contrast regression.
  Supported names are `solver` (`"C"`, the default, or `"R"`),
  `error.maxNR`, `max.iterNR`, `tune`, and `optim.control` (a named list
  passed to [`stats::optim`](https://rdrr.io/r/stats/optim.html) on
  fallback). List entries override the corresponding top-level NR
  arguments. Options are kept separate from `...`, which is forwarded to
  boosting. Only count outcomes support this list; other outcomes retain
  their existing solver.

- ...:

  Additional arguments for `gbm()`

## Value

Depending on what score.method is, the outputs is a combination of the
following: result.boosting: Results of boosting fit and best iteration,
for trt = 0 and trt = 1 separately result.poisson: Naive Poisson
estimator (beta1 - beta0); vector of length `p.cate` + 1 result.twoReg:
Two regression estimator (beta1 - beta0); vector of length `p.cate` + 1
result.contrastReg: A list of the contrast regression results with 3
elements: \$delta.contrastReg: Contrast regression DR estimator; vector
of length `p.cate` + 1 \$sigma.contrastReg: Variance covariance matrix
for delta.contrastReg; matrix of size `p.cate` + 1 by `p.cate` + 1
\$converge.contrastReg: Indicator that the Newton Raphson algorithm
converged for `delta_0`; boolean result.negBin: Negative binomial
estimator (beta1 - beta0); vector of length `p.cate` + 1 best.iter:
Largest best iterations for boosting (if used) fgam: Formula applied in
GAM (if used)
