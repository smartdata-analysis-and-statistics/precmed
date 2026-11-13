# Estimate the CATE model using specified scoring methods

Coefficients of the CATE estimated with boosting, linear regression, two
regression, contrast regression, random forest, generalized additive
model

## Usage

``` r
intxmean(
  y,
  trt,
  x.cate,
  x.init,
  x.ps,
  score.method = c("boosting", "gaussian", "twoReg", "contrastReg", "gam",
    "randomForest"),
  ps.method = "glm",
  minPS = 0.01,
  maxPS = 0.99,
  initial.predictor.method = "boosting",
  xvar.smooth.init,
  xvar.smooth.score,
  tree.depth = 2,
  n.trees.rf = 1000,
  n.trees.boosting = 200,
  B = 1,
  Kfold = 2,
  plot.gbmperf = TRUE,
  ps.control = list(),
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

- x.init:

  Matrix of `p.init` baseline covariates; dimension `n` by `p.init` It
  must be specified when `score.method = contrastReg` or `twoReg`.

- x.ps:

  Matrix of `p.ps` baseline covariates (plus a leading column of 1 for
  the intercept); dimension `n` by `p.ps + 1` (covariates in the
  propensity score model plus intercept)

- score.method:

  A vector of one or multiple methods to estimate the CATE score.
  Allowed values are: `'boosting'`, `'gaussian'`, `'twoReg'`,
  `'contrastReg'`, `'randomForest'`, `'gam'`. Default specifies all 6
  methods.

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
  one of `'gaussian'` (fastest), `'boosting'` (default) and `'gam'`.

- xvar.smooth.init:

  A vector of characters indicating the name of the variables used as
  the smooth terms if `initial.predictor.method = 'gam'`. The variables
  must be selected from the variables listed in `init.model`. Default is
  `NULL`, which uses all variables in `init.model`.

- xvar.smooth.score:

  A vector of characters indicating the name of the variables used as
  the smooth terms if `score.method = 'gam'`. The variables must be
  selected from the variables listed in `cate.model`. Default is `NULL`,
  which uses all variables in `cate.model`.

- tree.depth:

  A positive integer specifying the depth of individual trees in
  boosting (usually 2-3). Used only if `score.method = 'boosting'` or if
  `score.method = 'twoReg'` or `'contrastReg'` and
  `initial.predictor.method = 'boosting'`. Default is `2`.

- n.trees.rf:

  A positive integer specifying the number of trees. Used only if
  `score.method = 'randomForest'`. Default is `1000`.

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

- ps.control:

  Named options for LASSO propensity-score interaction screening.
  Supported entry: `solver = "C"` (default) or `solver = "R"`. Only
  affects screening when `ps.method = "lasso"`; the LASSO fit and
  interaction selection rule are unchanged. Kept separate from
  boosting's `...` and other numerical controls.

- ...:

  Additional arguments for `gbm()`

## Value

Depending on what score.method is, the outputs is a combination of the
following: result.boosting: Results of boosting fit and best iteration,
for trt = 0 and trt = 1 separately result.gaussian: Linear regression
estimator (beta1 - beta0); vector of length `p.cate` + 1 result.twoReg:
Two regression estimator (beta1 - beta0); vector of length `p.cate` + 1
result.contrastReg: A list of the contrast regression results with 3
elements: \$delta.contrastReg: Contrast regression DR estimator; vector
of length `p.cate` + 1 \$sigma.contrastReg: Variance covariance matrix
for delta.contrastReg; matrix of size `p.cate` + 1 by `p.cate` + 1
result.randomForest: Results of random forest fit and best iteration,
for trt = 0 and trt = 1 separately result.gam: Results of generalized
additive model fit and best iteration, for trt = 0 and trt = 1
separately best.iter: Largest best iterations for boosting (if used)
fgam: Formula applied in GAM when `initial.predictor.method = 'gam'`
warn.fit: Warnings occurred when fitting `score.method` err.fit:: Errors
occurred when fitting `score.method`
