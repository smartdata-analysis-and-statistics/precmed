# Estimation of the conditional average treatment effect (CATE) score for continuous data

Provides singly robust and doubly robust estimation of CATE score with
up to 6 scoring methods among the following: Linear regression,
boosting, two regressions, contrast regression, random forest and
generalized additive model.

## Usage

``` r
catefitmean(
  data,
  score.method,
  cate.model,
  ps.model,
  ps.method = "glm",
  init.model = NULL,
  initial.predictor.method = "boosting",
  minPS = 0.01,
  maxPS = 0.99,
  higher.y = TRUE,
  prop.cutoff = seq(0.5, 1, length = 6),
  xvar.smooth.score = NULL,
  xvar.smooth.init = NULL,
  tree.depth = 2,
  n.trees.rf = 1000,
  n.trees.boosting = 200,
  B = 3,
  Kfold = 6,
  plot.gbmperf = FALSE,
  error.maxNR = 0.001,
  tune = c(0.5, 2),
  seed = NULL,
  verbose = 0,
  ps.control = list(),
  ...
)
```

## Arguments

- data:

  A data frame containing the variables in the outcome and propensity
  score models; a data frame with `n` rows (1 row per observation).

- score.method:

  A vector of one or multiple methods to estimate the CATE score.
  Allowed values are: `'boosting'`, `'gaussian'`, `'twoReg'`,
  `'contrastReg'`, `'randomForest'`, `'gam'`.

- cate.model:

  A formula describing the outcome model to be fitted. The outcome must
  appear on the left-hand side.

- ps.model:

  A formula describing the propensity score model to be fitted. The
  treatment must appear on the left-hand side. The treatment must be a
  numeric vector coded as 0/1. If data are from a RCT, specify
  `ps.model` as an intercept-only model.

- ps.method:

  A character value for the method to estimate the propensity score.
  Allowed values include one of: `'glm'` for logistic regression with
  main effects only (default), or `'lasso'` for a logistic regression
  with main effects and LASSO penalization on two-way interactions
  (added to the model if not specified in `ps.model`).

- init.model:

  A formula describing the initial predictor model. The outcome must
  appear on the left-hand side. It must be specified when
  `score.method = contrastReg` or `twoReg`.

- initial.predictor.method:

  A character vector for the method used to get initial outcome
  predictions conditional on the covariates in `init.model` in
  `score.method = 'twoReg'` and `'contrastReg'`. Allowed values include
  one of `'gaussian'` (fastest), `'boosting'` and `'gam'`. Default is
  `'boosting'`.

- minPS:

  A numerical value (in \`\[0, 1\]\`) below which estimated propensity
  scores should be truncated. Default is `0.01`.

- maxPS:

  A numerical value (in \`(0, 1\]\`) above which estimated propensity
  scores should be truncated. Must be strictly greater than `minPS`.
  Default is `0.99`.

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

- xvar.smooth.score:

  A vector of characters indicating the name of the variables used as
  the smooth terms if `score.method = 'gam'`. The variables must be
  selected from the variables listed in `cate.model`. Default is `NULL`,
  which uses all variables in `cate.model`.

- xvar.smooth.init:

  A vector of characters indicating the name of the variables used as
  the smooth terms if `initial.predictor.method = 'gam'`. The variables
  must be selected from the variables listed in `init.model`. Default is
  `NULL`, which uses all variables in `init.model`.

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

- error.maxNR:

  A numerical value \> 0 indicating the minimum value of the mean
  absolute error in Newton Raphson algorithm. Used only if
  `score.method = 'contrastReg'`. Default is `0.001`.

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

- ps.control:

  Named options for LASSO propensity-score interaction screening.
  Supported entry: `solver = "C"` (default) or `solver = "R"`. Only
  affects screening when `ps.method = "lasso"`; the LASSO fit and
  interaction selection rule are unchanged. Kept separate from
  boosting's `...` and other numerical controls.

- ...:

  Additional arguments for `gbm()`

## Value

Returns a list containing the following components:

- `ate.gaussian`: A vector of numerical values of length `prop.cutoff`
  containing the estimated ATE in nested subgroups (defined by
  `prop.cutoff`) constructed based on the estimated CATE scores with
  Poisson regression. Only provided if `score.method` includes
  `'gaussian'`.

- `ate.boosting`: Same as `ate.gaussian`, but with the nested subgroups
  based the estimated CATE scores with boosting. Only provided if
  `score.method` includes `'boosting'`.

- `ate.twoReg`: Same as `ate.gaussian`, but with the nested subgroups
  based the estimated CATE scores with two regressions. Only provided if
  `score.method` includes `'twoReg'`.

- `ate.contrastReg`: Same as `ate.gaussian`, but with the nested
  subgroups based the estimated CATE scores with contrast regression.
  Only provided if `score.method` includes `'contrastReg'`.

- `ate.randomForest`: Same as `ate.gaussian`, but with the nested
  subgroups based the estimated CATE scores with random forest. Only
  provided if `score.method` includes `'gam'`.

- `ate.gam`: Same as `ate.gaussian`, but with the nested subgroups based
  the estimated CATE scores with generalized additive model. Only
  provided if `score.method` includes `'gam'`.

- `score.gaussian`: A vector of numerical values of length n (number of
  observations in `data`) containing the estimated CATE scores according
  to the linear regression. Only provided if `score.method` includes
  `'gaussian'`.

- `score.boosting`: Same as `score.gaussian`, but with estimated CATE
  score according to boosting. Only provided if `score.method` includes
  `'boosting'`.

- `score.twoReg`: Same as `score.gaussian`, but with estimated CATE
  score according to two regressions. Only provided if `score.method`
  includes `'twoReg'`.

- `score.contrastReg`: Same as `score.gaussian`, but with estimated CATE
  score according to contrast regression. Only provided if
  `score.method` includes `'contrastReg'`.

- `score.randomForest`: Same as `score.gaussian`, but with estimated
  CATE score according to random forest. Only provided if `score.method`
  includes `'randomForest'`.

- `score.gam`: Same as `score.gaussian`, but with estimated CATE score
  according to generalized additive model. Only provided if
  `score.method` includes `'gam'`.

- `fit`: Additional details on model fitting if `score.method` includes
  'boosting' or 'contrastReg':

  - `result.boosting`: Details on the boosting model fitted to
    observations with treatment = 0 `($fit0.boosting)` and to
    observations with treatment = 1 `($fit1.boosting)`. Only provided if
    `score.method` includes `'boosting'`.

  - `result.randomForest`: Details on the boosting model fitted to
    observations with treatment = 0 `($fit0.randomForest)` and to
    observations with treatment = 1 `($fit1.randomForest)`. Only
    provided if `score.method` includes `'randomForest'`.

  - `result.gam`: Details on the boosting model fitted to observations
    with treatment = 0 `($fit0.gam)` and to observations with treatment
    = 1 `($fit1.gam)`. Only provided if `score.method` includes `'gam'`.

  - `result.contrastReg$sigma.contrastReg`: Variance-covariance matrix
    of the estimated CATE coefficients in contrast regression. Only
    provided if `score.method` includes `'contrastReg'`.

- `coefficients`: A data frame with the coefficients of the estimated
  CATE score by `score.method`. The data frame has number of rows equal
  to the number of covariates in `cate.model` and number of columns
  equal to `length(score.method)`. If `score.method` includes
  `'contrastReg'`, the data frame has an additional column containing
  the standard errors of the coefficients estimated with contrast
  regression. `'boosting'`, `'randomForest'`, `'gam'` do not have
  coefficient results because these methods do not express the CATE as a
  linear combination of coefficients and covariates.

## Details

The CATE score represents an individual-level treatment effect,
estimated with either linear regression, boosting, random forest and
generalized additive model applied separately by treatment group or with
two doubly robust estimators, two regressions and contrast regression
(Yadlowsky, 2020) applied to the entire dataset.

`catefitmean()` provides the coefficients of the CATE score for each
scoring method requested through `score.method`. Currently, contrast
regression is the only method which allows for inference of the CATE
coefficients by providing standard errors of the coefficients. The
coefficients can be used to learn the effect size of each variable and
predict the CATE score for a new observation.

`catefitmean()` also provides the predicted CATE score of each
observation in the data set, for each scoring method. The predictions
allow ranking the observations from potentially high responders to the
treatment to potentially low or standard responders.

The estimated ATE among nested subgroups of high responders are also
provided by scoring method. Note that the ATEs in `catefitmean()` are
derived based on the CATE score which is estimated using the same data
sample. Therefore, overfitting may be an issue. `catefitmean()` is more
suitable to inspect the estimated ATEs across scoring methods as it
implements internal cross validation to reduce optimism.

## References

Yadlowsky, S., Pellegrini, F., Lionetto, F., Braune, S., & Tian, L.
(2020). *Estimation and validation of ratio-based conditional average
treatment effects using observational data. Journal of the American
Statistical Association, 1-18.* DOI: 10.1080/01621459.2020.1772080.

## See also

[`catecvmean()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecvmean.md)
function
