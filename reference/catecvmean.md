# Cross-validation of the conditional average treatment effect (CATE) score for continuous outcomes

Provides doubly robust estimation of the average treatment effect (ATE)
in nested and mutually exclusive subgroups of patients defined by an
estimated conditional average treatment effect (CATE) score via
cross-validation (CV). The CATE score can be estimated with up to 6
methods among the following: Linear regression, boosting, two
regressions, contrast regression, random forest and generalized additive
model (see `score.method`).

## Usage

``` r
catecvmean(
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
  prop.multi = c(0, 1/3, 2/3, 1),
  abc = TRUE,
  train.prop = 3/4,
  cv.n = 10,
  error.max = 0.1,
  max.iter = 5000,
  xvar.smooth.score = NULL,
  xvar.smooth.init = NULL,
  tree.depth = 2,
  n.trees.rf = 1000,
  n.trees.boosting = 200,
  B = 3,
  Kfold = 6,
  plot.gbmperf = TRUE,
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
  (added to the model if interactions are not specified in `ps.model`).
  Relevant only when `ps.model` has more than one variable.

- init.model:

  A formula describing the initial predictor model. The outcome must
  appear on the left-hand side. It must be specified when
  `score.method = contrastReg` or `twoReg`.

- initial.predictor.method:

  A character vector for the method used to get initial outcome
  predictions conditional on the covariates in `cate.model` in
  `score.method = 'twoReg'` and `'contrastReg'`. Allowed values include
  one of `'poisson'` (fastest), `'boosting'` and `'gam'`. Default is
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
  of the estimated CATE scores to define nested subgroups. Each element
  represents the cutoff to separate observations in nested subgroups
  (below vs above cutoff). The length of `prop.cutoff` is the number of
  nested subgroups. An equally-spaced sequence of proportions ending
  with 1 is recommended. Default is `seq(0.5, 1, length = 6)`.

- prop.multi:

  A vector of numerical values (in \`\[0, 1\]\`) specifying percentiles
  of the estimated CATE scores to define mutually exclusive subgroups.
  It should start with 0, end with 1, and be of
  `length(prop.multi) > 2`. Each element represents the cutoff to
  separate the observations into `length(prop.multi) - 1` mutually
  exclusive subgroups. Default is `c(0, 1/3, 2/3, 1)`.

- abc:

  A logical value indicating whether the area between curves (ABC)
  should be calculated at each cross-validation iterations, for each
  `score.method`. Default is `TRUE`.

- train.prop:

  A numerical value (in \`(0, 1)\`) indicating the proportion of total
  data used for training. Default is `3/4`.

- cv.n:

  A positive integer value indicating the number of cross-validation
  iterations. Default is `10`.

- error.max:

  A numerical value \> 0 indicating the tolerance (maximum value of
  error) for the largest standardized absolute difference in the
  covariate distributions or in the doubly robust estimated rate ratios
  between the training and validation sets. This is used to define a
  balanced training-validation splitting. Default is `0.1`.

- max.iter:

  A positive integer value indicating the maximum number of iterations
  when searching for a balanced training-validation split. Default is
  `5,000`.

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

  A positive integer specifying the maximum number of trees in random
  forest. Used if `score.method = 'ranfomForest'` or if
  `initial.predictor.method = 'randomForest'` with
  `score.method = 'twoReg'` or `'contrastReg'`. Only applies for
  survival outcomes. Default is `1000`.

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
  progress bar and run time. `2` means progress bar, run time, and all
  errors and warnings. Default is `0`.

- ps.control:

  Named options for LASSO propensity-score interaction screening.
  Supported entry: `solver = "C"` (default) or `solver = "R"`. Only
  affects screening when `ps.method = "lasso"`; the LASSO fit and
  interaction selection rule are unchanged. Kept separate from
  boosting's `...` and other numerical controls.

- ...:

  Additional arguments for `gbm()`

## Value

Returns a list containing the following components saved as a
`"precmed"` object:

- `ate.gaussian`: A list of results output if `score.method` includes
  `'gaussian'`:

  - `ate.est.train.high.cv`: A matrix of numerical values with
    `length(prop.cutoff)` rows and `cv.n` columns. The ith column/jth
    row cell contains the estimated ATE in the nested subgroup of high
    responders defined by CATE score above (if `higher.y = TRUE`) or
    below (if `higher.y = FALSE`) the `prop.cutoff[j]`x100% percentile
    of the estimated CATE score in the training set in the ith
    cross-validation iteration.

  - `ate.est.train.low.cv`: A matrix of numerical values with
    `length(prop.cutoff) - 1` rows and `cv.n` columns. The ith
    column/jth row cell contains the estimated ATE in the nested
    subgroup of low responders defined by CATE score below (if
    `higher.y = TRUE`) or above (if `higher.y = FALSE`) the
    `prop.cutoff[j]`x100% percentile of the estimated CATE score in the
    training set in the ith cross-validation iteration.

  - `ate.est.valid.high.cv`: Same as `ate.est.train.high.cv`, but in the
    validation set.

  - `ate.est.valid.low.cv`: Same as `ate.est.train.low.cv`, but in the
    validation set.

  - `ate.est.train.group.cv`: A matrix of numerical values with
    `length(prop.multi) - 1` rows and `cv.n` columns. The ith column
    contains the estimated ATE in `length(prop.multi) - 1` mutually
    exclusive subgroups defined by `prop.multi` in the training set in
    ith cross-validation iteration.

  - `ate.est.valid.group.cv`: Same as `ate.est.train.group.cv`, but in
    the validation set.

  - `abc.valid`: A vector of numerical values of length `cv.n`, The ith
    element returns the ABC of the validation curve in the ith
    cross-validation iteration. Only returned if `abc = TRUE`.

- `ate.boosting`: A list of results similar to `ate.gaussian` output if
  `score.method` includes `'boosting'`.

- `ate.twoReg`: A list of results similar to `ate.gaussian` output if
  `score.method` includes `'twoReg'`.

- `ate.contrastReg`: A list of results similar to `ate.gaussian` output
  if `score.method` includes `'contrastReg'`.

- `ate.randomForest`: A list of ATE output measured by the RMTL ratio if
  `score.method` includes `'randomForest'`:

- `ate.gam`: A list of results similar to `ate.gaussian` output if
  `score.method` includes `'gam'`.

- `props`: A list of 3 elements:

  - `prop.onlyhigh`: The original argument `prop.cutoff`, reformatted as
    necessary.

  - `prop.bi`: The original argument `prop.cutoff`, similar to
    `prop.onlyhigh` but reformatted to exclude 1.

  - `prop.multi`: The original argument `prop.multi`, reformatted as
    necessary.

- `overall.ate.train`: A vector of numerical values of length `cv.n`.
  The ith element contains the ATE in the training set of the ith
  cross-validation iteration, estimated with the doubly robust
  estimator.

- `overall.ate.valid`: A vector of numerical values of length `cv.n`.
  The ith element contains the ATE in the validation set of the ith
  cross-validation iteration, estimated with the doubly robust
  estimator.

- `higher.y`: The original `higher.y` argument.

- `abc`: The original `abc` argument.

- `cv.n`: The original `cv.n` argument.

- `response`: The type of response. Always 'continuous' for this
  function.

- `formulas`:A list of 3 elements: (1) `cate.model` argument, (2)
  `ps.model` argument and (3) original labels of the left-hand side
  variable in `ps.model` (treatment) if it was not 0/1.

## Details

The CATE score represents an individual-level treatment effect for
continuous data, estimated with boosting, linear regression, random
forest, generalized additive model and the doubly robust estimator (two
regressions, Yadlowsky, 2020) applied separately by treatment group or
with the other doubly robust estimators (contrast regression, Yadlowsky,
2020) applied to the entire data set.

Internal CV is applied to reduce optimism in choosing the CATE
estimation method that captures the most treatment effect heterogeneity.
The CV is applied by repeating the following steps `cv.n` times:

1.  Split the data into a training and validation set according to
    `train.prop`. The training and validation sets must be balanced with
    respect to covariate distributions and doubly robust rate ratio
    estimates (see `error.max`).

2.  Estimate the CATE score in the training set with the specified
    scoring method.

3.  Predict the CATE score in the validation set using the scoring model
    fitted from the training set.

4.  Build nested subgroups of treatment responders in the training and
    validation sets, separately, and estimate the ATE within each nested
    subgroup. For each element i of `prop.cutoff` (e.g.,
    `prop.cutoff[i]` = 0.6), take the following steps:

    1.  Identify high responders as observations with the 60% (i.e.,
        `prop.cutoff[i]`x100%) highest (if `higher.y = TRUE`) or lowest
        (if `higher.y = FALSE`) estimated CATE scores.

    2.  Estimate the ATE in the subgroup of high responders using a
        doubly robust estimator.

    3.  Conversely, identify low responders as observations with the 40%
        (i.e., 1 - `prop.cutoff[i]`x100%) lowest (if `higher.y` = TRUE)
        or highest (if `higher.y` = FALSE) estimated CATE scores.

    4.  Estimate the ATE in the subgroup of low responders using a
        doubly robust estimator.

5.  Build mutually exclusive subgroups of treatment responders in the
    training and validation sets, separately, and estimate the ATE
    within each subgroup. Mutually exclusive subgroups are built by
    splitting the estimated CATE scores according to `prop.multi`.

6.  If `abc` = TRUE, calculate the area between the ATE and the series
    of ATEs in nested subgroups of high responders in the validation
    set.

## References

Yadlowsky, S., Pellegrini, F., Lionetto, F., Braune, S., & Tian, L.
(2020). *Estimation and validation of ratio-based conditional average
treatment effects using observational data. Journal of the American
Statistical Association, 1-18.* DOI: 10.1080/01621459.2020.1772080.

## See also

[`plot.precmed()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/plot.precmed.md),
[`boxplot.precmed()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/boxplot.precmed.md),
[`abc()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/abc.md)
methods for `"precmed"` objects, and
[`catefitmean()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitmean.md)
function.

## Examples

``` r
# Not implemented yet!
```
