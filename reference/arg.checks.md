# Check arguments Catered to all types of outcome Apply at the beginning of `pmcount()`, `cvcount()`, `drcount.inference()`, `catefitsurv()`, `catecvsurv()`, and `drsurv.inference()`

Check arguments Catered to all types of outcome Apply at the beginning
of `pmcount()`, `cvcount()`, `drcount.inference()`,
[`catefitsurv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitsurv.md),
[`catecvsurv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecvsurv.md),
and `drsurv.inference()`

## Usage

``` r
arg.checks(
  fun,
  response,
  data,
  followup.time = NULL,
  tau0 = NULL,
  surv.min = NULL,
  ipcw.method = NULL,
  ps.method,
  minPS,
  maxPS,
  higher.y = NULL,
  score.method = NULL,
  abc = NULL,
  prop.cutoff = NULL,
  prop.multi = NULL,
  train.prop = NULL,
  cv.n = NULL,
  error.max = NULL,
  max.iter = NULL,
  initial.predictor.method = NULL,
  tree.depth = NULL,
  n.trees.rf = NULL,
  n.trees.boosting = NULL,
  B = NULL,
  Kfold = NULL,
  plot.gbmperf = NULL,
  error.maxNR = NULL,
  max.iterNR = NULL,
  tune = NULL,
  n.boot = NULL,
  plot.boot = NULL,
  interactions = NULL
)
```

## Arguments

- fun:

  A function for which argument check is needed; "catefit" for
  [`catefitcount()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitcount.md)
  and
  [`catefitsurv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitsurv.md),
  "crossv" for
  [`catecvcount()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecvcount.md)
  and
  [`catecvsurv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecvsurv.md),
  and "drinf" for `drcount.inference()` and `drsurv.inference()`. No
  default.

- response:

  The type of response. Always 'survival' for this function.

- data:

  A data frame containing the variables in the outcome and propensity
  score models; a data frame with `n` rows (1 row per observation).

- followup.time:

  Follow-up time, interpreted as the potential censoring time. If the
  potential censoring time is known, followup.time is the name of a
  corresponding column in the data. Otherwise, set
  `followup.time == NULL`.

- tau0:

  The truncation time for defining restricted mean time lost.

- surv.min:

  Lower truncation limit for probability of being censored (positive and
  very close to 0).

- ipcw.method:

  The censoring model. Allowed values are: `'breslow'` (Cox regression
  with Breslow estimator of the baseline survivor function),
  `'aft (exponential)'`, `'aft (weibull)'`, `'aft (lognormal)'` or
  `'aft (loglogistic)'`. Default is `'breslow'`.

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

  A numerical value (in \`(0, 1\]\`) above which estimated propensity
  scores should be truncated. Must be strictly greater than `minPS`.
  Default is `0.99`.

- higher.y:

  A logical value indicating whether higher (`TRUE`) or lower (`FALSE`)
  values of the outcome are more desirable. Default is `TRUE`.

- score.method:

  A vector of one or multiple methods to estimate the CATE score.
  Allowed values are: `'boosting'`, `'poisson'`, `'twoReg'`,
  `'contrastReg'`, `'negBin'`. Default specifies all 5 methods.

- abc:

  A logical value indicating whether the area between curves (ABC)
  should be calculated at each cross-validation iterations, for each
  `score.method`. Default is `TRUE`.

- prop.cutoff:

  A vector of numerical values (in \`(0, 1\]\`) specifying percentiles
  of the estimated log CATE scores to define nested subgroups. Each
  element represents the cutoff to separate observations in nested
  subgroups (below vs above cutoff). The length of `prop.cutoff` is the
  number of nested subgroups. An equally-spaced sequence of proportions
  ending with 1 is recommended. Default is `seq(0.5, 1, length = 6)`.

- prop.multi:

  A vector of numerical values (in \`\[0, 1\]\`) specifying percentiles
  of the estimated log CATE scores to define mutually exclusive
  subgroups. It should start with 0, end with 1, and be of
  `length(prop.multi) > 2`. Each element represents the cutoff to
  separate the observations into `length(prop.multi) - 1` mutually
  exclusive subgroups. Default is `c(0, 1/3, 2/3, 1)`.

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
  `initial.predictor.method = 'boosting'`. Default is `2`.

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

- n.boot:

  A numeric value indicating the number of bootstrap samples used. This
  is only relevant if `inference = TRUE`. Default is `500`.

- plot.boot:

  A logic value indicating whether histograms of the bootstrapped
  log(rate ratio) should be produced at every `n.boot/10`-th iteration
  and whether the final histogram should be outputted. Default is
  `FALSE`.

- interactions:

  A logical value indicating whether the outcome model should assume
  interactions between `x` and `trt`. If `TRUE`, interactions will be
  assumed only if at least 10 patients received each treatment option.
  Default is `TRUE`.

## Value

Nothing. Will stop if arguments are incorrect.
