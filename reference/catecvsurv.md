# Cross-validation of the conditional average treatment effect (CATE) score for survival outcomes

Provides doubly robust estimation of the average treatment effect (ATE)
by the RMTL (restricted mean time lost) ratio in nested and mutually
exclusive subgroups of patients defined by an estimated conditional
average treatment effect (CATE) score via cross-validation (CV). The
CATE score can be estimated with up to 5 methods among the following:
Random forest, boosting, poisson regression, two regressions, and
contrast regression (see `score.method`).

## Usage

``` r
catecvsurv(
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
  prop.multi = c(0, 1/3, 2/3, 1),
  abc = TRUE,
  train.prop = 3/4,
  cv.n = 10,
  error.max = 0.1,
  max.iter = 5000,
  surv.min = 0.025,
  tree.depth = 2,
  n.trees.rf = 1000,
  n.trees.boosting = 200,
  B = 3,
  Kfold = 5,
  error.maxNR = 0.001,
  max.iterNR = 150,
  tune = c(0.5, 2),
  seed = NULL,
  plot.gbmperf = TRUE,
  verbose = 0,
  rmst.control = list(),
  contrast.control = list(),
  ps.control = list()
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

- prop.multi:

  A vector of numerical values (in \`\[0, 1\]\`) specifying percentiles
  of the estimated log CATE scores to define mutually exclusive
  subgroups. It should start with 0, end with 1, and be of
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

- plot.gbmperf:

  A logical value indicating whether to plot the performance measures in
  boosting. Used only if `score.method = 'boosting'` or if
  `score.method = 'twoReg'` or `'contrastReg'` and
  `initial.predictor.method = 'boosting'`. Default is `TRUE`.

- verbose:

  An integer value indicating what kind of intermediate progress
  messages should be printed. `0` means no outputs. `1` means only
  progress bar and run time. `2` means progress bar, run time, and all
  errors and warnings. Default is `0`.

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
  affects screening when `ps.method = "lasso"`; the LASSO fit and
  interaction selection rule are unchanged. Kept separate from
  boosting's `...` and other numerical controls.

## Value

Returns a list containing the following components saved as a
`"precmed"` object:

- `ate.randomForest`: A list of ATE output measured by the RMTL ratio if
  `score.method` includes `'randomForest'`:

  - `ate.est.train.high.cv`: A matrix of numerical values with
    `length(prop.cutoff)` rows and `cv.n` columns. The ith column/jth
    row cell contains the estimated ATE in the nested subgroup of high
    responders defined by CATE score above (if `higher.y = FALSE`) or
    below (if `higher.y = TRUE`) the `prop.cutoff[j]`x100% percentile of
    the estimated CATE score in the training set in the ith
    cross-validation iteration.

  - `ate.est.train.low.cv`: A matrix of numerical values with
    `length(prop.cutoff) - 1` rows and `cv.n` columns. TThe ith
    column/jth row cell contains the estimated ATE in the nested
    subgroup of low responders defined by CATE score below (if
    `higher.y = FALSE`) or above (if `higher.y = TRUE`) the
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

- `ate.boosting`: A list of results similar to `ate.randomForest` output
  if `score.method` includes `'boosting'`.

- `ate.poisson`: A list of results similar to `ate.randomForest` output
  if `score.method` includes `'poisson'`.

- `ate.twoReg`: A list of results similar to `ate.randomForest` output
  if `score.method` includes `'twoReg'`.

- `ate.contrastReg`: A list of results similar to `ate.randomForest`
  output if `score.method` includes `'contrastReg'`. This method has an
  additional element in the list of results:

  - `converge.contrastReg.cv`: A vector of logical value of length
    `cv.n`. The ith element indicates whether the algorithm converged in
    the ith cross-validation iteration.

- `hr.randomForest`: A list of adjusted hazard ratio if `score.method`
  includes `'randomForest'`:

  - `hr.est.train.high.cv`: A matrix of numerical values with
    `length(prop.cutoff)` rows and `cv.n` columns. The ith column/jth
    row cell contains the estimated HR in the nested subgroup of high
    responders defined by CATE score above (if `higher.y = FALSE`) or
    below (if `higher.y = TRUE`) the `prop.cutoff[j]`x100% percentile of
    the estimated CATE score in the training set in the ith
    cross-validation iteration.

  - `hr.est.train.low.cv`: A matrix of numerical values with
    `length(prop.cutoff) - 1` rows and `cv.n` columns. TThe ith
    column/jth row cell contains the estimated HR in the nested subgroup
    of low responders defined by CATE score below (if
    `higher.y = FALSE`) or above (if `higher.y = TRUE`) the
    `prop.cutoff[j]`x100% percentile of the estimated CATE score in the
    training set in the ith cross-validation iteration.

  - `hr.est.valid.high.cv`: Same as `hr.est.train.high.cv`, but in the
    validation set.

  - `hr.est.valid.low.cv`: Same as `hr.est.train.low.cv`, but in the
    validation set.

  - `hr.est.train.group.cv`: A matrix of numerical values with
    `length(prop.multi) - 1` rows and `cv.n` columns. The ith column
    contains the estimated HR in `length(prop.multi) - 1` mutually
    exclusive subgroups defined by `prop.multi` in the training set in
    ith cross-validation iteration.

  - `hr.est.valid.group.cv`: Same as `hr.est.train.group.cv`, but in the
    validation set.

- `hr.boosting`: A list of results similar to `hr.randomForest` output
  if `score.method` includes `'boosting'`.

- `hr.poisson`: A list of results similar to `hr.randomForest` output if
  `score.method` includes `'poisson'`.

- `hr.twoReg`: A list of results similar to `hr.randomForest` output if
  `score.method` includes `'twoReg'`.

- `hr.contrastReg`: A list of results similar to `hr.randomForest`
  output if `score.method` includes `'contrastReg'`.

- `props`: A list of 3 elements:

  - `prop.onlyhigh`: The original argument `prop.cutoff`, reformatted as
    necessary.

  - `prop.bi`: The original argument `prop.cutoff`, similar to
    `prop.onlyhigh` but reformatted to exclude 1.

  - `prop.multi`: The original argument `prop.multi`, reformatted as
    necessary to include 0 and 1.

- `overall.ate.train`: A vector of numerical values of length `cv.n`.
  The ith element contains the ATE (RMTL ratio) in the training set of
  the ith cross-validation iteration, estimated with the doubly robust
  estimator.

- `overall.hr.train`: A vector of numerical values of length `cv.n`. The
  ith element contains the ATE (HR) in the training set of the ith
  cross-validation iteration.

- `overall.ate.valid`: A vector of numerical values of length `cv.n`.
  The ith element contains the ATE (RMTL ratio) in the validation set of
  the ith cross-validation iteration, estimated with the doubly robust
  estimator.

- `overall.hr.valid`: A vector of numerical values of length `cv.n`. The
  ith element contains the ATE (HR) in the validation set of the ith
  cross-validation iteration.

- `errors/warnings`: A nested list of errors and warnings that were
  wrapped during the calculation of ATE. Errors and warnings are
  organized by `score.method` and position in the CV flow.

- `higher.y`: The original `higher.y` argument.

- `abc`: The original `abc` argument.

- `cv.n`: The original `cv.n` argument.

- `response`: The type of response. Always 'survival' for this function.

- `formulas`:A list of 3 elements: (1) `cate.model` argument, (2)
  `ps.model` argument and (3) original labels of the left-hand side
  variable in `ps.model` (treatment) if it was not 0/1.

## Details

The CATE score represents an individual-level treatment effect expressed
as the restricted mean survival time (RMTL) ratio) for survival
outcomes. It can be estimated with boosting, Poisson regression, random
forest, and the doubly robust estimator two regressions (Yadlowsky,
2020) applied separately by treatment group or with the other doubly
robust estimator contrast regression (Yadlowsky, 2020) applied to the
entire data set.

Internal CV is applied to reduce optimism in choosing the CATE
estimation method that captures the most treatment effect heterogeneity.
The CV is applied by repeating the following steps `cv.n` times:

1.  Split the data into a training and validation set according to
    `train.prop`. The training and validation sets must be balanced with
    respect to covariate distributions and doubly robust RMTL ratio
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
        `prop.cutoff[i]`x100%) highest (if `higher.y = FALSE`) or lowest
        (if `higher.y = TRUE`) estimated CATE scores.

    2.  Estimate the ATE in the subgroup of high responders using a
        doubly robust estimator.

    3.  Conversely, identify low responders as observations with the 40%
        (i.e., 1 - `prop.cutoff[i]`x100%) lowest (if `higher.y` = FALSE)
        or highest (if `higher.y` = TRUE) estimated CATE scores.

    4.  Estimate the ATE in the subgroup of low responders using a
        doubly robust estimator.

5.  If `abc` = TRUE, calculate the area between the ATE and the series
    of ATEs in nested subgroups of high responders in the validation
    set.

6.  Build mutually exclusive subgroups of treatment responders in the
    training and validation sets, separately, and estimate the ATE
    within each subgroup. Mutually exclusive subgroups are built by
    splitting the estimated CATE scores according to `prop.multi`.

## References

Yadlowsky, S., Pellegrini, F., Lionetto, F., Braune, S., & Tian, L.
(2020). *Estimation and validation of ratio-based conditional average
treatment effects using observational data. Journal of the American
Statistical Association, 1-18.*. DOI: 10.1080/01621459.2020.1772080.

## See also

[`catefitsurv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitsurv.md)
function and [`boxplot()`](https://rdrr.io/r/graphics/boxplot.html),
[`abc`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/abc.md)
methods for `"precmed"` objects.

## Examples

``` r
# \donttest{
library(survival)

tau0 <- with(survivalExample,
             min(quantile(y[trt == "drug1"], 0.95), quantile(y[trt == "drug0"], 0.95)))

catecv <- catecvsurv(data = survivalExample,
                     score.method = "poisson",
                     cate.model = Surv(y, d) ~ age + female + previous_cost +
                                               previous_number_relapses,
                     ps.model = trt ~ age + previous_treatment,
                     initial.predictor.method = "logistic",
                     ipcw.model = ~ age + previous_cost + previous_treatment,
                     tau0 = tau0,
                     higher.y = TRUE,
                     cv.n = 5, seed = 999, verbose = 1)
#> Warning: Variable trt was recoded to 0/1 with drug0->0 and drug1->1.
#>   |                                                                              |                                                                      |   0%
#> cv = 1 
#>   splitting the data..
#>   training..
#>   validating..
#>    Fri Nov 13 16:22:40 2026 
#>   |                                                                              |==============                                                        |  20%
#> cv = 2 
#>   splitting the data..
#>   training..
#>   validating..
#>    Fri Nov 13 16:22:45 2026 
#>   |                                                                              |============================                                          |  40%
#> cv = 3 
#>   splitting the data..
#>   training..
#>   validating..
#>    Fri Nov 13 16:22:49 2026 
#>   |                                                                              |==========================================                            |  60%
#> cv = 4 
#>   splitting the data..
#>   training..
#>   validating..
#>    Fri Nov 13 16:22:53 2026 
#>   |                                                                              |========================================================              |  80%
#> cv = 5 
#>   splitting the data..
#>   training..
#>   validating..
#>    Fri Nov 13 16:22:58 2026 
#>   |                                                                              |======================================================================| 100%
#> Total runtime : 21.91 secs 

# Try:
plot(catecv, ylab = "RMTL ratio of drug1 vs drug0 in each subgroup")

boxplot(catecv, ylab = "RMTL ratio of drug1 vs drug0 in each subgroup")
#> Warning: Removed 4 rows containing non-finite outside the scale range
#> (`stat_boxplot()`).

abc(catecv)
#>               cv1       cv2       cv3       cv4        cv5
#> poisson 0.2165482 0.1988858 0.2184554 0.2111615 0.09635128
# }
```
