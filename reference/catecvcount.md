# Cross-validation of the conditional average treatment effect (CATE) score for count outcomes

Provides doubly robust estimation of the average treatment effect (ATE)
in nested and mutually exclusive subgroups of patients defined by an
estimated conditional average treatment effect (CATE) score via
cross-validation (CV). The CATE score can be estimated with up to 5
methods among the following: Poisson regression, boosting, two
regressions, contrast regression, and negative binomial (see
`outcome_method`).

## Usage

``` r
catecvcount(
  data,
  outcome_model,
  outcome_method = "poisson",
  benefit_model = NULL,
  benefit_method = "t_learner",
  ps_model,
  ps_method = "glm",
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
  xvar.smooth = NULL,
  tree.depth = 2,
  n.trees.boosting = 200,
  B = 3,
  Kfold = 5,
  error.maxNR = 0.001,
  max.iterNR = 150,
  tune = c(0.5, 2),
  seed = NULL,
  plot.gbmperf = FALSE,
  verbose = 0,
  contrast.control = list(),
  ps.control = list(),
  ...
)
```

## Arguments

- data:

  A data frame containing the variables in the outcome and treatment
  model; a data frame with `n` rows (1 row per observation).

- outcome_model:

  A formula specifying the outcome and baseline covariates. Use a `Surv`
  response for survival and include an exposure offset for count
  outcomes.

- outcome_method:

  Outcome-model method for each `benefit_method`. For a T-learner,
  choose `'poisson'` (default), `'negative_binomial'`, or `'boosting'`.
  For two-regression or contrast regression, choose `'poisson'`,
  `'gam'`, or `'boosting'`.

- benefit_model:

  Optional one-sided formula for the covariates and terms used in
  two-regression and contrast-regression scores. The default uses the
  terms in `outcome_model`. T-learners use `outcome_model` directly.

- benefit_method:

  Score-construction method, with one entry per `outcome_method`. The
  default is `'t_learner'`. See
  [`catefit`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
  for supported combinations.

- ps_model:

  A formula for the treatment model, with treatment on the left. It
  estimates propensity scores. For a common randomization probability,
  use an intercept-only formula such as `trt ~ 1`.

- ps_method:

  Treatment-model fitting method. The default `'glm'` uses logistic
  regression. `'lasso'` uses a logistic model with LASSO screening for
  two-way interactions.

- minPS:

  A numerical value between 0 and 1 below which estimated propensity
  scores should be truncated. Default is `0.01`.

- maxPS:

  A numerical value between 0 and 1 above which estimated propensity
  scores should be truncated. Must be strictly greater than `minPS`.
  Default is `0.99`.

- higher.y:

  A logical value indicating whether higher (`TRUE`) or lower (`FALSE`)
  values of the outcome are more desirable. Default is `TRUE`.

- prop.cutoff:

  A vector of numerical values between 0 and 1 specifying percentiles of
  the estimated log CATE scores to define nested subgroups. Each element
  represents the cutoff to separate observations in nested subgroups
  (below vs above cutoff). The length of `prop.cutoff` is the number of
  nested subgroups. An equally-spaced sequence of proportions ending
  with 1 is recommended. Default is `seq(0.5, 1, length = 6)`.

- prop.multi:

  A vector of numerical values between 0 and 1 specifying percentiles of
  the estimated log CATE scores to define mutually exclusive subgroups.
  It should start with 0, end with 1, and be of
  `length(prop.multi) > 2`. Each element represents the cutoff to
  separate the observations into `length(prop.multi) - 1` mutually
  exclusive subgroups. Default is `c(0, 1/3, 2/3, 1)`.

- abc:

  A logical value indicating whether the area between curves (ABC)
  should be calculated at each cross-validation iterations, for each
  scoring method. Default is `TRUE`.

- train.prop:

  A numerical value between 0 and 1 indicating the proportion of total
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

- xvar.smooth:

  Names of outcome-model variables to use as smooth terms with
  `outcome_method = 'gam'`. Defaults to all eligible variables.

- tree.depth:

  Maximum depth of boosting trees. Defaults to 2.

- n.trees.boosting:

  Maximum number of trees for boosting. Defaults to 200.

- B:

  Number of cross-fitting repetitions used by two-regression and
  contrast regression. Defaults to 3.

- Kfold:

  Number of folds used for cross-fitting within each training set.
  Defaults to 5.

- error.maxNR:

  Newton–Raphson convergence tolerance for contrast regression. Defaults
  to 0.001.

- max.iterNR:

  Maximum number of Newton–Raphson iterations for contrast regression.
  Defaults to 150.

- tune:

  Two positive Newton–Raphson tuning values for contrast regression.
  Defaults to `c(0.5, 2)`.

- seed:

  An optional integer specifying an initial randomization seed for
  reproducibility. Default is `NULL`, corresponding to no seed.

- plot.gbmperf:

  Whether to plot boosting performance. Defaults to FALSE.

- verbose:

  An integer value indicating what kind of intermediate progress
  messages should be printed. `0` means no outputs. `1` means only
  progress bar and run time. `2` means progress bar, run time, and all
  errors and warnings. Default is `0`.

- contrast.control:

  A named list of options for count or survival contrast regression.
  Supported names are `solver` (`"C"`, the default, or `"R"`),
  `error.maxNR`, `max.iterNR`, `tune`, and `optim.control` (a named list
  passed to [`stats::optim`](https://rdrr.io/r/stats/optim.html) on
  fallback). List entries override the corresponding top-level NR
  arguments. Options are kept separate from `...`, which is forwarded to
  boosting. Count and survival outcomes support this list; continuous
  outcomes retain their existing solver.

- ps.control:

  Named options for LASSO propensity-score interaction screening.
  Supported entry: `solver = "C"` (default) or `solver = "R"`. Only
  affects screening when `ps_method = "lasso"`; the LASSO fit and
  interaction selection rule are unchanged. Kept separate from
  boosting's `...` and other numerical controls.

- ...:

  Additional arguments. Deprecated argument names from precmed 1.1.0
  remain supported for backward compatibility and emit warnings.

## Value

A `"precmed"` object with the components below. When doubly robust
scores use different outcome methods, `outcome_methods` records the
method used for each score.

- `ate.poisson`: A list of results output if the selected scoring
  methods includes `'poisson'`:

  - `ate.est.train.high.cv`: A matrix of numerical values with
    `length(prop.cutoff)` rows and `cv.n` columns. The ith row/jth
    column cell contains the estimated ATE in the nested subgroup of
    high responders defined by CATE score above (if `higher.y = TRUE`)
    or below (if `higher.y = FALSE`) the `prop.cutoff[i]`x100%
    percentile of the estimated CATE score in the training set in the
    jth cross-validation iteration.

  - `ate.est.train.low.cv`: A matrix of numerical values with
    `length(prop.cutoff) - 1` rows and `cv.n` columns. The ith row/jth
    column cell contains the estimated ATE in the nested subgroup of low
    responders defined by CATE score below (if `higher.y = TRUE`) or
    above (if `higher.y = FALSE`) the `prop.cutoff[i]`x100% percentile
    of the estimated CATE score in the training set in the jth
    cross-validation iteration.

  - `ate.est.valid.high.cv`: Same as `ate.est.train.high.cv`, but in the
    validation set.

  - `ate.est.valid.low.cv`: Same as `ate.est.train.low.cv`, but in the
    validation set.

  - `ate.est.train.group.cv`: A matrix of numerical values with
    `length(prop.multi) - 1` rows and `cv.n` columns. The jth column
    contains the estimated ATE in `length(prop.multi) - 1` mutually
    exclusive subgroups defined by `prop.multi` in the training set in
    jth cross-validation iteration.

  - `ate.est.valid.group.cv`: Same as `ate.est.train.group.cv`, but in
    the validation set.

  - `abc.valid`: A vector of numerical values of length `cv.n`. The ith
    element returns the ABC of the validation curve in the ith
    cross-validation iteration. Only returned if `abc = TRUE`.

- `ate.boosting`: A list of results similar to `ate.poisson` output if
  the selected scoring methods includes `'boosting'`.

- `ate.twoReg`: A list of results similar to `ate.poisson` output if the
  selected scoring methods includes `'twoReg'`.

- `ate.contrastReg`: A list of results similar to `ate.poisson` output
  if the selected scoring methods includes `'contrastReg'`. This method
  has additional convergence information:

  - `converge.contrastReg.cv`: A vector of logical value of length
    `cv.n`. The ith element is true if at least one cross-fitting
    repetition converged in the ith cross-validation iteration.

  - `converge.repetitions.cv`: A list of logical vectors recording
    convergence of each cross-fitting repetition in each CV iteration.

- `ate.negBin`: A list of results similar to `ate.poisson` output if the
  selected scoring methods includes `'negBin'`.

- `props`: A list of 3 elements:

  - `prop.onlyhigh`: The original argument `prop.cutoff`, reformatted as
    necessary.

  - `prop.bi`: The original argument `prop.cutoff`, similar to
    `prop.onlyhigh` but reformatted to exclude 1.

  - `prop.multi`: The original argument `prop.multi`, reformatted as
    necessary to include 0 and 1.

- `overall.ate.valid`: A vector of numerical values of length `cv.n`.
  The ith element contains the ATE in the validation set of the ith
  cross-validation iteration, estimated with the doubly robust
  estimator.

- `overall.ate.train`: A vector of numerical values of length `cv.n`.
  The ith element contains the ATE in the training set of the ith
  cross-validation iteration, estimated with the doubly robust
  estimator.

- `fgam`: The formula used in GAM if `outcome_method = 'gam'`.

- `higher.y`: The original `higher.y` argument.

- `abc`: The original `abc` argument.

- `cv.n`: The original `cv.n` argument.

- `response`: The type of response. Always 'count' for this function.

- `formulas`:A list of 3 elements: (1) `outcome_model` argument, (2)
  `ps_model` argument and (3) original labels of the left-hand side
  variable in `ps_model` (treatment) if it was not 0/1.

## Details

The CATE score represents an individual-level treatment effect expressed
as a rate ratio for count outcomes. It can be estimated with boosting,
Poisson regression, negative binomial regression, and the doubly robust
estimator two regressions (Yadlowsky, 2020) applied separately by
treatment group or with the other doubly robust estimator contrast
regression (Yadlowsky, 2020) applied to the entire data set.

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

[`plot.precmed()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/plot.precmed.md),
[`boxplot.precmed()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/boxplot.precmed.md),
[`abc()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/abc.md)
methods for `"precmed"` objects, and
[`catefitcount()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitcount.md)
function.

## Examples

``` r
# \donttest{
catecv <- catecvcount(data = countExample,
                      outcome_method = "poisson",
                      outcome_model = y ~ age + female + previous_treatment +
                                   previous_cost + previous_number_relapses +
                                   offset(log(years)),
                      ps_model = trt ~ age + previous_treatment,
                      higher.y = FALSE,
                      cv.n = 5,
                      seed = 999,
                      plot.gbmperf = FALSE,
                      verbose = 1)
#> Variable trt was recoded to 0/1 with drug0->0 and drug1->1.
#>   |                                                                              |                                                                      |   0%
#> cv = 1 
#>   splitting the data..
#>   training..
#>   validating..
#> 
#> cv = 2 
#>   splitting the data..
#>   training..
#>   validating..
#> 
#> cv = 3 
#>   splitting the data..
#>   training..
#>   validating..
#> 
#> cv = 4 
#>   splitting the data..
#>   training..
#>   validating..
#> 
#> cv = 5 
#>   splitting the data..
#>   training..
#>   validating..
#> 
#> Total runtime : 6.21 secs 

plot(catecv, ylab = "Rate ratio of drug1 vs drug0 in each subgroup")

boxplot(catecv, ylab = "Rate ratio of drug1 vs drug0 in each subgroup")

abc(catecv)
#>               cv1       cv2       cv3        cv4       cv5
#> poisson 0.1375899 0.1172791 0.1120313 0.06136285 0.1877103
# }
```
