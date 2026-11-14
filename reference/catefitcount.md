# Estimation of the conditional average treatment effect (CATE) score for count data

Provides singly robust and doubly robust estimation of CATE score with
up to 5 scoring methods among the following: Poisson regression,
boosting, two regressions, contrast regression, and negative binomial.

## Usage

``` r
catefitcount(
  data,
  score.method,
  cate.model,
  ps.model,
  ps.method = "glm",
  initial.predictor.method = "boosting",
  minPS = 0.01,
  maxPS = 0.99,
  higher.y = TRUE,
  prop.cutoff = seq(0.5, 1, length = 6),
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

  A data frame containing the variables in the outcome and propensity
  score model; a data frame with `n` rows (1 row per observation).

- score.method:

  A vector of one or multiple methods to estimate the CATE score.
  Allowed values are: `'boosting'`, `'poisson'`, `'twoReg'`,
  `'contrastReg'`, and `'negBin'`.

- cate.model:

  A formula describing the outcome model to be fitted. The outcome must
  appear on the left-hand side.

- ps.model:

  A formula describing the propensity score model to be fitted. The
  treatment must appear on the left-hand side. The treatment must be a
  numeric vector coded as 0/1. If data are from a randomized trial,
  specify `ps.model` as an intercept-only model.

- ps.method:

  A character value for the method to estimate the propensity score.
  Allowed values include one of: `'glm'` for logistic regression with
  main effects only (default), or `'lasso'` for a logistic regression
  with main effects and LASSO penalization on two-way interactions
  (added to the model if interactions are not specified in `ps.model`).
  Relevant only when `ps.model` has more than one variable.

- initial.predictor.method:

  A character vector for the method used to get initial outcome
  predictions conditional on the covariates in `cate.model`. Only
  applies when `score.method` includes `'twoReg'` or `'contrastReg'`.
  Allowed values include one of `'poisson'` (fastest), `'boosting'` and
  `'gam'`. Default is `'boosting'`.

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

- xvar.smooth:

  A vector of characters indicating the name of the variables used as
  the smooth terms if `initial.predictor.method = 'gam'`. The variables
  must be selected from the variables listed in `cate.model`. Default is
  `NULL`, which uses all variables in `cate.model`.

- tree.depth:

  A positive integer specifying the depth of individual trees in
  boosting (usually 2-3). Used only if `score.method = 'boosting'` or if
  `initial.predictor.method = 'boosting'` with `score.method = 'twoReg'`
  or `'contrastReg'`. Default is 2.

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
  progress and run time. `2` means progress, run time, and all errors
  and warnings. Default is `0`.

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
  affects screening when `ps.method = "lasso"`; the LASSO fit and
  interaction selection rule are unchanged. Kept separate from
  boosting's `...` and other numerical controls.

- ...:

  Additional boosting arguments; see
  [`boosting-control`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/boosting-control.md).

## Value

Returns a list containing the following components:

- `ate.poisson`: A vector of numerical values of length `prop.cutoff`
  containing the estimated ATE in nested subgroups (defined by
  `prop.cutoff`) constructed based on the estimated CATE scores with
  poisson regression. Only provided if `score.method` includes
  `'poisson'`.

- `ate.boosting`: Same as `ate.poisson`, but with the nested subgroups
  based the estimated CATE scores with boosting. Only provided if
  `score.method` includes `'boosting'`.

- `ate.twoReg`: Same as `ate.poisson`, but with the nested subgroups
  based the estimated CATE scores with two regressions. Only provided if
  `score.method` includes `'twoReg'`.

- `ate.contrastReg`: Same as `ate.poisson`, but with the nested
  subgroups based the estimated CATE scores with contrast regression.
  Only provided if `score.method` includes `'contrastReg'`.

- `ate.negBin`: Same as `ate.poisson`, but with the nested subgroups
  based the estimated CATE scores with negative binomial regression.
  Only provided if `score.method` includes `'negBin'`.

- `score.poisson`: A vector of numerical values of length n (number of
  observations in `data`) containing the estimated log-CATE scores
  according to the Poisson regression. Only provided if `score.method`
  includes `'poisson'`.

- `score.boosting`: Same as `score.poisson`, but with estimated log-CATE
  score according to boosting. Only provided if `score.method` includes
  `'boosting'`.

- `score.twoReg`: Same as `score.poisson`, but with estimated log-CATE
  score according to two regressions. Only provided if `score.method`
  includes `'twoReg'`.

- `score.contrastReg`: Same as `score.poisson`, but with estimated
  log-CATE score according to contrast regression. Only provided if
  `score.method` includes `'contrastReg'`.

- `score.negBin`: Same as `score.poisson`, but with estimated log-CATE
  score according to negative binomial regression. Only provided if
  `score.method` includes `'negBin'`.

- `fit`: Additional details on model fitting if `score.method` includes
  'boosting' or 'contrastReg':

  - `result.boosting`: Details on the boosting model fitted to
    observations with treatment = 0 `($fit0.boosting)` and to
    observations with treatment = 1 `($fit1.boosting)`. Only provided if
    `score.method` includes `'boosting'`.

  - `result.contrastReg$sigma.contrastReg`: Variance-covariance matrix
    of the estimated log-CATE coefficients in contrast regression. Only
    provided if `score.method` includes `'contrastReg'`.

- `coefficients`: A data frame with the coefficients of the estimated
  log-CATE score by `score.method`. The data frame has number of rows
  equal to the number of covariates in `cate.model` and number of
  columns equal to `length(score.method)`. If `score.method` includes
  `'contrastReg'`, the data frame has an additional column containing
  the standard errors of the coefficients estimated with contrast
  regression. `'boosting'` does not have coefficient results because
  tree-based methods typically do not express the log-CATE as a linear
  combination of coefficients and covariates.

## Details

The CATE score represents an individual-level treatment effect,
estimated with either Poisson regression, boosting or negative binomial
regression applied separately by treatment group or with two doubly
robust estimators, two regressions and contrast regression (Yadlowsky,
2020) applied to the entire dataset.

`catefitcount()` provides the coefficients of the CATE score for each
scoring method requested through `score.method`. Currently, contrast
regression is the only method which allows for inference of the CATE
coefficients by providing standard errors of the coefficients. The
coefficients can be used to learn the effect size of each variable and
predict the CATE score for a new observation.

`catefitcount()` also provides the predicted CATE score of each
observation in the data set, for each scoring method. The predictions
allow ranking the observations from potentially high responders to the
treatment to potentially low or standard responders.

The estimated ATE among nested subgroups of high responders are also
provided by scoring method. Note that the ATEs in `catefitcount()` are
derived based on the CATE score which is estimated using the same data
sample. Therefore, overfitting may be an issue.
[`catecvcount()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecvcount.md)
is more suitable to inspect the estimated ATEs across scoring methods as
it implements internal cross validation to reduce optimism.

## References

Yadlowsky, S., Pellegrini, F., Lionetto, F., Braune, S., & Tian, L.
(2020). *Estimation and validation of ratio-based conditional average
treatment effects using observational data. Journal of the American
Statistical Association, 1-18.* DOI: 10.1080/01621459.2020.1772080.

## See also

[`catecvcount()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecvcount.md)

## Examples

``` r
fit <- catefitcount(data = countExample,
                    score.method = "poisson",
                    cate.model = y ~ age + female + previous_treatment +
                                 previous_cost + previous_number_relapses +
                                 offset(log(years)),
                    ps.model = trt ~ age + previous_treatment,
                    higher.y = FALSE,
                    seed = 999, verbose = 1)
#> Warning: Variable trt was recoded to 0/1 with drug0->0 and drug1->1.
#> Total runtime : 0.16 secs 
```
