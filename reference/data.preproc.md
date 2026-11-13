# Data preprocessing Apply at the beginning of `pmcount()` and `cvcount()`, after `arg.checks()`

Data preprocessing Apply at the beginning of `pmcount()` and
`cvcount()`, after
[`arg.checks()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/arg.checks.md)

## Usage

``` r
data.preproc(
  fun,
  cate.model,
  ps.model,
  data,
  prop.cutoff = NULL,
  prop.multi = NULL,
  ps.method,
  initial.predictor.method = NULL
)
```

## Arguments

- fun:

  A function for which argument check is needed; "pm" for `pmcount()`,
  "cv" for `cvcount()`, and "drinf" for `drcount.inference()`. No
  default.

- cate.model:

  A formula describing the outcome model to be fitted. The outcome must
  appear on the left-hand side.

- ps.model:

  A formula describing the propensity score model to be fitted. The
  treatment must appear on the left-hand side. The treatment must be a
  numeric vector coded as 0/1. If data are from a RCT, specify
  `ps.model` as an intercept-only model.

- data:

  A data frame containing the variables in the outcome and propensity
  score models; a data frame with `n` rows (1 row per observation).

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

- ps.method:

  A character value for the method to estimate the propensity score.
  Allowed values include one of: `'glm'` for logistic regression with
  main effects only (default), or `'lasso'` for a logistic regression
  with main effects and LASSO penalization on two-way interactions
  (added to the model if interactions are not specified in `ps.model`).
  Relevant only when `ps.model` has more than one variable.

- initial.predictor.method:

  A character vector for the method used to get initial outcome
  predictions conditional on the covariates. Only applies when
  `score.method` includes `'twoReg'` or `'contrastReg'`. Allowed values
  include one of `'boosting'`, `'poisson'` (fast), and `'gam'`. Default
  is `NULL`, which assigns `'boosting'` for count outcomes.

## Value

A list of 6 elements: - y: outcome; vector of length `n`
(observations) - trt: binary treatment; vector of length `n` - x.ps:
matrix of `p.ps` baseline covariates (plus intercept); dimension `n` by
`p.ps + 1` - x.cate: matrix of `p.cate` baseline covariates; dimension
`n` by `p.cate` - time: offset; vector of length `n` - if
`fun = "pm"`: - prop: formatted `prop.cutoff` - if `fun = "cv"` -
prop.onlyhigh: formatted `prop.cutoff` with 0 removed if applicable -
prop.bi; formatted `prop.cutoff` with 0 and 1 removed if applicable -
prop.multi: formatted `prop.multi`, starting with 0 and ending with 1
