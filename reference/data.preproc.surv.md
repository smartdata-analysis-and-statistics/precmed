# Data preprocessing Apply at the beginning of `catefitcount()`, `catecvcount()`, `catefitsurv()`, and `catecvsurv()`, after `arg.checks()`

Data preprocessing Apply at the beginning of
[`catefitcount()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitcount.md),
[`catecvcount()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecvcount.md),
[`catefitsurv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitsurv.md),
and
[`catecvsurv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecvsurv.md),
after
[`arg.checks()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/arg.checks.md)

## Usage

``` r
data.preproc.surv(
  fun,
  cate.model,
  ps.model,
  ipcw.model = NULL,
  tau0 = NULL,
  data,
  prop.cutoff = NULL,
  prop.multi = NULL,
  ps.method,
  initial.predictor.method = NULL,
  response = "count"
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

- cate.model:

  A formula describing the outcome model to be fitted. The outcome must
  appear on the left-hand side.

- ps.model:

  A formula describing the propensity score model to be fitted. The
  treatment must appear on the left-hand side. The treatment must be a
  numeric vector coded as 0/1. If data are from a RCT, specify
  `ps.model` as an intercept-only model.

- ipcw.model:

  A formula describing inverse probability of censoring weighting(IPCW)
  model to be fitted. If covariates are the same as outcome model, set
  `ipcw.model = NULL`. Otherwise, the left-hand side must be empty and
  the right-hand side is a covariates model.

- tau0:

  The truncation time for defining restricted mean time lost. Default is
  `NULL`, which corresponds to setting the truncation time as the
  maximum survival time in the data

- data:

  A data frame containing the variables in the outcome, propensity
  score, and IPCW models; a data frame with `n` rows (1 row per
  observation).

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
  include one of `'randomForest'` (survival outcomes only),
  `'boosting'`, `'logistic'` (survival outcomes only, fast), `'poisson'`
  (count outcomes only, fast), and `'gam'` (count outcomes only).
  Default is `NULL`, which assigns `'boosting'` for count outcomes and
  `'randomForest'` for survival outcomes.

- response:

  The type of response variables; `count` (default) or `survival`.

## Value

A list of elements: - y: outcome; vector of length `n` (observations) -
d : the event indicator; vector of length `n`; only if
`respone = "survival"` - trt: binary treatment; vector of length `n` -
x.ps: matrix of `p.ps` baseline covariates specified in the propensity
score model (plus intercept); dimension `n` by `p.ps + 1` - x.cate:
matrix of `p.cate` baseline covariates specified in the outcome model;
dimension `n` by `p.cate` - x.ipcw: matrix of `p.ipw` baseline covarites
specified in inverse probability of censoring weighting model; dimension
`n` by `p.ipw` - time: offset; vector of length `n`; only if
`response = "count"` - if `fun = "catefit"`: - prop: formatted
`prop.cutoff` - prop.no1: formatted `prop.cutoff` with 1 removed if
applicable; otherwise prop.no1 is the same as prop - if
`fun = "crossv"` - prop.onlyhigh: formatted `prop.cutoff` with 0 removed
if applicable - prop.bi; formatted `prop.cutoff` with 0 and 1 removed if
applicable - prop.multi: formatted `prop.multi`, starting with 0 and
ending with 1
