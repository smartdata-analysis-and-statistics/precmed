# Estimate the ATE of the RMTL ratio and unadjusted hazard ratio in multiple bi-level subgroups defined by the proportions

If only care about the higher subgroup (above cutoff), only need
ate.rmtl.high and hr.high so set "onlyhigh" to be TRUE Scores are
adjusted to the opposite sign if `higher.y` == FALSE; scores stay the
same if `higher.y` == FALSE; this is because estsurv() function always
takes the subgroup of the top highest adjusted scores, and higher
adjusted scores should always represent high responders of trt=1

## Usage

``` r
estsurv.bilevel.subgroups(
  y,
  d,
  x.cate,
  x.ps,
  x.ipcw,
  trt,
  yf,
  tau0 = tau0,
  score,
  higher.y,
  prop,
  onlyhigh,
  surv.min = 0.025,
  ps.method = "glm",
  minPS = 0.01,
  maxPS = 0.99,
  ipcw.method = "breslow"
)
```

## Arguments

- y:

  Observed survival or censoring time; vector of size `n`.

- d:

  The event indicator, normally `1 = event, 0 = censored`; vector of
  size `n`.

- x.cate:

  Matrix of `p.cate` baseline covariates specified in the outcome model;
  dimension `n` by `p.cate`.

- x.ps:

  Matrix of `p.ps` baseline covariates specified in the propensity score
  model; dimension `n` by `p.ps`.

- x.ipcw:

  Matrix of `p.ipw` baseline covariate specified in inverse probability
  of censoring weighting; dimension `n` by `p.ipw`.

- trt:

  Treatment received; vector of size `n` with treatment coded as 0/1.

- yf:

  Follow-up time, interpreted as the potential censoring time; vector of
  size `n` if the potential censoring time is known.

- tau0:

  The truncation time for defining restricted mean time lost.

- score:

  Estimated log CATE scores for all `n` observations from one of the
  five methods (random forest, boosting, naive Poisson, two regressions,
  contrast regression); vector of size `n`.

- higher.y:

  A logical value indicating whether higher (`TRUE`) or lower (`FALSE`)

- prop:

  Proportions corresponding to percentiles in the estimated log CATE
  scores that define subgroups to calculate ATE for; vector of floats in
  \`(0, 1\]\` (if `onlyhigh=TRUE`) or in \`(0, 1)\` (if
  `onlyhigh=FALSE`): Each element of `prop` represents the high/low
  cutoff in each bi-level subgroup and the length of `prop` is number of
  bi-level subgroups

- onlyhigh:

  Indicator of returning only the ATEs in the higher-than-cutoff
  category of the bi-level subgroups; boolean.

- surv.min:

  Lower truncation limit for probability of being censored (positive and
  very close to 0).

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

- ipcw.method:

  The censoring model. Allowed values are: `'breslow'` (Cox regression
  with Breslow estimator of the baseline survivor function),
  `'aft (exponential)'`, `'aft (weibull)'`, `'aft (lognormal)'` or
  `'aft (loglogistic)'`. Default is `'breslow'`.

## Value

ate.rmtl.high: estimated ATEs (ratio of RMTL) in the multiple bi-level
subgroups that are in the higher-than-cutoff category; vector of size
equal to the length of prop; always returned. ate.rmtl.low: estimated
ATEs (ratio of RMTL) in the multiple bi-level subgroups that are in the
lower-than-cutoff category; vector of size equal to the length of prop;
returned only when `onlyhigh = TRUE`. hr.high: unadjusted hazard ratio
in the multiple bi-level subgroups that are in the higher-than-cutoff
category; vector of size equal to the length of prop; always returned.
hr.low: unadjusted hazard ratio in the multiple bi-level subgroups that
are in the lower-than-cutoff category; vector of size equal to the
length of prop; returned only when `onlyhigh = TRUE`
