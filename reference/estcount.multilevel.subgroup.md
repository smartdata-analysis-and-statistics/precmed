# Estimate the ATE of the log RR ratio in one multilevel subgroup defined by the proportions

Scores are adjusted to the opposite sign if `higher.y` == FALSE; scores
stay the same if `higher.y` == TRUE; this is because subgroups defined
in estcount.multilevel.subgroup() start from the lowest to the highest
adjusted scores, and higher adjusted scores should always represent high
responders of trt=1

## Usage

``` r
estcount.multilevel.subgroup(
  y,
  x.cate,
  x.ps,
  time,
  trt,
  score,
  higher.y,
  prop,
  ps.method = "glm",
  minPS = 0.01,
  maxPS = 0.99,
  ps.control = list()
)
```

## Arguments

- y:

  Observed outcome; vector of size `n` (observations)

- x.cate:

  Matrix of `p.cate` baseline covariates; dimension `n` by `p.cate`
  (covariates in the outcome model)

- x.ps:

  Matrix of `p.ps` baseline covariates (plus a leading column of 1 for
  the intercept); dimension `n` by `p.ps + 1` (covariates in the
  propensity score model plus intercept)

- time:

  Log-transformed person-years of follow-up; vector of size `n`

- trt:

  Treatment received; vector of size `n` units with treatment coded as
  0/1

- score:

  Estimated log CATE scores for all `n` observations from one of the
  four methods (boosting, naive Poisson, two regressions, contrast
  regression); vector of size `n`

- higher.y:

  A logical value indicating whether higher (`TRUE`) or lower (`FALSE`)
  values of the outcome are more desirable. Default is `TRUE`.

- prop:

  Proportions corresponding to percentiles in the estimated log CATE
  scores that define subgroups to calculate ATE for; vector of floats in
  \`\[0, 1\]\` always starting with 0 and ending with 1: Each element of
  `prop` represents inclusive cutoffs in the multilevel subgroup and the
  length of `prop` is number of categories in the multilevel subgroup

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

- ps.control:

  Named options for LASSO propensity-score interaction screening.
  Supported entry: `solver = "C"` (default) or `solver = "R"`. Only
  affects screening when `ps.method = "lasso"`; the LASSO fit and
  interaction selection rule are unchanged. Kept separate from
  boosting's `...` and other numerical controls.

## Value

estimated ATEs of all categories in the one multilevel subgroup; vector
of size equal to the length of categories in the multilevel subgroup
