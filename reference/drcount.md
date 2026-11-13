# Doubly robust estimator of the average treatment effect for count data

Doubly robust estimator of the average treatment effect between two
treatments, which is the rate ratio of treatment 1 over treatment 0 for
count outcomes.

## Usage

``` r
drcount(
  y,
  trt,
  x.cate,
  x.ps,
  time,
  ps.method = "glm",
  minPS = 0.01,
  maxPS = 0.99,
  interactions = TRUE,
  ps.control = list()
)
```

## Arguments

- y:

  A numeric vector of size `n` with each element representing the
  observed count outcome for each subject.

- trt:

  A numeric vector (in {0, 1}) of size `n` with each element
  representing the treatment received for each subject.

- x.cate:

  A numeric matrix of dimension `n` by `p.cate` with each column
  representing each baseline covariate specified in the outcome model
  for all subjects.

- x.ps:

  A numeric matrix of dimension `n` by `p.ps + 1` with a leading column
  of 1 as the intercept and each remaining column representing each
  baseline covariate specified in the propensity score model for all
  subjects.

- time:

  A numeric vector of size `n` with each element representing the
  log-transformed person-years of follow-up for each subject.

- ps.method:

  A character value for the method to estimate the propensity score.
  Allowed values include one of: `'glm'` for logistic regression with
  main effects only (default), or `'lasso'` for a logistic regression
  with main effects and LASSO penalization on two-way interactions
  (added to the model if interactions are not specified in `ps.model`).
  Relevant only when `ps.model` has more than one variable.

- minPS:

  A numerical value between 0 and 1 below which estimated propensity
  scores should be truncated. Default is `0.01`.

- maxPS:

  A numerical value between 0 and 1 above which estimated propensity
  scores should be truncated. Must be strictly greater than `minPS`.
  Default is `0.99`.

- interactions:

  A logical value indicating whether the outcome model should allow for
  treatment-covariate interaction by `x`. If `TRUE`, interactions will
  be assumed only if at least 10 patients received each treatment
  option. Default is `TRUE`.

- ps.control:

  Named options for LASSO propensity-score interaction screening.
  Supported entry: `solver = "C"` (default) or `solver = "R"`. Only
  affects screening when `ps.method = "lasso"`; the LASSO fit and
  interaction selection rule are unchanged. Kept separate from
  boosting's `...` and other numerical controls.

## Value

Return a list of 4 elements:

- `log.rate.ratio`: A numeric value of the estimated log rate ratio.

- `rate0`: A numeric value of the estimated rate in the group trt=0.

- `rate1`: A numeric value of the estimated rate in the group trt=1.
