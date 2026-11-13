# Doubly robust estimator of the average treatment effect for continuous data

Doubly robust estimator of the average treatment effect between two
treatments, which is the mean difference of treatment 1 over treatment 0
for continuous outcomes.

## Usage

``` r
drmean(
  y,
  trt,
  x.cate,
  x.ps,
  ps.method = "glm",
  minPS = 0.01,
  maxPS = 0.99,
  interactions = TRUE
)
```

## Arguments

- y:

  A numeric vector of size `n` with each element representing the
  observed continuous outcome for each subject.

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
  subjects

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

  A logical value indicating whether the outcome model should assume
  interactions between `x` and `trt`. If `TRUE`, interactions will be
  assumed only if at least 10 patients received each treatment option.
  Default is `TRUE`.

## Value

Return a list of 4 elements:

- `mean.diff`: A numeric value of the estimated mean difference.

- `mean.diff0`: A numeric value of the estimated mean difference in
  treatment group 0.

- `mean.diff1`: A numeric value of the estimated mean difference in
  treatment group 1.
