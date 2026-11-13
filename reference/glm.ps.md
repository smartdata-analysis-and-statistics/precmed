# Propensity score estimation with LASSO

Propensity score based on a multivariate logistic regression with LASSO
penalization on the two-way interactions

## Usage

``` r
glm.ps(trt, x.ps, xnew = NULL, minPS = 0.01, maxPS = 0.99)
```

## Arguments

- trt:

  Treatment received; vector of size `n` (observations) with treatment
  coded as 0/1

- x.ps:

  Matrix of `p.ps` baseline covariates (plus a leading column of 1 for
  the intercept); dimension `n` by `p.ps + 1` (covariates in the
  propensity score model plus intercept)

- xnew:

  Matrix of `p.ps` baseline covariates (plus a leading column of 1 for
  the intercept) for which we want propensity scores predictions;
  dimension `m` (observations in the new data set) by `p.ps + 1`

- minPS:

  A numerical value (in \`\[0, 1\]\`) below which estimated propensity
  scores should be truncated. Default is `0.01`.

- maxPS:

  A numerical value (in \`(0, 1\]\`) above which estimated propensity
  scores should be truncated. Must be strictly greater than `minPS`.
  Default is `0.99`.

## Value

The trimmed propensity score for each unit; vector of size `n` (if
`xnew` is NULL) or `m`
