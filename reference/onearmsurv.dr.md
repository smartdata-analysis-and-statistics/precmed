# Doubly robust estimators of the coefficients in the two regression

Doubly robust estimators of the coefficients in the two regression

## Usage

``` r
onearmsurv.dr(ynew, dnew, trt, x.cate, tau0, weightsurv, ps, f.predictor)
```

## Arguments

- ynew:

  Truncated survival or censoring time; vector of size `n`.

- dnew:

  The event indicator after truncation,
  `1 = event or censored after truncation, 0 = censored before truncation`;
  vector of size `n`.

- trt:

  Treatment received; vector of size `n` with treatment coded as 0/1.

- x.cate:

  Matrix of `p.cate` baseline covariates specified in the outcome model;
  dimension `n` by `p.cate`.

- tau0:

  The truncation time for defining restricted mean time lost.

- weightsurv:

  Estimated inverse probability of censoring weights with truncation for
  all observations; vector of size `n`.

- ps:

  Estimated propensity scores for all observations; vector of size `n`

- f.predictor:

  Initial prediction of the outcome (restricted mean time loss)
  conditioned on the covariates `x.cate` for one treatment group `r`;
  `mu_r(x.cate)`, step 1 in the two regression; vector of size `n`

## Value

Doubly robust estimators of the two regression coefficients `beta_r`
where `r = 0, 1` is treatment received; vector of size `p.cate` + 1
(intercept included)
