# Doubly robust estimators of the coefficients in the two regression

Doubly robust estimators of the coefficients in the two regression

## Usage

``` r
onearmglmcount.dr(y, x.cate, time, trt, ps, f.predictor)
```

## Arguments

- y:

  Observed outcome; vector of size `n`

- x.cate:

  Matrix of `p` baseline covariates; dimension `n` by `p`

- time:

  Log-transformed person-years of follow-up; vector of size `n`

- trt:

  Treatment received; vector of size `n` units with treatment coded as
  0/1

- ps:

  Estimated propensity scores for all observations; vector of size `n`

- f.predictor:

  Initial prediction of the outcome (expected number of relapses for one
  unit of exposure time) conditioned on the covariates `x` for one
  treatment group `r`; `mu_r(x)`, step 1 in the two regression; vector
  of size `n`

## Value

Doubly robust estimators of the regression coefficients `beta_r` in the
doubly robust estimating equation where `r = 0, 1` is treatment
received; vector of size `p` + 1 (intercept included)
