# Estimate restricted mean survival time (RMST) based on Cox regression model

Estimate restricted mean survival time (RMST) based on Cox regression
model

## Usage

``` r
cox.rmst(y, d, x.cate, xnew, tau0)
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

- xnew:

  Matrix of `p.cate` baseline covariates for which we want an estimate
  of the RMST; dimension `m` (observations in the new data set) by
  `p.cate`

- tau0:

  The truncation time for defining restricted mean time lost.

## Value

The estimated RMST for new subjects with covariates `xnew`; vector of
size `m`.
