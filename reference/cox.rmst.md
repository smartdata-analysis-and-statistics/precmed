# Estimate restricted mean survival time (RMST) based on Cox regression model

Estimate restricted mean survival time (RMST) based on Cox regression
model

## Usage

``` r
cox.rmst(y, d, x.cate, xnew, tau0, solver = c("C", "R"), return_survival = FALSE)
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

- solver:

  Integration backend: `"C"` (default) or `"R"`. Both use the same
  fitted Cox model and step-function integration; the C implementation
  avoids the patient-by-event-time survival matrix.

- return_survival:

  Whether to also return predicted survival at `tau0`.

## Value

The estimated RMST for new subjects with covariates `xnew`; a vector of
size `m`. If `return_survival = TRUE`, a list with `rmst` and `survival`
vectors.
