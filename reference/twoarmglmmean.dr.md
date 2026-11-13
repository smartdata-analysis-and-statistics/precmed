# Doubly robust estimators of the coefficients in the contrast regression as well as their covariance matrix

Solving the estimating equation `bar S_n (delta) = 0`

## Usage

``` r
twoarmglmmean.dr(y, x.cate, trt, ps, f1.predictor, f0.predictor)
```

## Arguments

- y:

  Observed outcome; vector of size `n`

- x.cate:

  Matrix of `p.cate` baseline covariates; dimension `n` by `p.cate`

- trt:

  Treatment received; vector of size `n` units with treatment coded as
  0/1

- ps:

  Estimated propensity scores for all observations; vector of size `n`

- f1.predictor:

  Initial predictions of the outcome (expected number of relapses for
  one unit of exposure time) conditioned on the covariates `x` for
  treatment group trt = 1; `mu_1(x)`, step 1 in the two regression;
  vector of size `n`

- f0.predictor:

  Initial predictions of the outcome (expected number of relapses for
  one unit of exposure time) conditioned on the covariates `x` for
  treatment group trt = 0; `mu_0(x)`, step 1 in the two regression;
  vector of size `n`

## Value

coef: Doubly robust estimators of the regression coefficients `delta_0`;
vector of size `p` + 1 (intercept included) vcov: Variance-covariance
matrix of the estimated coefficient `delta_0`; matrix of size `p` + 1 by
`p` + 1
