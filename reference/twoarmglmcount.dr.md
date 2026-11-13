# Doubly robust estimators of the coefficients in the contrast regression as well as their covariance matrix and convergence information

Newton-Raphson algorithm is used to solve the estimating equation
`bar S_n (delta) = 0`

## Usage

``` r
twoarmglmcount.dr(
  y,
  x.cate,
  time,
  trt,
  ps,
  f1.predictor,
  f0.predictor,
  error.maxNR = 0.001,
  max.iterNR = 150,
  tune = c(0.5, 2)
)
```

## Arguments

- y:

  Observed outcome; vector of size `n`

- x.cate:

  Matrix of `p.cate` baseline covariates; dimension `n` by `p.cate`

- time:

  Log-transformed person-years of follow-up; vector of size `n`

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

- error.maxNR:

  A numerical value \> 0 indicating the minimum value of the mean
  absolute error in Newton Raphson algorithm. Used only if
  `score.method = 'contrastReg'`. Default is `0.001`.

- max.iterNR:

  A positive integer indicating the maximum number of iterations in the
  Newton Raphson algorithm. Used only if `score.method = 'contrastReg'`.
  Default is `150`.

- tune:

  A vector of 2 numerical values \> 0 specifying tuning parameters for
  the Newton Raphson algorithm. `tune[1]` is the step size, `tune[2]`
  specifies a quantity to be added to diagonal of the slope matrix to
  prevent singularity. Used only if `score.method = 'contrastReg'`.
  Default is `c(0.5, 2)`.

## Value

coef: Doubly robust estimators of the regression coefficients `delta_0`;
vector of size `p` + 1 (intercept included) vcov: Variance-covariance
matrix of the estimated coefficient `delta_0`; matrix of size `p` + 1 by
`p` + 1 converge: Indicator that the Newton Raphson algorithm converged
for `delta_0`; boolean
