# Doubly robust estimators of the coefficients in the contrast regression as well as their covariance matrix and convergence information

Newton-Raphson algorithm is used to solve the estimating equation
`bar S_n (delta) = 0`

## Usage

``` r
twoarmsurv.dr(
  ynew,
  dnew,
  trt,
  x.cate,
  tau0,
  weightsurv,
  ps,
  f1.predictor,
  f0.predictor,
  error.maxNR = 0.001,
  max.iterNR = 100,
  tune = c(0.5, 2),
  solver = c("C", "R"),
  optim.control = list()
)
```

## Arguments

- ynew:

  Truncated survival time; vector of size `n`

- dnew:

  Event indicator after truncation; vector of size `n`

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

- f1.predictor:

  Initial predictions of the outcome (restricted mean time loss)
  conditioned on the covariates `x.cate` for treatment group trt = 1;
  `mu_1(x.cate)`, step 1 in the two regression; vector of size `n`

- f0.predictor:

  Initial predictions of the outcome (restricted mean time loss)
  conditioned on the covariates `x.cate` for treatment group trt = 0;
  `mu_0(x.cate)`, step 1 in the two regression; vector of size `n`

- error.maxNR:

  A numerical value \> 0 indicating the minimum value of the mean
  absolute error in Newton Raphson algorithm. Used only if
  `score.method = 'contrastReg'`. Default is `0.001`.

- max.iterNR:

  A positive integer indicating the maximum number of iterations in the
  Newton Raphson algorithm. Used only if `score.method = 'contrastReg'`.
  Default is `100`.

- tune:

  A vector of 2 numerical values \> 0 specifying tuning parameters for
  the Newton Raphson algorithm. `tune[1]` is the step size, `tune[2]`
  specifies a quantity to be added to diagonal of the slope matrix to
  prevent singularity. Used only if `score.method = 'contrastReg'`.
  Default is `c(0.5, 2)`.

- solver:

  Numerical backend: `"C"` (default) or `"R"`. Both preserve the
  weighted estimating equation, first-iteration regularization and eta
  stopping rule. Small rounding differences are possible.

- optim.control:

  Named controls passed to
  [`stats::optim`](https://rdrr.io/r/stats/optim.html) if NR fails to
  converge, for example `list(maxit = 1000, reltol = 1e-10)`.
  `max.iterNR` governs NR; `optim.control$maxit` governs fallback only.

## Value

coef: Doubly robust estimators of the contrast regression coefficients
`delta_0`; vector of size `p.cate` + 1 (intercept included) converge:
Indicator that the Newton Raphson algorithm converged for `delta_0`;
boolean
