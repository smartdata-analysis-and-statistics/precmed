# Probability of being censored

Probability of being censored which is used to correct the effect of
right censoring.

## Usage

``` r
ipcw.surv(
  y,
  d,
  x.ipcw,
  yf = NULL,
  ipcw.method = "breslow",
  tau0,
  surv.min = 0.025
)
```

## Arguments

- y:

  Observed survival or censoring time; vector of size `n`.

- d:

  The event indicator, normally `1 = event, 0 = censored`; vector of
  size `n`.

- x.ipcw:

  Matrix of `p.ipw` baseline covariate specified in inverse probability
  of censoring weighting; dimension `n` by `p.ipw`.

- yf:

  Follow-up time, interpreted as the potential censoring time; vector of
  size `n` if the potential censoring time is known. If unknown, set
  `yf == NULL` and `yf` will be taken as `y` in the function.

- ipcw.method:

  The censoring model. Allowed values are: `'breslow'` (Cox regression
  with Breslow estimator of the baseline survivor function),
  `'aft (exponential)'`, `'aft (weibull)'`, `'aft (lognormal)'` or
  `'aft (loglogistic)'`. Default is `'breslow'`.

- tau0:

  The truncation time for defining restricted mean time lost.

- surv.min:

  Lower truncation limit for probability of being censored (positive and
  very close to 0).

## Value

A vector of size `n` with the estimated probabilities
`Pr(C > min(y, tau0) | x.ipcw)`
