# Doubly robust estimator of the average treatment effect with Cox model for survival data

Doubly robust estimator of the average treatment effect between two
treatments, which is the restricted mean time lost (RMTL) ratio of
treatment 1 over treatment 0 for survival outcomes.

## Usage

``` r
drsurv(
  y,
  d,
  x.cate,
  x.ps,
  x.ipcw,
  trt,
  yf = NULL,
  tau0,
  surv.min = 0.025,
  ps.method = "glm",
  minPS = 0.01,
  maxPS = 0.99,
  ipcw.method = "breslow",
  rmst.control = list(),
  ps.control = list()
)
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

- x.ps:

  Matrix of `p.ps` baseline covariates specified in the propensity score
  model; dimension `n` by `p.ps`.

- x.ipcw:

  Matrix of `p.ipw` baseline covariate specified in inverse probability
  of censoring weighting; dimension `n` by `p.ipw`.

- trt:

  Treatment received; vector of size `n` with treatment coded as 0/1.

- yf:

  Follow-up time, interpreted as the potential censoring time; vector of
  size `n` if the potential censoring time is known.

- tau0:

  The truncation time for defining restricted mean time lost.

- surv.min:

  Lower truncation limit for probability of being censored (positive and
  very close to 0).

- ps.method:

  A character value for the method to estimate the propensity score.
  Allowed values include one of: `'glm'` for logistic regression with
  main effects only (default), or `'lasso'` for a logistic regression
  with main effects and LASSO penalization on two-way interactions
  (added to the model if interactions are not specified in `ps.model`).
  Relevant only when `ps.model` has more than one variable.

- minPS:

  A numerical value (in \`\[0, 1\]\`) below which estimated propensity
  scores should be truncated. Default is `0.01`.

- maxPS:

  A numerical value (in \`(0, 1\]\`) above which estimated propensity
  scores should be truncated. Must be strictly greater than `minPS`.
  Default is `0.99`.

- ipcw.method:

  The censoring model. Allowed values are: `'breslow'` (Cox regression
  with Breslow estimator of the baseline survivor function),
  `'aft (exponential)'`, `'aft (weibull)'`, `'aft (lognormal)'` or
  `'aft (loglogistic)'`. Default is `'breslow'`.

- rmst.control:

  Named options for Cox restricted mean survival integration. Supported
  entry: `solver = "C"` (default) or `solver = "R"`. The C backend
  accumulates each patient's integral without allocating a
  patient-by-event-time matrix. Cox model fitting is unchanged.
  Available only for survival outcomes; independent of
  `contrast.control`.

- ps.control:

  Named options for LASSO propensity-score interaction screening.
  Supported entry: `solver = "C"` (default) or `solver = "R"`. Only
  affects screening when `ps.method = "lasso"`; the LASSO fit and
  interaction selection rule are unchanged. Kept separate from
  boosting's `...` and other numerical controls.

## Value

Return a list of 4 elements:

- `rmst1`: A numeric value of the estimated restricted mean survival
  time n the group `trt = 1`.

- `rmst0`: A numeric value of the estimated restricted mean survival
  time n the group `trt = 0`.

- `log.rmtl.ratio`: A numeric value of the estimated log rmtl ratio.

- `log.hazard.ratio`: A numeric value of the estimated log hazard ratio.
