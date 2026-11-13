# Doubly robust estimator of and inference for the average treatment effect for continuous data

Doubly robust estimator of the average treatment effect between two
treatments, which is the rate ratio of treatment 1 over treatment 0 for
count outcomes. Bootstrap is used for inference.

## Usage

``` r
atefitmean(
  data,
  cate.model,
  ps.model,
  ps.method = "glm",
  minPS = 0.01,
  maxPS = 0.99,
  interactions = TRUE,
  n.boot = 500,
  plot.boot = FALSE,
  seed = NULL,
  verbose = 0,
  ps.control = list()
)
```

## Arguments

- data:

  A data frame containing the variables in the outcome and propensity
  score models; a data frame with `n` rows (1 row per observation).

- cate.model:

  A formula describing the outcome model to be fitted. The outcome must
  appear on the left-hand side.

- ps.model:

  A formula describing the propensity score model to be fitted. The
  treatment must appear on the left-hand side. The treatment must be a
  numeric vector coded as 0/1. If data are from a RCT, specify
  `ps.model` as an intercept-only model.

- ps.method:

  A character value for the method to estimate the propensity score.
  Allowed values include one of: `'glm'` for logistic regression with
  main effects only (default), or `'lasso'` for a logistic regression
  with main effects and LASSO penalization on two-way interactions
  (added to the model if interactions are not specified in `ps.model`).
  Relevant only when `ps.model` has more than one variable.

- minPS:

  A numerical value between 0 and 1 below which estimated propensity
  scores should be truncated. Default is `0.01`.

- maxPS:

  A numerical value between 0 and 1 above which estimated propensity
  scores should be truncated. Must be strictly greater than `minPS`.
  Default is `0.99`.

- interactions:

  A logical value indicating whether the outcome model should be fitted
  separately by treatment arm with the variables in `cate.model`, which
  is equivalent to assuming treatment-covariate interaction by all of
  the variables in `cate.model`. If `TRUE`, the outcome model will be
  fitted separately by treatment arms only if at least 10 patients
  received each treatment option. Default is `TRUE`.

- n.boot:

  A numeric value indicating the number of bootstrap samples used.
  Default is `500`.

- plot.boot:

  A logical value indicating whether histograms of the bootstrapped
  treatment effect estimates should be produced at every `n.boot/10`-th
  iteration and whether the final histogram should be outputted. Default
  is `FALSE`.

- seed:

  An optional integer specifying an initial randomization seed for
  reproducibility. Default is `NULL`, corresponding to no seed.

- verbose:

  An integer value indicating whether intermediate progress messages and
  histograms should be printed. `1` indicates messages are printed and
  `0` otherwise. Default is `0`.

- ps.control:

  Named options for LASSO propensity-score interaction screening.
  Supported entry: `solver = "C"` (default) or `solver = "R"`. Only
  affects screening when `ps.method = "lasso"`; the LASSO fit and
  interaction selection rule are unchanged. Kept separate from
  boosting's `...` and other numerical controls.

## Value

Return a list of 8 elements:

- `log.rate.ratio`: A numeric value of the estimated log rate ratio.

- `se.boot.log.rate.ratio`: A numeric value of the bootstrap standard
  error of log rate ratio.

- `rate.ratio`: A numeric value of the estimated rate ratio.

- `rate.ratio0`: A numeric value of the estimated rate in the group
  trt=0.

- `rate.ratio1`: A numeric value of the estimated rate in the group
  trt=1.

- `rate.ratio.CIl`: A numeric value of the lower limit 95% bootstrap
  confidence interval for estimated rate ratio.

- `rate.ratio.CIu`: A numeric value of the upper limit 95% bootstrap
  confidence interval for estimated rate ratio.

- `pvalue`: A numeric value of the p-value derived from the bootstrapped
  values based on a Chi-squared distribution.

- `warning`: A warning message produced if the treatment variable was
  not coded as 0/1. The key to map the original coding of the variable
  to a 0/1 key is displayed in the warning to facilitate the
  interpretation of the remaining of the output.

- `plot`: If `plot.boot` is `TRUE`, a histogram displaying the
  distribution of the bootstrapped log rate ratios. The red vertical
  reference line in the histogram represents the estimated log rate
  ratio.

## Details

This helper function estimates the average treatment effect (ATE)
between two treatment groups in a given dataset specified by
`y, trt, x.cate, x.ps, time`. The ATE is estimated with a doubly robust
estimator that accounts for imbalances in covariate distributions
between the two treatment groups with inverse probability treatment
weighting. For count outcomes, the estimated ATE is the estimated rate
ratio between treatment 1 versus treatment 0. Both original and
log-transformed ATEs are returned, as well as the rate in either
treatment group. If `inference = TRUE`, the variability of the estimated
rate ratio is also calculated using bootstrap. Additional variability
outputs include standard error of the log rate ratio, 95% confidence
interval of the rate ratio, p-value, and a histogram of the log rate
ratio.

## Examples

``` r
# This module is not implemented yet!
```
