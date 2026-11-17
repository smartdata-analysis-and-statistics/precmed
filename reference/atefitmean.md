# Doubly robust estimator of and inference for the average treatment effect for continuous data

Continuous-outcome average treatment effect estimation is not yet
implemented.

## Usage

``` r
atefitmean(
  data,
  outcome_model,
  ps_model,
  ps_method = "glm",
  minPS = 0.01,
  maxPS = 0.99,
  n.boot = 500,
  plot.boot = FALSE,
  seed = NULL,
  verbose = 0,
  ps.control = list(),
  ...
)
```

## Arguments

- data:

  A data frame containing the variables in the outcome and treatment
  models; a data frame with `n` rows (1 row per observation).

- outcome_model:

  A formula specifying the outcome model, with the continuous outcome on
  the left and baseline covariates on the right. Do not include
  treatment. Continuous-outcome fitting is not yet implemented.

- ps_model:

  A formula describing the treatment model (propensity-score model) to
  be fitted. The treatment must appear on the left-hand side. The
  treatment must be a numeric vector coded as 0/1. If data are from a
  RCT, specify `ps_model` as an intercept-only model.

- ps_method:

  A character value for the method to estimate the propensity score.
  Allowed values include one of: `'glm'` for logistic regression with
  main effects only (default), or `'lasso'` for a logistic regression
  with main effects and LASSO penalization on two-way interactions
  (added to the model if interactions are not specified in `ps_model`).
  Relevant only when `ps_model` has more than one variable.

- minPS:

  A numerical value between 0 and 1 below which estimated propensity
  scores should be truncated. Default is `0.01`.

- maxPS:

  A numerical value between 0 and 1 above which estimated propensity
  scores should be truncated. Must be strictly greater than `minPS`.
  Default is `0.99`.

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
  affects screening when `ps_method = "lasso"`; the LASSO fit and
  interaction selection rule are unchanged. Kept separate from
  boosting's `...` and other numerical controls.

- ...:

  Additional arguments. The deprecated named argument `cate.model` is
  accepted as an alias for `outcome_model`, with a warning. Supply only
  one of these arguments. The deprecated aliases `ps.model` and
  `ps.method` are also accepted with warnings; use `ps_model` and
  `ps_method` instead. Other arguments are currently unsupported.

## Value

Continuous-outcome fitting is not implemented; this function stops with
an error and does not return an object.

## Details

Continuous-outcome estimation is not currently available.

## Examples

``` r
# This module is not implemented yet!
```
