# Doubly robust estimator of and inference for the average treatment effect for count data

Doubly robust estimator of the average treatment effect between two
treatments, which is the rate ratio for count outcomes. Bootstrap is
used for inference.

## Usage

``` r
atefitcount(
  data,
  outcome_model,
  outcome_method = "poisson",
  ps_model,
  ps_method = "glm",
  minPS = 0.01,
  maxPS = 0.99,
  n.boot = 500,
  seed = NULL,
  verbose = 0,
  ps.control = list(),
  ...
)
```

## Arguments

- data:

  A data frame containing the variables in the outcome, treatment and
  censoring models (if specified); a data frame with `n` rows (1 row per
  observation).

- outcome_model:

  A formula specifying the count outcome model, fitted separately in
  each treatment arm by default. Put the count outcome on the left and
  baseline covariates on the right; do not include treatment. Include an
  exposure offset where applicable. See `outcome_method` for the model
  family.

- outcome_method:

  Outcome-model method: `"poisson"` (default) or `"negative_binomial"`
  (also accepted as `"negbin"`). Negative binomial models estimate
  dispersion separately by treatment arm. The method is refitted in each
  bootstrap sample. Failed fits are recorded in bootstrap diagnostics.

- ps_model:

  A formula describing the treatment model (propensity-score model) to
  be fitted. The treatment must appear on the left-hand side. The
  treatment must be a numeric vector coded as 0 or 1. If data are from a
  randomized controlled trial, specify `ps_model = trt ~ 1` as an
  intercept-only model.

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

- seed:

  An optional integer specifying an initial randomization seed for
  reproducibility. Default is `NULL`, corresponding to no seed.

- verbose:

  An integer value indicating whether intermediate progress messages
  should be printed. `1` indicates messages are printed and `0`
  otherwise. Default is `0`.

- ps.control:

  Named options for LASSO propensity-score interaction screening.
  Supported entry: `solver = "C"` (default) or `solver = "R"`. Only
  affects screening when `ps_method = "lasso"`; the LASSO fit and
  interaction selection rule are unchanged. Kept separate from
  boosting's `...` and other numerical controls.

- ...:

  Additional arguments. Deprecated argument names from precmed 1.1.0
  remain supported for backward compatibility and emit warnings.

## Value

An object of class `atefit`. See
[`atefit`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
for available methods and the contents of the fitted object.

## Details

This helper function estimates the average treatment effect (ATE)
between two treatment groups in a given dataset. The ATE is estimated
with a doubly robust estimator that accounts for imbalances in covariate
distributions between the two treatment groups with inverse probability
treatment weighting. For count outcomes, the estimated ATE is the
estimated rate ratio between treatment 1 versus treatment 0.

Failed or non-finite bootstrap estimates are recorded in
`bootstrap.diagnostics`, with replicate number, affected measure and
reason. A warning summarizes failures; inference for each affected
measure is unavailable rather than based on a subset of replicates. See
also
[`diagnostics()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/diagnostics.md)
to inspect these records.

## Examples

``` r
output <- atefitcount(data = countExample,
                      outcome_model = y ~ age + female + previous_treatment +
                                   previous_cost + previous_number_relapses +
                                   offset(log(years)),
                      ps_model = trt ~ age + previous_treatment,
                      verbose = 1, n.boot = 50, seed = 999)
#> Variable trt was recoded to 0/1 with drug0->0 and drug1->1.
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   2%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   6%  |                                                                              |======                                                                |   8%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  12%  |                                                                              |==========                                                            |  14%  |                                                                              |===========                                                           |  16%  |                                                                              |=============                                                         |  18%  |                                                                              |==============                                                        |  20%  |                                                                              |================                                                      |  22%  |                                                                              |=================                                                     |  24%  |                                                                              |===================                                                   |  27%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  31%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  35%  |                                                                              |==========================                                            |  37%  |                                                                              |===========================                                           |  39%  |                                                                              |=============================                                         |  41%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  45%  |                                                                              |=================================                                     |  47%  |                                                                              |==================================                                    |  49%  |                                                                              |====================================                                  |  51%  |                                                                              |=====================================                                 |  53%  |                                                                              |=======================================                               |  55%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  59%  |                                                                              |===========================================                           |  61%  |                                                                              |============================================                          |  63%  |                                                                              |==============================================                        |  65%  |                                                                              |===============================================                       |  67%  |                                                                              |=================================================                     |  69%  |                                                                              |==================================================                    |  71%  |                                                                              |===================================================                   |  73%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  78%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  82%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  88%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  92%  |                                                                              |==================================================================    |  94%  |                                                                              |===================================================================   |  96%  |                                                                              |===================================================================== |  98%  |                                                                              |======================================================================| 100%
output
#> Average treatment effect:
#> 
#>                  estimate        SE   CI.lower  CI.upper    pvalue
#> log.rate.ratio 0.06311746 0.1624482 -0.2552752 0.3815101 0.6976173
#> 
#> Estimated event rates:
#> 
#>        estimate
#> rate1 0.2860635
#> rate0 0.2685660
#> 
#> Variable trt was recoded to 0/1 with drug0->0 and drug1->1.
plot(output)
```
