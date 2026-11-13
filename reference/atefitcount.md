# Doubly robust estimator of and inference for the average treatment effect for count data

Doubly robust estimator of the average treatment effect between two
treatments, which is the rate ratio for count outcomes. Bootstrap is
used for inference.

## Usage

``` r
atefitcount(
  data,
  cate.model,
  ps.model,
  ps.method = "glm",
  minPS = 0.01,
  maxPS = 0.99,
  interactions = TRUE,
  n.boot = 500,
  seed = NULL,
  verbose = 0
)
```

## Arguments

- data:

  A data frame containing the variables in the outcome, propensity
  score, and inverse probability of censoring models (if specified); a
  data frame with `n` rows (1 row per observation).

- cate.model:

  A formula describing the outcome model to be fitted. The outcome must
  appear on the left-hand side.

- ps.model:

  A formula describing the propensity score (PS) model to be fitted. The
  treatment must appear on the left-hand side. The treatment must be a
  numeric vector coded as 0 or 1. If data are from a randomized
  controlled trial, specify `ps.model = ~1` as an intercept-only model.

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

  A logical value indicating whether the outcome model should assume
  treatment-covariate interaction by `x`. If `TRUE`, interactions will
  be assumed only if at least 10 patients received each treatment
  option. Default is `TRUE`.

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

## Value

Return an item of the class `atefit` with the following elements:

- `log.rate.ratio`: A vector of numeric values of the estimated ATE
  (expressed as a log rate ratio of `trt=1` over `trt=0`), the bootstrap
  standard error, the lower and upper limits of 95% confidence interval,
  and the p-value.

- `rate0`: A numeric value of the estimated rate in the group `trt=0`.

- `rate1`: A numeric value of the estimated rate in the group `trt=1`.

- `trt.boot`: Estimated log rate ratios in each bootstrap sample.

- `warning`: A warning message produced if the treatment variable was
  not coded as 0 or 1. The key to map the original coding of the
  variable to a 0-1 coding is displayed in the warning to facilitate the
  interpretation of the remaining of the output.

## Details

This helper function estimates the average treatment effect (ATE)
between two treatment groups in a given dataset. The ATE is estimated
with a doubly robust estimator that accounts for imbalances in covariate
distributions between the two treatment groups with inverse probability
treatment weighting. For count outcomes, the estimated ATE is the
estimated rate ratio between treatment 1 versus treatment 0.

## Examples

``` r
output <- atefitcount(data = countExample,
                      cate.model = y ~ age + female + previous_treatment +
                                   previous_cost + previous_number_relapses +
                                   offset(log(years)),
                      ps.model = trt ~ age + previous_treatment,
                      verbose = 1, n.boot = 50, seed = 999)
#> Warning: Variable trt was recoded to 0/1 with drug0->0 and drug1->1.
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
#> Warning: Variable trt was recoded to 0/1 with drug0->0 and drug1->1.
plot(output)
```
