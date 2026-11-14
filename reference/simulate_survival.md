# Simulate an illustrative survival dataset

Generate fully synthetic data for demonstrating survival
treatment-effect methods. Baseline covariates are generated
independently of patient records. Treatment is assigned using a logistic
propensity model; the numbers in each treatment group arise from those
probabilities and vary by sample. Event times follow a Weibull
proportional hazards model with treatment interactions with age and
previous relapses. Censoring is independent of covariates and event
times. These data are illustrative and are not clinical estimates.

## Usage

``` r
simulate_survival(
  n = 4000L,
  seed = 20261114L,
  ps_weights = c(intercept = 0.11, age = 0.25, drugB = 0.25, drugC = -0.2),
  lp_weights = c(intercept = 0, age = 0.2, female = 0.12, previous_cost = 0.1,
    previous_number_relapses = 0.16, treatment_age = -0.28, treatment_relapses = -0.16),
  conditional_hr = exp(-0.35)
)
```

## Arguments

- n:

  Positive integer giving the total sample size. Default is 4,000.
  Treatment group sizes are determined by the propensity model.

- seed:

  An integer seed for reproducibility, or `NULL` to use the current
  random-number state. Supplying a seed sets R's random-number state.

- ps_weights:

  Named numeric coefficients for the logistic propensity model.
  Recognized names are `intercept`, `age`, `drugB`, and `drugC`; named
  subsets override the defaults. Age is measured in decades from 48
  years; previous-treatment drugA is the reference category. The default
  intercept yields approximately 54 split in the package's original
  survival example.

- lp_weights:

  Named numeric coefficients for the log event hazard. Recognized names
  are `intercept`, `age`, `female`, `previous_cost`,
  `previous_number_relapses`, and `treatment_age`, and
  `treatment_relapses`; named subsets override the defaults. Age is
  measured in decades from 48 years. The interaction coefficients
  describe changes in the log HR per decade and per previous relapse,
  respectively. Previous cost is centered at 14,362 and divided by
  24,266 dollars.

- conditional_hr:

  Positive conditional instantaneous hazard ratio for drug1 versus drug0
  at age 48 with zero previous relapses. For other values, the log HR is
  `log(conditional_hr) + treatment_age * (age_in_years - 48)/10 + treatment_relapses * previous_number_relapses`,
  where `age_in_years` is the uncentered age.

## Value

A data frame with `n` rows and nine columns: centered `age`, `female`,
`previous_treatment`, scaled `previous_cost`,
`previous_number_symptoms`, `previous_number_relapses`, `trt`, observed
time `y`, and event indicator `d`. The age and cost columns retain
centering and scaling attributes. Times are in days; 1 indicates relapse
and 0 censoring.

## Details

Event times use a Weibull baseline with shape 1.15 and scale 650 days.
Censoring times are exponential with mean 450 days and administrative
censoring at 900 days. Extreme propensity coefficients can produce poor
treatment overlap or a sample with no patients in one treatment group.

## Examples

``` r
toy <- simulate_survival(n = 200, seed = 1)
table(toy$trt)
#> 
#> drug0 drug1 
#>    88   112 
toy_unequal <- simulate_survival(n = 200, seed = 1,
  ps_weights = c(intercept = -1), conditional_hr = 0.8)
table(toy_unequal$trt)
#> 
#> drug0 drug1 
#>   143    57 
```
