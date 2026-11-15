# Simulate an illustrative count-outcome dataset

Generate fully synthetic data for demonstrating count treatment-effect
methods. Treatment is assigned by a logistic propensity model. Relapse
counts follow a Poisson or negative binomial distribution with a
log-rate model, treatment interactions and varying follow-up.
Coefficients are illustrative, not clinical estimates or a
reconstruction of the MarketScan data.

## Usage

``` r
simulate_count(
  n = 4000L,
  seed = 20261114L,
  ps_weights = c(intercept = 0.11, age = 1.15, drugB = 0.25, drugC = -0.2),
  lp_weights = c(intercept = log(0.5), age = 0.25, female = 0.12, previous_cost = 0.1,
    previous_number_relapses = 0.16, treatment_age = -0.2, treatment_relapses = -0.12),
  conditional_rr = 0.9,
  followup_range = c(0.5, 3),
  dispersion = Inf
)
```

## Arguments

- n:

  Positive integer sample size. Default is 4,000.

- seed:

  Integer seed, or `NULL` to use the current random-number state.
  Supplying a seed sets R's random-number state.

- ps_weights:

  Named coefficients for treatment assignment. Recognized names are
  `intercept`, `age`, `drugB`, and `drugC`; named subsets override the
  defaults. Age is measured in decades from 48 years, and drugA is the
  reference previous treatment.

- lp_weights:

  Named coefficients for the log relapse rate. Recognized names are
  `intercept`, `age`, `female`, `previous_cost`,
  `previous_number_relapses`, `treatment_age`, and `treatment_relapses`;
  named subsets override the defaults. Age is measured in decades from
  48 years and cost is centered at 13,824 and divided by 20,172 US
  dollars. The intercept is the log annual rate under drug0 at the
  reference covariate values.

- conditional_rr:

  Positive rate ratio for drug1 versus drug0 at age 48 and zero previous
  relapses. Elsewhere its log is
  `log(conditional_rr) + treatment_age * (age_in_years - 48)/10 + treatment_relapses * previous_number_relapses`.

- followup_range:

  Two finite positive values specifying the minimum and maximum
  follow-up in years. Follow-up is sampled uniformly and independently
  of covariates, treatment and the outcome process; equal values give
  constant follow-up.

- dispersion:

  Positive negative-binomial size parameter, or `Inf` (the default) for
  Poisson counts. For finite size, the count variance is
  `mean + mean^2 / dispersion`.

## Value

A data frame with centered `age`, `female`, `previous_treatment`, scaled
`previous_cost`, `previous_number_symptoms`, `previous_number_relapses`,
`trt`, count outcome `y`, and follow-up in years `years`. Age and cost
retain centering and scaling attributes. Treatment labels are `drug0`
and `drug1`.

## Details

The model assumes a constant individual event rate over follow-up;
conditional on covariates, treatment and follow-up, expected counts
equal the annual rate multiplied by follow-up. Extreme coefficients may
cause poor overlap or invalid rates. The generator is for illustration,
not a calibrated representation of any patient population.

## Author

Thomas Debray <tdebray@fromdatatowisdom.com>

## Examples

``` r
toy <- simulate_count(n = 500, seed = 1)
table(toy$trt)
#> 
#> drug0 drug1 
#>   224   276 
toy_overdispersed <- simulate_count(n = 500, seed = 1,
  conditional_rr = 0.8, dispersion = 2)
```
