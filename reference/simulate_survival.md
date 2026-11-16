# Simulate an illustrative survival dataset

Generate fully synthetic data for demonstrating survival
treatment-effect methods. Baseline covariates are generated
independently of patient records. Treatment is assigned using a logistic
propensity model; the numbers in each treatment group arise from those
probabilities and vary by sample. Event times follow a Weibull model
with treatment interactions with age and previous relapses. By default,
censoring is independent of covariates, treatment, and event times.
Default coefficients are chosen to illustrate confounding and
heterogeneous treatment effects; they are not clinical estimates and do
not establish the relative performance of methods on other data.

## Usage

``` r
simulate_survival(
  n = 4000L,
  seed = 20261114L,
  ps_weights = c(intercept = 0.11, age = 1.15, drugB = 0.25, drugC = -0.2),
  lp_weights = c(intercept = 0, age = 0.45, female = 0.12, previous_cost = 0.1,
    previous_number_relapses = 0.16, treatment_age = -0.28, treatment_relapses = -0.16),
  conditional_hr = 0.9,
  censor_weights = c(intercept = log(1/450), age = 0, female = 0, previous_cost = 0,
    previous_number_relapses = 0, drugB = 0, drugC = 0, treatment = 0),
  shape = 1.15,
  scale = 650,
  treatment_shape = shape,
  hr_reference_time = scale
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
  at age 48 with zero previous relapses, at `hr_reference_time`. At that
  time, for other covariate values, the log HR is
  `log(conditional_hr) + treatment_age * (age_in_years - 48)/10 + treatment_relapses * previous_number_relapses`,
  where `age_in_years` is the uncentered age.

- censor_weights:

  Named coefficients for an exponential censoring model with log hazard
  per day. Recognized names are `intercept`, `age`, `female`,
  `previous_cost`, `previous_number_relapses`, `drugB`, `drugC`, and
  `treatment`. Named subsets override the defaults. Age is measured in
  decades from 48 years; cost uses the scaling described above; drugA
  and drug0 are the reference groups. The default intercept is
  `log(1/450)` and all other coefficients are zero.

- shape:

  Positive Weibull shape for drug0. Defaults to 1.15.

- scale:

  Positive Weibull baseline scale in days for drug0. Defaults to 650.

- treatment_shape:

  Positive Weibull shape for drug1. Defaults to `shape`, giving
  proportional hazards. Different shapes produce a time-varying hazard
  ratio.

- hr_reference_time:

  Positive time in days at which `conditional_hr` is defined. Defaults
  to `scale`.

## Value

A data frame with `n` rows and nine columns: centered `age`, `female`,
`previous_treatment`, scaled `previous_cost`,
`previous_number_symptoms`, `previous_number_relapses`, `trt`, observed
time `y`, and event indicator `d`. The age and cost columns retain
centering and scaling attributes. Times are in days; 1 indicates relapse
and 0 censoring.

## Details

At positive time `t`, the covariate-specific hazard ratio is its value
at `hr_reference_time` multiplied by
`(t / hr_reference_time)^(treatment_shape - shape)`. It increases over
time if the treated shape is larger, and decreases if it is smaller. By
default, random censoring times are exponential with mean 450 days;
administrative censoring occurs at 900 days. Censoring and event times
are generated independently conditional on covariates and treatment,
including when censoring depends on those variables. Extreme propensity
coefficients can produce poor treatment overlap or a sample with no
patients in one treatment group.

## Author

Thomas Debray <tdebray@fromdatatowisdom.com>

## Examples

``` r
toy <- simulate_survival(n = 200, seed = 1)
table(toy$trt)
#> 
#> drug0 drug1 
#>    85   115 
toy_unequal <- simulate_survival(n = 200, seed = 1,
  ps_weights = c(intercept = -1), conditional_hr = 0.8)
table(toy_unequal$trt)
#> 
#> drug0 drug1 
#>   145    55 
toy_censor <- simulate_survival(n = 200, seed = 1,
  censor_weights = c(treatment = log(2), age = 0.3))
toy_nonph <- simulate_survival(n = 200, seed = 1, treatment_shape = 1.65)
```
