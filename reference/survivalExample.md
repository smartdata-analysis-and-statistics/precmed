# Simulated data with survival outcome

A dataset containing a time-to-event outcome, an event indicator,
treatment group, and 6 baseline covariates

## Usage

``` r
data(survivalExample)
```

## Format

A dataframe with 4000 rows (patients) and 9 variables:

- age:

  age at baseline, centered to 48 years old, in years

- female:

  sex, 0 for male, 1 for female

- previous_treatment:

  previous treatment, "drugA", "drugB", or "drugC"

- previous_cost:

  previous medical cost, in US dollars

- previous_number_symptoms:

  previous number of symptoms, "0", "1", or "\>=2"

- previous_number_relapses:

  previous number of relapses

- trt:

  current treatment, "drug0" or "drug1"

- y:

  time to first relapse or censoring

- d:

  event indicator, 1: relapse, 0: censored

## Examples

``` r
data(survivalExample)
str(survivalExample)
#> 'data.frame':    4000 obs. of  9 variables:
#>  $ age                     : num [1:4000, 1] -26.83 -2.83 -12.83 -20.83 6.17 ...
#>   ..- attr(*, "scaled:center")= num 45.8
#>  $ female                  : int  1 1 1 0 0 1 1 1 1 1 ...
#>  $ previous_treatment      : Factor w/ 3 levels "drugA","drugB",..: 3 1 3 1 3 1 3 3 3 1 ...
#>  $ previous_cost           : num [1:4000, 1] -0.341 -0.239 -0.171 -0.392 -0.518 ...
#>   ..- attr(*, "scaled:center")= num 14362
#>   ..- attr(*, "scaled:scale")= num 24266
#>  $ previous_number_symptoms: Factor w/ 3 levels "0","1",">=2": 2 2 2 3 3 2 3 2 2 2 ...
#>  $ previous_number_relapses: int  0 1 2 0 1 0 0 2 0 0 ...
#>  $ trt                     : Factor w/ 2 levels "drug0","drug1": 1 1 1 1 2 1 2 2 2 2 ...
#>  $ y                       : num  102.1 51.3 104.9 61.5 35.6 ...
#>  $ d                       : num  0 0 0 0 0 0 0 0 0 1 ...
```
