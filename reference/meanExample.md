# Simulated data with a continuous outcome

A dataset containing a continuous outcome and 6 baseline covariates

## Usage

``` r
data(meanExample)
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

  count outcome, current number of relapses

## Examples

``` r
data(meanExample)
str(meanExample)
#> 'data.frame':    4000 obs. of  8 variables:
#>  $ age                     : num [1:4000, 1] 11.8 -8.2 -3.2 -16.2 -8.2 ...
#>   ..- attr(*, "scaled:center")= num 46.2
#>  $ female                  : int  1 0 0 1 1 0 1 1 1 1 ...
#>  $ previous_treatment      : Factor w/ 3 levels "drugA","drugB",..: 3 3 1 3 2 3 3 3 3 1 ...
#>  $ previous_cost           : num [1:4000, 1] -0.331 -0.575 0.519 -0.519 -0.567 ...
#>   ..- attr(*, "scaled:center")= num 13788
#>   ..- attr(*, "scaled:scale")= num 20150
#>  $ previous_number_symptoms: Factor w/ 3 levels "0","1",">=2": 2 2 3 2 2 2 2 2 1 3 ...
#>  $ previous_status_measure : num  6.15 2.003 3.786 1.889 -0.516 ...
#>  $ trt                     : chr  "drug1" "drug0" "drug0" "drug0" ...
#>  $ y                       : num [1:4000, 1] 21.93 -8.26 -3.32 -22.13 -7.91 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:4000] "1" "2" "3" "4" ...
#>   .. ..$ : NULL
```
