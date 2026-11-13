# Simulated data with count outcome

A dataset containing a count outcome, a length of follow-up and 6
baseline covariates

## Usage

``` r
data(countExample)
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

- years:

  length of follow-up, in years

## Examples

``` r
data(countExample)
str(countExample)
#> 'data.frame':    4000 obs. of  9 variables:
#>  $ age                     : num [1:4000, 1] -21.22 -2.22 -10.22 6.78 4.78 ...
#>   ..- attr(*, "scaled:center")= num 46.2
#>  $ female                  : int  0 1 1 0 1 0 1 0 1 1 ...
#>  $ previous_treatment      : Factor w/ 3 levels "drugA","drugB",..: 1 1 1 1 3 3 3 3 3 1 ...
#>  $ previous_cost           : num [1:4000, 1] 0.451 -0.194 -0.534 -0.337 -0.271 ...
#>   ..- attr(*, "scaled:center")= num 13824
#>   ..- attr(*, "scaled:scale")= num 20172
#>  $ previous_number_symptoms: Factor w/ 3 levels "0","1",">=2": 2 2 2 2 2 2 2 3 2 1 ...
#>  $ previous_number_relapses: int  0 0 1 1 0 0 0 0 1 1 ...
#>  $ trt                     : Factor w/ 2 levels "drug0","drug1": 2 1 1 2 2 1 2 2 2 2 ...
#>  $ y                       : int  0 0 0 0 1 0 0 0 0 0 ...
#>  $ years                   : num  0.2847 0.6105 2.653 0.0383 2.2697 ...
rate <- countExample$y / countExample$years
```
