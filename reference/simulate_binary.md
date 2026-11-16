# Simulate a binary response from relapse counts

Generate synthetic relapse counts over a common follow-up period with
[`simulate_count()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/simulate_count.md),
then classify patients as responders when their count is at most a
chosen threshold. The data-generating model and its coefficients are
illustrative, not clinical estimates.

## Usage

``` r
simulate_binary(n = 4000L, seed = 20261114L, followup = 1, threshold = 0L, ...)
```

## Arguments

- n:

  Positive integer sample size. Default is 4,000.

- seed:

  Integer seed, or `NULL` to use the current random-number state.
  Supplying a seed sets R's random-number state.

- followup:

  Positive finite common follow-up in years. Default is 1.

- threshold:

  Non-negative integer maximum number of relapses for a response. The
  default, 0, defines a relapse-free responder.

- ...:

  Additional arguments to
  [`simulate_count()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/simulate_count.md),
  such as `conditional_rr`, `ps_weights`, `lp_weights`, or `dispersion`.
  Do not supply `followup_range`; follow-up is fixed by `followup`.

## Value

A data frame with the same covariates, treatment and `years` as
[`simulate_count()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/simulate_count.md),
a 0/1 response `y`, and the underlying count `relapse_count`. A value of
1 means the count is at most the threshold. Treatment labels are `drug0`
and `drug1`.

## Details

The binary outcome is fully observed at the common horizon; there is no
censoring. A treatment effect on the relapse rate generally implies a
different effect on the response probability, depending on follow-up,
the threshold, covariates and count distribution. The model does not
reconstruct any patient population or individual patient records.

## Author

Thomas Debray <tdebray@fromdatatowisdom.com>

## Examples

``` r
d <- simulate_binary(n = 500, seed = 1, followup = 1, threshold = 0)
table(d$trt, d$y)
#>        
#>           0   1
#>   drug0 110 114
#>   drug1 124 152
```
