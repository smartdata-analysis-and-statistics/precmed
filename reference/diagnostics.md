# Summarize recorded fitting diagnostics

Collect treatment-recoding information from
[`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md),
or recorded subgroup errors and warnings from count or survival
[`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
and
[`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md),
into a data frame. Each row represents one recorded message. Survival
subgroup messages also remain in `object[["errors/warnings"]]`.

## Usage

``` r
diagnostics(object)
```

## Arguments

- object:

  A count or survival `atefit`, `catefit`, or `catecv` result.

## Value

A data frame with columns `method`, `split`, `iteration`, `group`,
`subgroup`, `type`, `measure`, and `message`. Count mutually exclusive
groups use lower and upper boundaries; survival uses the upper
percentile label. Method, split, iteration, group, subgroup, and measure
are `NA` for `atefit`; split and iteration are `NA` for `catefit`.
Otherwise, measure is `NA` unless the recorded message identifies it.
Returns a zero-row data frame with the same columns when there are no
recorded messages.

## Details

This reports only diagnostics stored in the result, not every warning
emitted during fitting. In particular,
[`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
stores the treatment-recoding information, but not all bootstrap
warnings. Count fits record subgroup support in
`object$subgroup.support`. Continuous results are not supported by this
helper.

## Author

Thomas Debray <tdebray@fromdatatowisdom.com>

## Examples

``` r
if (FALSE) { # \dontrun{
diagnostics(output_atefit)
diagnostics(output_catecv)
subset(diagnostics(output_catecv), split == "validation")
} # }
```
