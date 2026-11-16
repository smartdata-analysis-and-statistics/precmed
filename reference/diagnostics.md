# Summarize recorded fitting diagnostics

Collect treatment-recoding information from
[`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md),
or recorded subgroup errors and warnings from count or survival
[`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
and binary, count or survival
[`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md),
into a data frame. Each row represents one recorded message. Survival
subgroup messages also remain in `object[["errors/warnings"]]`.

## Usage

``` r
diagnostics(object)
```

## Arguments

- object:

  A binary, count or survival `atefit`, a count or survival `catefit`,
  or a binary, count or survival `catecv` result.

## Value

A data frame with columns `method`, `split`, `iteration`, `group`,
`subgroup`, `type`, `measure`, and `message`. For count mutually
exclusive groups, `subgroup` gives the lower and upper boundaries;
survival uses the upper percentile label recorded by the estimator. For
`atefit`, bootstrap failures have split `bootstrap`, the replicate
number in iteration, and the affected estimate in measure. Split and
iteration are `NA` for `catefit`. Otherwise, measure is `NA` unless the
recorded message identifies it. Returns a zero-row data frame with the
same columns when there are no recorded messages.

## Details

This reports only diagnostics stored in the result, not every warning
emitted during fitting. In particular,
[`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
stores the treatment-recoding information and failed or undefined
bootstrap estimates, but not all bootstrap warnings. Count fits record
subgroup messages and support in `object$subgroup.support`. Continuous
results are not supported by this helper.

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
