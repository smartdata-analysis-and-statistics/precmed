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
diagnostics(
  object,
  type = c("messages", "propensity", "overlap"),
  layout = c("formatted", "arms", "metrics")
)
```

## Arguments

- object:

  A binary, count or survival `atefit`, a count or survival `catefit`,
  or a binary, count or survival `catecv` result.

- type:

  Choose `"messages"` (the default) for recorded issues or
  `"propensity"` for full-sample treatment-model diagnostics from a
  binary, count or survival `atefit` result. `"overlap"` is an alias for
  `"propensity"`.

- layout:

  For overlap diagnostics, `"formatted"` (the default) returns a
  compact, rounded character table; `"arms"` returns one numeric row per
  treatment; `"metrics"` returns one numeric column per treatment and
  one row per diagnostic. Not used for messages.

## Value

A data frame with columns `method`, `split`, `iteration`, `group`,
`subgroup`, `type`, `measure`, and `message`. For count mutually
exclusive groups, `subgroup` gives the lower and upper boundaries;
survival uses the upper percentile label recorded by the estimator. For
`atefit`, bootstrap failures have split `bootstrap`, the replicate
number in iteration, and the affected estimate in measure. Split and
iteration are `NA` for `catefit`. Otherwise, measure is `NA` unless the
recorded message identifies it. Returns a zero-row data frame with the
same columns when there are no recorded messages. For
`type = "propensity"`, the default formatted data frame has one row per
diagnostic and one character column per treatment. Its percentages use
the number of patients in the corresponding treatment group. Use
`layout = "arms"` or `layout = "metrics"` for unrounded numeric results.
The underlying diagnostics report sample size, untruncated propensity
score range and median, number outside the other arm's observed score
range, number truncated below and above the specified limits, and the
median, maximum and effective sample size of the inverse treatment
probability weights after truncation. The effective sample size is
\\(\sum w)^2 / \sum w^2\\. These are treatment weights only, not
censoring weights or the complete weights of a particular estimator.
With `layout = "formatted"` or `layout = "metrics"`, treatment columns
use the original labels when available, or treatment 0 and treatment 1
otherwise. Older fitted objects do not contain this information.

## Details

Overlap and weights describe the full analysis sample and the fitted
treatment model; they do not establish covariate balance, positivity, or
absence of unmeasured confounding. Extreme weights or little common
support merit investigation, but no universal cutoff guarantees
validity. The message view reports only diagnostics stored in the
result, not every warning emitted during fitting. In particular,
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
diagnostics(output_atefit, type = "propensity")
diagnostics(output_atefit, type = "propensity", layout = "arms")
diagnostics(output_catecv)
subset(diagnostics(output_catecv), split == "validation")
} # }
```
