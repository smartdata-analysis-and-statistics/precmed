# Extract estimates for score-defined subgroups

Extract treatment effects or treatment-specific outcomes for the nested
subgroups in a binary, count or survival `catefit` object.

## Usage

``` r
subgroup_estimates(object, method = NULL, measure = NULL)
```

## Arguments

- object:

  A binary, count or survival `catefit` object.

- method:

  One scoring method used to fit the object, such as
  `"contrast_regression"` or a T-learner outcome method such as
  `"poisson"`. The same method selectors as
  [`predict()`](https://rdrr.io/r/stats/predict.html) are accepted; see
  [`catefit-extract`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit-extract.md).
  If omitted, it is inferred when exactly one method is available.

- measure:

  For binary outcomes, `"risk_difference"` (risk under treatment 1 minus
  treatment 0). For count outcomes, `"rate_ratio"`. For survival
  outcomes, one of `"rmtl"` (RMTL ratio), `"cumulative_hr"` (marginal
  cumulative hazard ratio), `"conditional_hr"` (conditional Cox hazard
  ratio), `"risk0"`, `"risk1"`, `"risk_difference"` (risk under
  treatment 1 minus treatment 0), `"absolute_risk_reduction"` (the
  opposite contrast), `"rmst0"`, or `"rmst1"` at `tau0`. The default is
  `"risk_difference"` for binary, `"rate_ratio"` for count, and `"rmtl"`
  for survival.

## Value

A named numeric vector, with one estimate per nested subgroup in
increasing order of `prop.cutoff`. Treatment-specific risks and RMSTs
are outcomes, not treatment effects. Undefined estimates are `NA`.

## Details

The subgroups are ranked by the selected score and estimated using the
same data used to fit it; these results do not validate the score. The
cumulative hazard ratio and risk estimates use outcome, treatment, and
censoring models; a causal interpretation requires appropriate treatment
and censoring assumptions. Older fits without the requested estimates
must be refitted. No confidence intervals are returned. The measure
names are also accepted in dotted form.

## Author

Thomas Debray
