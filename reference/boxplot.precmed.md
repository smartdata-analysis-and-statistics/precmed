# A set of box plots of estimated ATEs from the `"precmed"` object

Provides box plots which depict distributions of estimated ATEs for each
multi-category subgroup in the validation set across all
cross-validation iterations. The subgroups are mutually exclusive and
are categorized by the CATE score percentiles (`prop.multi` specified in
[`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)
or
[`catecvmean()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecvmean.md)).
Box plots of mutually exclusive subgroups are constructed separately by
scoring method specified in
[`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md).
This should be run only after results of
[`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)
or
[`catecvmean()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecvmean.md))
have been obtained.

## Usage

``` r
# S3 method for class 'precmed'
boxplot(
  x,
  ylab = NULL,
  plot.hr = FALSE,
  title = waiver(),
  theme = theme_classic(),
  hazard.type = c("cumulative", "conditional"),
  measure = "rmtl",
  ...
)
```

## Arguments

- x:

  An object of class `"precmed"`.

- ylab:

  A character value for the y-axis label to describe what the ATE is.
  Default is `NULL`, which creates a default y-axis label based on
  available data.

- plot.hr:

  Deprecated; use `measure` instead. Whether to plot a hazard-based
  ratio for survival outcomes instead of the RMTL ratio. When `TRUE`,
  the default is the marginal cumulative hazard ratio.

- title:

  The text for the title

- theme:

  Defaults to `theme_classic()`. Other options include `theme_grey()`,
  `theme_bw()`, `theme_light()`, `theme_dark()`, and `theme_void()`

- hazard.type:

  Deprecated; use `measure` instead. Select `"cumulative"` (default) or
  the supplementary `"conditional"` Cox HR when `plot.hr = TRUE`. Older
  objects without cumulative hazard ratios require
  `hazard.type = "conditional"` or refitting.

- measure:

  For survival outcomes, choose `"rmtl"` (default), `"cumulative_hr"`
  (cumulative hazard ratio), or `"conditional_hr"` (conditional Cox
  hazard ratio fitted through `tau0`). Earlier measure names remain
  supported as aliases.

- ...:

  Other parameters

## Value

Returns sets of box plots, one set for each scoring method, over each of
the multi-category subgroups. A gray horizontal dashed line of the
overall ATE is included as a reference.

## Details

[`boxplot()`](https://rdrr.io/r/graphics/boxplot.html) takes in outputs
from
[`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)
and generates the box plots of estimated ATEs for multi-category
subgroups of the validation set. The box plots together with the overall
ATE reference line can help compare the scoring methods' ability to
distinguish subgroups of patients with different treatment effects.

For a given scoring method, box plots showing increasing or decreasing
trends across the multi-category subgroups indicate presence of
treatment effect heterogeneity (and the ability of the scoring method to
capture it). On the contrary, box plots which are relatively aligned
across the multi-category subgroups indicate absence of treatment effect
heterogeneity (or the inability of the scoring method to capture it).

## References

Yadlowsky, S., Pellegrini, F., Lionetto, F., Braune, S., & Tian, L.
(2020). *Estimation and validation of ratio-based conditional average
treatment effects using observational data. Journal of the American
Statistical Association, 1-18.* DOI: 10.1080/01621459.2020.1772080.

Jiang, X., Simoneau, G., Pellegrini, F., Shen, C., de Moor, C., & Tian,
L. (2026). Visualization and interpretation of individualized treatment
rule results. In T. P. A. Debray, T.-L. Nguyen, & R. W. Platt (Eds.),
*Comparative Effectiveness and Personalized Medicine Research Using
Real-World Data* (Chapter 20, pp. 513–536). Chapman & Hall/CRC.
[doi:10.1201/9781003300809-20](https://doi.org/10.1201/9781003300809-20)

## See also

[`plot`](https://rdrr.io/r/graphics/plot.default.html) and
[`abc()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/abc.md)
for `"precmed"` objects.

## Examples

``` r
data(catecvCountExample)
boxplot(catecvCountExample, ylab = "Rate ratio in each subgroup")

```
