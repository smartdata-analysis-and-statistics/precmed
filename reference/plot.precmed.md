# Two side-by-side line plots of validation curves from the `"precmed"` object

Provides validation curves in two side-by-side plots, visualizing the
estimated ATEs in a series of nested subgroups in the training set and
validation set separately, where each line represents one scoring method
specified in
[`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)
or
[`catecvmean()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecvmean.md).
This should be run only after results of
[`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)
or
[`catecvmean()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecvmean.md)
have been obtained.

## Usage

``` r
# S3 method for class 'precmed'
plot(
  x,
  cv.i = NULL,
  combine = "mean",
  show.abc = TRUE,
  valid.only = FALSE,
  plot.hr = FALSE,
  ylab = NULL,
  legend.position = "bottom",
  xlim = NULL,
  title = waiver(),
  theme = precmed_plot_theme(),
  hazard.type = c("cumulative", "conditional"),
  measure = "rmtl",
  show.overall = TRUE,
  layout = "set",
  ...
)
```

## Arguments

- x:

  An object of class `"precmed"`.

- cv.i:

  A positive integer indicating the index of the CV iteration results to
  be plotted. Allowed values are: a positive integer \\\<=\\ `cv.n` in
  [`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)
  or `NULL`. If `cv.i = NULL`, the results across all CV iterations are
  combined according to `combine` and then plotted. Default is `NULL`.

- combine:

  A character value indicating how to combine the estimated ATEs across
  all CV iterations into a validation curve for each nested subgroup,
  separately for the training and validation results. Allowed values
  are: `'mean'` or `'median'`. Used only if `cv.i = NULL`. Default is
  `'mean'`.

- show.abc:

  A logical value indicating whether to show the ABC statistics in the
  validation set. Used only if `x$abc = TRUE` and `xlim` is not limited
  to a smaller range (i.e., `xlim = NULL` or equal to the entire
  `x$prop.onlyhigh` range). If `cv.i` is NULL, ABC statistics will be
  based on the combined CV iterations. If `cv.i` is an integer, ABC
  statistics will be based solely on that CV iteration. Default is
  `TRUE`.

- valid.only:

  A logical value indicating whether only the validation curves in the
  validation set should be plotted (`TRUE`). Otherwise, the validation
  curves in both the training and validation sets are plotted
  side-by-side (`FALSE`). Default is `FALSE`.

- plot.hr:

  Deprecated; use `measure` instead. Whether to plot a hazard-based
  ratio for survival outcomes instead of the RMTL ratio. When `TRUE`,
  the default is the marginal cumulative hazard ratio.

- ylab:

  A character value for the y-axis label to describe what the ATE is.
  Default is `NULL`, which creates a default y-axis label based on
  available data.

- legend.position:

  A character value for the legend position argument to be passed to
  `ggplot` object. Default is `'bottom'`.

- xlim:

  A numeric value for the range of the x-axis. Default is `NULL`, which
  means there is no range specified.

- title:

  The text for the title

- theme:

  A ggplot2 theme. By default, a customized minimal theme is used;
  supply, for example,
  [`ggplot2::theme_classic()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)
  to use another theme.

- hazard.type:

  Deprecated; use `measure` instead. Select `"cumulative"` (default) or
  the supplementary `"conditional"` Cox HR when `plot.hr = TRUE`. Older
  objects without cumulative hazard ratios require
  `hazard.type = "conditional"` or refitting.

- measure:

  For survival outcomes, a vector of `"rmtl"` (default),
  `"cumulative_hr"` (cumulative hazard ratio), and `"conditional_hr"`
  (conditional Cox hazard ratio fitted through `tau0`). Multiple
  measures appear in separate panels. The earlier names `"chr"`,
  `"cph"`, `"cumulative_hazard_ratio"`, and `"conditional_hazard_ratio"`
  remain supported as aliases.

- show.overall:

  Whether to show the dashed overall treatment effect reference line in
  each set. Defaults to `TRUE`.

- layout:

  Use `"set"` (default) for panels by dataset and curves by method, or
  `"method"` for panels by scoring method and training and validation
  curves in each panel. Overall-effect reference lines in the latter
  layout are dotted and colored by dataset.

- ...:

  Other parameters

## Value

A ggplot of training and validation curves, faceted by set or method.
Overall-effect reference lines are included when `show.overall = TRUE`.
ABC statistics appear in the legend or method headings when
`show.abc = TRUE`.

## Details

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) takes in
outputs from
[`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)
and generates two plots of validation curves side-by-side, one for the
training set and one for validation set. Separate validation curves are
produced for each scoring method specified via `score.method` in
[`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)
or
[`catecvmean()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecvmean.md).

The validation curves (and ABC statistics, if applicable) can help
compare the performance of different scoring methods in terms of
discerning potential treatment heterogeneity in subgroups with internal
validation. Steeper validation curves in the validation set suggest
presence of treatment effect heterogeneity (and the ability of the
scoring methods to capture it) while flat validation curves indicate
absence of treatment effect heterogeneity (or inability of the scoring
method to capture it).

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

[`abc()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/abc.md)
and [`boxplot()`](https://rdrr.io/r/graphics/boxplot.html) for
`"precmed"` objects.

## Examples

``` r
data(catecvCountExample)
plot(catecvCountExample)

plot(catecvCountExample, show.abc = FALSE)

plot(catecvCountExample, cv.i = 2)

```
