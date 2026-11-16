# Plot estimated conditional treatment-effect scores

Display the distribution of scores fitted by
[`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md).
These are in-sample predictions; their spread alone does not establish
treatment-effect heterogeneity or validate the score.

## Usage

``` r
# S3 method for class 'catefit'
plot(
  x,
  methods = NULL,
  type = c("histogram", "density"),
  bins = 30,
  null.value = 0,
  theme = precmed_plot_theme(),
  ...
)
```

## Arguments

- x:

  An object of class `"catefit"`.

- methods:

  Optional character vector of fitted scoring methods to show, in the
  desired order. The default shows all available scores.

- type:

  Either `"histogram"` (default) or `"density"`. Density curves require
  at least two distinct finite values per method.

- bins:

  Positive integer number of histogram bins; default 30.

- null.value:

  Position of a dashed reference line on the score scale; default 0. Use
  `NULL` or `NA` to omit the line.

- theme:

  A ggplot2 theme, defaulting to the package plot theme.

- ...:

  Currently unused.

## Value

A ggplot object. Non-finite scores are omitted with a warning.

## Details

Binary scores are risk differences, count scores are log rate ratios,
survival scores are log restricted mean time lost ratios, and continuous
scores are mean differences. The treatment contrast is treatment 1
versus treatment 0 (or treatment 1 minus treatment 0 for differences).
The methods may impose different model assumptions; see
[`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md).
Original treatment labels are shown when stored in the fit; older
objects without labels use treatment 1 and treatment 0. Histogram panels
show patient counts for each method. Density curves are normalized
separately and overlaid for comparison.

## Author

Thomas Debray <tdebray@fromdatatowisdom.com>

## Examples

``` r
# \donttest{
d <- simulate_binary(n = 300, seed = 1)
fit <- catefit("binary", d, "twoReg", y ~ age, trt ~ age,
               prop.cutoff = 1, B = 1, Kfold = 2, seed = 1, higher.y = FALSE)
plot(fit)

# }
```
