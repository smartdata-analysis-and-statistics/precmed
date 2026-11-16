# Histogram of bootstrap estimates

Histogram of bootstrap estimates

## Usage

``` r
# S3 method for class 'atefit'
plot(
  x,
  bins,
  alpha = 1,
  title = waiver(),
  theme = precmed_plot_theme(),
  null.value = 0,
  metrics = NULL,
  ...
)
```

## Arguments

- x:

  An object of class `"atefit"`.

- bins:

  Number of bins

- alpha:

  Histogram bar opacity, between 0 and 1. Defaults to 1.

- title:

  The text for the title

- theme:

  A ggplot2 theme. By default, a customized minimal theme is used;
  supply, for example,
  [`ggplot2::theme_classic()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)
  to use another theme.

- null.value:

  Reference value on the plotted effect scale. Defaults to 0 (no risk
  difference or a ratio of 1); use `NULL` or `NA` to omit the null line.

- metrics:

  Character vector of metrics to plot, in display order. For binary
  fits, choose from `"risk0"`, `"risk1"`, and `"risk.difference"`. For
  count fits, choose from `"rate0"`, `"rate1"`, and `"log.rate.ratio"`;
  for survival fits, choose from `"rmst0"`, `"rmst1"`,
  `"log.rmtl.ratio"`, and `"log.cumulative.hazard.ratio"`. Older
  survival fits may instead contain `"log.hazard.ratio"`. The default
  `NULL` plots all metrics available in the object.

- ...:

  Other parameters

## Value

A plot of the class `ggplot`, displaying the estimated ATE across the
bootstrap samples

## Details

Create a histogram displaying the distribution of the bootstrap
estimates. The solid red vertical line represents the full-sample
estimate. A dashed black line marks `null.value` in risk-difference and
log-ratio panels; it is not shown in treatment-specific risk, rate or
RMST panels. Treatment-specific facets use the original treatment labels
when available. Older binary or count objects without treatment-specific
bootstrap estimates show only the risk difference or log rate ratio,
respectively.

## Author

Thomas Debray
