# Histogram of bootstrap estimates

Histogram of bootstrap estimates

## Usage

``` r
# S3 method for class 'atefit'
plot(x, bins, alpha = 0.7, title = waiver(), theme = theme_classic(), ...)
```

## Arguments

- x:

  An object of class `"atefit"`.

- bins:

  Number of bins

- alpha:

  Opacity

- title:

  The text for the title

- theme:

  Defaults to `theme_classic()`. Other options include `theme_grey()`,
  `theme_bw()`, `theme_light()`, `theme_dark()`, and `theme_void()`

- ...:

  Other parameters

## Value

A plot of the class `ggplot`, displaying the estimated ATE across the
bootstrap samples

## Details

Create a histogram displaying the distribution of the bootstrap
estimates. The red vertical reference line represents the final
estimate.

## Author

Thomas Debray
