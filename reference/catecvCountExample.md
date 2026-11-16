# Pre-fitted count outcome validation example

A reproducible result of applying
[`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)
to `countExample` using the Poisson score and five validation
iterations. This saved result can be used to illustrate plots and area
between curves without refitting the model. The regeneration script is
in `data-raw/catecvCountExample.R` in the package source repository.

## Usage

``` r
data(catecvCountExample)
```

## Format

An object of class `precmed` returned by
[`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md).

## Examples

``` r
data(catecvCountExample)
abc(catecvCountExample)
#>               cv1       cv2       cv3        cv4       cv5
#> poisson 0.1375899 0.1172791 0.1120313 0.06136285 0.1877103
```
