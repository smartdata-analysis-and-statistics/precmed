# Compute the area between curves from the `"precmed"` object

Compute the area between curves (ABC) for each scoring method in the
`"precmed"` object. This should be run only after results of
[`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)
have been obtained.

## Usage

``` r
abc(x, ...)

# Default S3 method
abc(x, ...)
```

## Arguments

- x:

  An object of class `"precmed"`.

- ...:

  Additional arguments (currently unused).

## Value

Returns a matrix of numeric values with number of columns equal to the
number cross-validation iteration and number of rows equal to the number
of scoring methods in `x`.

## Details

The ABC is the area between a validation curve and the overall ATE in
the validation set. It is calculated for each scoring method separately.
Higher ABC values are preferable as they indicate that more treatment
effect heterogeneity is captured by the scoring method. Negative values
of ABC are possible if segments of the validation curve cross the
overall ATE line. The ABC is calculated with the
[`auc()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/auc.md)
in `utility.R` with a natural cubic spline interpolation. The
calculation of the ABC is always based on validation curves based on 100
proportions equally spaced from `min(prop.cutoff)` to
`max(prop.cutoff)`.

The ABC is a metric to help users select the best scoring method in
terms of capturing treatment effect heterogeneity in the data. It should
be used in complement to the visual inspection of the validation curves
in the validation set in
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## References

Zhao, L., Tian, L., Cai, T., Claggett, B., & Wei, L. J. (2013).
*Effectively selecting a target population for a future comparative
study. Journal of the American Statistical Association, 108(502),
527-539.*

## See also

[`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)
function and [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`boxplot()`](https://rdrr.io/r/graphics/boxplot.html) methods for
`"precmed"` objects.

## Examples

``` r
# \donttest{
# Count outcome
cv_count <- catecv(response = "count",
                   data = countExample,
                   score.method = "poisson",
                   cate.model = y ~ age + female + previous_treatment +
                                previous_cost + previous_number_relapses +
                                offset(log(years)),
                   ps.model = trt ~ age + previous_treatment,
                   higher.y = FALSE, cv.n = 5, verbose = 1)
#> Warning: Variable trt was recoded to 0/1 with drug0->0 and drug1->1.
#>   |                                                                              |                                                                      |   0%
#> cv = 1 
#>   splitting the data..
#>   training..
#>   validating..
#> 
#> cv = 2 
#>   splitting the data..
#>   training..
#>   validating..
#> 
#> cv = 3 
#>   splitting the data..
#>   training..
#>   validating..
#> 
#> cv = 4 
#>   splitting the data..
#>   training..
#>   validating..
#> 
#> cv = 5 
#>   splitting the data..
#>   training..
#>   validating..
#> 
#> Total runtime : 18.45 secs 

# ABC of the validation curves for each method and each CV iteration
abc(cv_count)
#>               cv1       cv2        cv3        cv4       cv5
#> poisson 0.1475806 0.1155643 0.09207781 0.09615332 0.1261193

# Survival outcome
library(survival)
cv_surv <- catecv(response = "survival",
                  data = survivalExample,
                  score.method = c("poisson", "randomForest"),
                  cate.model = Surv(y, d) ~ age + female + previous_cost +
                               previous_number_relapses,
                  ps.model = trt ~ age + previous_treatment,
                  higher.y = FALSE,
                  cv.n = 5)
#> Warning: No value supplied for tau0. Default sets tau0 to the maximum survival time.
#> Warning: Variable trt was recoded to 0/1 with drug0->0 and drug1->1.

# ABC of the validation curves for each method and each CV iteration
abc(cv_surv)
#>                    cv1       cv2       cv3       cv4       cv5
#> poisson      0.1857850 0.1420101 0.1904015 0.1879060 0.1942061
#> randomForest 0.1806077 0.1394080 0.1891241 0.1822773 0.1817726

# }
```
