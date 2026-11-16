# Compute the area between curves from the `"precmed"` object

Compute the area between curves (ABC) for each scoring method in the
`"precmed"` object. This should be run only after results of
[`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)
have been obtained.

## Usage

``` r
# S3 method for class 'precmed'
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

Jiang, X., Simoneau, G., Pellegrini, F., Shen, C., de Moor, C., & Tian,
L. (2026). Visualization and interpretation of individualized treatment
rule results. In T. P. A. Debray, T.-L. Nguyen, & R. W. Platt (Eds.),
*Comparative Effectiveness and Personalized Medicine Research Using
Real-World Data* (Chapter 20, pp. 513–536). Chapman & Hall/CRC.
[doi:10.1201/9781003300809-20](https://doi.org/10.1201/9781003300809-20)

## See also

[`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)
function and [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`boxplot()`](https://rdrr.io/r/graphics/boxplot.html) methods for
`"precmed"` objects.

## Examples

``` r
# The result was created with:
# catecvCountExample <- catecv(
#   response = "count", data = countExample, score.method = "poisson",
#   cate.model = y ~ age + female + previous_treatment + previous_cost +
#     previous_number_relapses + offset(log(years)),
#   ps.model = trt ~ age + previous_treatment, higher.y = FALSE,
#   cv.n = 5, seed = 999, verbose = 0)
# Load the pre-fitted result instead of repeating the calculation.
data(catecvCountExample)
# Rows are scoring methods; columns are validation iterations.
abc(catecvCountExample)
#>               cv1       cv2       cv3        cv4       cv5
#> poisson 0.1375899 0.1172791 0.1120313 0.06136285 0.1877103
```
