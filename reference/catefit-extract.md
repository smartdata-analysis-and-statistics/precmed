# Extract CATE scores and score coefficients

Extract scores for the original analysis data, coefficients, and, where
available, confidence intervals for score coefficients from a `catefit`.

## Usage

``` r
# S3 method for class 'catefit'
predict(object, newdata = NULL, method = NULL, ...)

# S3 method for class 'catefit'
coef(object, method = NULL, ...)

# S3 method for class 'catefit'
confint(object, parm, level = 0.95, method = "contrastReg", ...)
```

## Arguments

- object:

  An object of class `catefit`.

- newdata:

  New data are not currently supported; omit this argument to retrieve
  scores for the original analysis sample in its original order.

- method:

  A single scoring method used to fit the object, such as
  `"contrastReg"`. For
  [`predict()`](https://rdrr.io/r/stats/predict.html), it may be omitted
  if exactly one method is available. For
  [`coef()`](https://rdrr.io/r/stats/coef.html), omission returns the
  stored coefficient table. For
  [`confint()`](https://rdrr.io/r/stats/confint.html), the default is
  `"contrastReg"`.

- ...:

  Additional arguments; currently unused.

- parm:

  Coefficient names or positions; omitted selects all coefficients.

- level:

  Confidence level between 0 and 1; defaults to 0.95.

## Value

[`predict()`](https://rdrr.io/r/stats/predict.html) returns a numeric
vector of stored scores. [`coef()`](https://rdrr.io/r/stats/coef.html)
returns the stored coefficient table or a named numeric vector for a
selected method. [`confint()`](https://rdrr.io/r/stats/confint.html)
returns a matrix of Wald confidence limits for the selected
coefficients.

## Details

Scores for count and survival outcomes are on the log ratio scale;
binary scores are on the risk-difference scale. They describe fitted
scores, not validated individual treatment effects. Confidence intervals
require stored coefficient standard errors, currently available for
contrast regression in count and continuous fits only. They are not
available for survival or binary fits. Methods that do not fit linear
scores, such as random forests and boosting, do not provide score
coefficients.

## Author

Thomas Debray
