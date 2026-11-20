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

  A single score to extract from the fitted object. Use
  `"two_regression"` or `"contrast_regression"` for those procedures, or
  the outcome method (such as `"logistic"`, `"boosting"`, `"poisson"`,
  `"coxph"`, `"weibull"`, `"log_normal"`, `"negative_binomial"`, or
  `"random_forest"`) for a T-learner. The established stored names
  `"twoReg"`, `"contrastReg"`, `"negBin"`, and `"randomForest"` are also
  accepted; contrast regression can also be selected with
  `"contrastreg"`. If exactly one T-learner score is present,
  `"t_learner"` selects it. If exactly one score of any kind is present,
  `method` may be omitted from
  [`predict()`](https://rdrr.io/r/stats/predict.html). For
  [`coef()`](https://rdrr.io/r/stats/coef.html), omission returns the
  coefficient table. For
  [`confint()`](https://rdrr.io/r/stats/confint.html), the default is
  `"contrastReg"`. Coefficients and confidence intervals are available
  only for certain score methods.

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
selected method. When score models use different covariates, the table
has the union of coefficient names and missing entries for terms not
used by a method; selecting a method returns only its terms.
[`confint()`](https://rdrr.io/r/stats/confint.html) returns a matrix of
Wald confidence limits for the selected coefficients.

## Details

Scores for count and survival outcomes are on the log ratio scale;
binary scores are on the risk-difference scale. They describe fitted
scores, not validated individual treatment effects. Pass the fitted
`catefit` object to [`predict()`](https://rdrr.io/r/stats/predict.html),
not an individual outcome model in its `fit` component. When multiple
T-learners are fitted, select one by its outcome method. Confidence
intervals require stored coefficient standard errors, currently
available for contrast regression in count and continuous fits only.
They are not available for survival or binary fits. Methods that do not
fit linear scores, such as random forests, boosting, Cox, Weibull and
log-normal T-learners, do not provide linear score coefficients.

## Author

Thomas Debray
