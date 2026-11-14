# Extract average treatment effect estimates and confidence intervals

Extract treatment effects from an `atefit` object on their native
estimation scale: the log rate ratio for count outcomes, or the log RMTL
ratio, conditional log hazard ratio, and log cumulative hazard ratio for
survival outcomes. Exponentiate the returned estimates and interval
limits to obtain ratios.

## Usage

``` r
# S3 method for class 'atefit'
coef(object, ...)

# S3 method for class 'atefit'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  An object of class `atefit`.

- ...:

  Additional arguments; currently unused.

- parm:

  A character vector of effect names or a numeric vector of effect
  positions. If omitted, all treatment effects are included.

- level:

  Confidence level, a single number strictly between 0 and 1. Default is
  `0.95`.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html) returns a named numeric
vector of treatment effect estimates.
[`confint()`](https://rdrr.io/r/stats/confint.html) returns a matrix
with one row per selected treatment effect and columns for the lower and
upper confidence limits.

## Details

Confidence intervals are normal (Wald) intervals calculated using the
bootstrap standard errors stored in the fitted object. Changing `level`
does not repeat the bootstrap. Treatment-specific event rates and
restricted mean survival times are not included.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- atefitcount(data = countExample,
  cate.model = y ~ age + offset(log(years)), ps.model = trt ~ age)
exp(coef(fit))
exp(confint(fit))
exp(confint(fit, level = 0.90))
} # }
```
