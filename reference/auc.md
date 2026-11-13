# Compute the area under the curve using linear or natural spline interpolation

This function computes the area under the curve for two vectors where
one corresponds to the x values and the other corresponds to the y
values. It supports both linear and spline interpolation.

## Usage

``` r
auc(
  x,
  y,
  from = min(x, na.rm = TRUE),
  to = max(x, na.rm = TRUE),
  type = c("linear", "spline"),
  subdivisions = 100,
  ...
)
```

## Arguments

- x:

  A numeric vector of x values.

- y:

  A numeric vector of y values of the same length as x.

- from:

  The value from where to start calculating the area under the curve.
  Defaults to the smallest x value.

- to:

  The value from where to end the calculation of the area under the
  curve. Defaults to the greatest x value.

- type:

  The type of interpolation: "linear" or "spline". Defaults to "linear".

- subdivisions:

  An integer indicating how many subdivisions to use for \`integrate\`
  (for spline-based approximations).

- ...:

  Additional arguments passed on to \`approx\` (for linear
  interpolations).

## Value

A numeric value representing the area under the curve.
