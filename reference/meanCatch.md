# Catch errors and warnings when estimating the ATEs in the nested subgroup for continuous data

Storing the errors and warnings that occurred when estimating the ATEs
in the nested subgroups. If there are no errors and no warnings, the
estimated mean difference is provided. If there are warnings but no
errors, the estimated mean difference is provided with a warning
attribute set. If there are errors, the NA values are returned for mean
difference. A error attribute set is also provided.

## Usage

``` r
meanCatch(fun)
```

## Arguments

- fun:

  The drsurv function...

## Value

A list containing
