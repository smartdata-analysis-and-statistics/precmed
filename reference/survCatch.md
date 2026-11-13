# Catch errors and warnings when estimating the ATEs in the nested subgroup

Storing the errors and warnings that occurred when estimating the ATEs
in the nested subgroups. If there are no errors and no warnings, the
estimated log.rmtl.ratio and log.hazard.ratio are provided. If there are
warnings but no errors, the estimated log.rmtl.ratio and
log.hazard.ratio are provided with a warning attribute set. If there are
errors, the NA values are returned for log.rmtl.ratio and
log.hazard.ratio. A error attribute set is also provided.

## Usage

``` r
survCatch(fun)
```

## Arguments

- fun:

  The drsurv function...

## Value

A list containing
