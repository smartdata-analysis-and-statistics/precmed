# Estimate an average treatment effect for a binary outcome

Estimate the population-average risk difference between two treatments
using augmented inverse probability weighting (AIPW). Bootstrap
resampling provides a standard error, confidence interval, and p-value.

## Usage

``` r
atefitbinary(
  data,
  outcome_model,
  outcome_method = "logistic",
  ps_model,
  ps_method = "glm",
  minPS = 0.01,
  maxPS = 0.99,
  n.boot = 500,
  seed = NULL,
  verbose = 0,
  ps.control = list(),
  ...
)
```

## Arguments

- data:

  A data frame containing the outcome, treatment, and covariates.

- outcome_model:

  A formula specifying the logistic outcome model, fitted separately in
  each treatment arm by default. Put the binary outcome on the left and
  baseline covariates on the right; do not include treatment. The
  outcome must be numeric or logical and coded 0/1. Missing values and
  offsets are not supported.

- outcome_method:

  Outcome-model method; only `"logistic"` (default) is supported.

- ps_model:

  A formula with treatment on the left and treatment-model covariates on
  the right. The model estimates the propensity score. Both treatment
  groups must be present. Numeric 0/1 or logical FALSE/TRUE map to
  treatment 0/1; other two-valued labels are sorted and mapped to 0/1 in
  that order. Use `trt ~ 1` for an intercept-only model.

- ps_method:

  Propensity score method: `"glm"` (default) or `"lasso"`, using the
  package's propensity-score LASSO procedure.

- minPS, maxPS:

  Lower and upper truncation limits for propensity scores. Defaults are
  0.01 and 0.99; `minPS` must be smaller than `maxPS`.

- n.boot:

  Number of bootstrap samples, an integer of at least 2. Default is 500.

- seed:

  An optional integer seed; `NULL` uses the current random number state.
  Supplying a seed sets R's random-number state.

- verbose:

  Accepted for compatibility with the other ATE functions; the binary
  implementation does not currently display progress messages.

- ps.control:

  Named options for LASSO propensity-score interaction screening.
  Supported entry: `solver = "C"` (default) or `"R"`.

- ...:

  Additional arguments. Deprecated argument names from precmed 1.1.0
  remain supported for backward compatibility and emit warnings.

## Value

An object of class `atefit`. See
[`atefit`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
for available methods and the contents of the fitted object.

## Details

The estimator combines predicted risks from logistic outcome models with
inverse-propensity-weighted residual corrections. Under appropriate
identification and regularity assumptions, it is consistent if either
the outcome model or treatment model is correctly specified. A causal
interpretation also requires consistency, no unmeasured confounding and
adequate treatment overlap. Propensity score truncation can affect the
guarantee based on the treatment model. AIPW risk estimates need not lie
in \[0, 1\]. Outcomes should be fully observed for the intended
follow-up horizon; this function does not adjust for censoring.

Each bootstrap replicate refits both models. If any replicate for the
risk difference fails or is non-finite, inference is unavailable;
failures are summarized in a warning and stored in
`bootstrap.diagnostics`.
[`diagnostics()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/diagnostics.md)
can be used to inspect the records. The point estimate and standard
error are returned by [`coef()`](https://rdrr.io/r/stats/coef.html) and
through the result object;
[`confint()`](https://rdrr.io/r/stats/confint.html) calculates Wald
intervals at a requested level.

## References

Hoogland, J., IntHout, J., Belias, M., et al. (2021). A tutorial on
individualized treatment effect prediction from randomized trials with a
binary endpoint. *Statistics in Medicine*, 40(26), 5961–5981.
[doi:10.1002/sim.9154](https://doi.org/10.1002/sim.9154)

## See also

[`atefit`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md),
[`simulate_binary`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/simulate_binary.md),
[`diagnostics`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/diagnostics.md)

## Author

Thomas Debray <tdebray@fromdatatowisdom.com>

## Examples

``` r
d <- simulate_binary(n = 500, seed = 42)
fit <- atefitbinary(d, outcome_model = y ~ age + female,
                    ps_model = trt ~ age, outcome_method = "logistic",
                    n.boot = 10, seed = 42)
coef(fit)
#> risk.difference 
#>      -0.0287478 
confint(fit)
#>                       2.5 %     97.5 %
#> risk.difference -0.09705359 0.03955799
```
