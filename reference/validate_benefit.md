# Validate predicted treatment effects

This experimental internal function is under development and issues a
warning when called. It is not called by
[`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md).
Summarize calibration of predicted treatment effects and the population
effect of a treatment rule using independent validation data.

## Usage

``` r
validate_benefit(
  y,
  trt,
  benefit,
  groups = 5,
  threshold = 0,
  higher.y = TRUE,
  propensity = NULL
)
```

## Arguments

- y:

  Numeric observed binary or continuous outcome.

- trt:

  Numeric treatment assignment, coded 0 for control and 1 for treatment.

- benefit:

  Numeric predicted treatment effect for each participant, on the
  outcome difference scale (treatment minus control).

- groups:

  Positive integer number of quantile groups for calibration. Tied
  predictions are kept together, so fewer groups may result.

- threshold:

  Numeric treatment-effect threshold on the same scale as `benefit`.
  Treatment is recommended when predicted benefit exceeds it.

- higher.y:

  Whether larger outcomes are desirable. If `FALSE`, the treatment rule
  recommends treatment when the predicted difference is less than
  `threshold`. Effects are still reported as treatment minus control.

- propensity:

  Optional numeric vector of treatment probabilities, one per
  participant, strictly between 0 and 1. If supplied, inverse
  probability weighting is used; if NULL, unadjusted arm means are used.

## Value

A list with `mean.bias` (observed treatment difference minus mean
prediction), `calibration` (group sizes, arm counts, mean predictions,
observed differences and standard errors), `calibration.fit`
(group-level intercept, slope, R-squared and RMSE), and `policy` (rule
treatment proportion and outcome differences for the rule versus its
opposite, control for all, and treatment for all). Policy differences
are oriented so positive values indicate better outcomes under the rule.

## Details

Unadjusted estimators assume randomized treatment assignment with a
common, positive probability of each treatment. With propensity scores,
identification requires conditional exchangeability given the
treatment-model covariates, positivity and a correctly specified
treatment model. Both approaches require consistency and no interference
for causal interpretation. Predictions should be obtained from
independent data or appropriate cross-validation. Censoring and
treatment costs are not adjusted for. A threshold changes the treatment
rule, but policy results remain outcome differences and do not subtract
treatment costs.

Calibration groups are formed from empirical quantile cut points; ties
are never split. Group summaries are differences in arm means, using
normalized inverse probability weights when propensity scores are
provided. The fitted line and RMSE give equal weight to each group and
describe group-level approximations, not individual treatment effect
prediction error. A group lacking either treatment arm has an
unavailable observed difference; in that case the calibration fit is
unavailable. The policy estimates standardize arm-specific outcomes
within the two rule categories to their observed proportions, allowing
unequal randomization probabilities. A required empty treatment arm
makes the corresponding policy result unavailable. Unadjusted standard
errors require at least two observations in each arm. Weighted standard
errors are unavailable: uncertainty in estimated treatment models is not
accounted for by the unadjusted formula. Truncated or misspecified
propensity scores can bias weighted estimates.

## References

Hoogland et al. (2021),
[doi:10.1002/sim.9154](https://doi.org/10.1002/sim.9154) ; Efthimiou et
al. (2023), [doi:10.1002/sim.9665](https://doi.org/10.1002/sim.9665) .

## Author

Thomas Debray
