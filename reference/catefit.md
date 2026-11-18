# Estimate conditional treatment effect scores

Estimate CATE scores for binary, count and survival outcomes. Binary
scores include arm-specific logistic regression, doubly robust
two-regression estimation and boosting. Count and survival methods
include both outcome-regression and doubly robust scores. Continuous
outcomes are not yet implemented.

## Usage

``` r
catefit(
  response,
  data,
  outcome_model,
  outcome_method = NULL,
  benefit_method,
  ps.model,
  ps.method = "glm",
  init.model = NULL,
  initial.predictor.method = NULL,
  ipcw.model = NULL,
  ipcw.method = "breslow",
  minPS = 0.01,
  maxPS = 0.99,
  followup.time = NULL,
  tau0 = NULL,
  higher.y = TRUE,
  prop.cutoff = seq(0.5, 1, length = 6),
  surv.min = 0.025,
  xvar.smooth.score = NULL,
  xvar.smooth.init = NULL,
  tree.depth = 2,
  n.trees.rf = 1000,
  n.trees.boosting = 200,
  B = 3,
  Kfold = 5,
  error.maxNR = 0.001,
  max.iterNR = 150,
  tune = c(0.5, 2),
  seed = NULL,
  plot.gbmperf = FALSE,
  verbose = 0,
  contrast.control = list(),
  rmst.control = list(),
  ps.control = list(),
  rf.control = list(),
  ...
)
```

## Arguments

- response:

  A string describing the type of outcome in the data. Allowed values
  include "binary", "count" (see
  [`catecvcount()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecvcount.md)),
  "survival" (see
  [`catecvsurv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecvsurv.md))
  and "continuous" (see
  [`catecvmean()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecvmean.md)).

- data:

  A data frame containing the variables in the outcome, treatment, and
  censoring models (if specified); a data frame with `n` rows (1 row per
  observation).

- outcome_model:

  A formula describing the outcome and score covariates. The outcome
  must appear on the left-hand side. For survival outcomes, a `Surv`
  object must be used to describe the outcome.

- outcome_method:

  A character vector with one entry per `benefit_method`. For T-learners
  the available methods are `"logistic"` or `"boosting"` for binary
  outcomes, `"poisson"`, `"negative_binomial"`, or `"boosting"` for
  count outcomes, and `"poisson"`, `"boosting"`, or `"random_forest"`
  for survival outcomes. For two-regression or contrast regression it
  selects the initial outcome models: binary supports only `"logistic"`,
  count supports `"poisson"`, `"boosting"`, or `"gam"`, and survival
  supports `"logistic"`, `"boosting"`, or `"random_forest"`. When both
  doubly robust methods are requested, they must use the same initial
  outcome method. An explicitly supplied `initial.predictor.method` must
  agree with it. Repeated configurations that would produce the same
  score name are not supported.

- benefit_method:

  A vector of score-construction methods: `"t_learner"`,
  `"two_regression"`, or `"contrast_regression"`. Each is paired with
  the corresponding entry in `outcome_method`. Contrast regression is
  not available for binary outcomes. Continuous outcomes are not
  implemented. The previous method names `"twoReg"` and `"contrastReg"`
  are accepted as aliases for `"two_regression"` and
  `"contrast_regression"`, respectively.

- ps.model:

  A formula describing the treatment model (propensity-score model) to
  be fitted. The treatment must appear on the left-hand side. For binary
  outcomes, logical treatment values `FALSE` and `TRUE` are mapped to 0
  and 1, respectively. Numeric 0/1 values or two distinct treatment
  labels are also accepted. If data are from a randomized controlled
  trial, specify `ps.model = trt ~ 1` as an intercept-only model.

- ps.method:

  A character value for the method to estimate the propensity score.
  Allowed values include one of: `'glm'` for logistic regression with
  main effects only (default), or `'lasso'` for a logistic regression
  with main effects and LASSO penalization on two-way interactions
  (added to the model if interactions are not specified in `ps.model`).
  Relevant only when `ps.model` has more than one variable.

- init.model:

  Optional formula for the initial treatment-specific outcome models
  used by `twoReg` and `contrastReg` (binary outcomes support `twoReg`
  only). The outcome must match `outcome_model`. The default `NULL` uses
  `outcome_model`. For count outcomes the exposure offset is taken from
  `outcome_model`; an offset in `init.model`, if supplied, must match
  it. Other score methods do not use this argument. Continuous-outcome
  fitting and validation are not currently implemented.

- initial.predictor.method:

  A character vector for the method used to get initial outcome
  predictions conditional on the covariates specified in `init.model`
  (or `outcome_model` when omitted). Only applies when `score.method`
  includes `'twoReg'` or `'contrastReg'`.Allowed values include one of
  `'randomForest'` (survival outcomes only), `'boosting'`, `'logistic'`
  (survival outcomes only, fast), `'poisson'` (count outcomes only,
  fast), `'gaussian'` (continuous outcomes only) and `'gam'` (count and
  continuous outcomes only). Default is `NULL`, which assigns
  `'boosting'` for count outcomes and `'randomForest'` for survival
  outcomes.

- ipcw.model:

  A formula describing the censoring model used for inverse probability
  of censoring weighting (IPCW). The left-hand side must be empty. Only
  applies for survival outcomes. Default is `NULL`, which uses a
  censoring model with the same covariates as the outcome model
  `outcome_model`, plus the treatment.

- ipcw.method:

  A character value for the censoring model. Only applies for survival
  outcomes. Allowed values are: `'breslow'` (Cox regression with Breslow
  estimator of t he baseline survivor function), `'aft (exponential)'`,
  `'aft (weibull)'`, `'aft (lognormal)'` or `'aft (loglogistic)'`
  (accelerated failure time model with different distributions for y
  variable). Default is `'breslow'`.

- minPS:

  A numerical value (in \`\[0, 1\]\`) below which estimated propensity
  scores should be truncated. Default is `0.01`.

- maxPS:

  A numerical value (in \`(0, 1\]\`) above which estimated propensity
  scores should be truncated. Must be strictly greater than `minPS`.
  Default is `0.99`.

- followup.time:

  A column name in `data` specifying the maximum follow-up time,
  interpreted as the potential censoring time. Only applies for survival
  outcomes. Default is `NULL`, which corresponds to unknown potential
  censoring time.

- tau0:

  The truncation time for defining restricted mean time lost. Only
  applies for survival outcomes. Default is `NULL`, which corresponds to
  setting the truncation time as the maximum survival time in the data.

- higher.y:

  A logical value indicating whether higher (`TRUE`) or lower (`FALSE`)
  values of the outcome are more desirable. Default is `TRUE`.

- prop.cutoff:

  A vector of numerical values (in \`(0, 1\]\`) specifying percentiles
  of the estimated log CATE scores to define nested subgroups. Each
  element represents the cutoff to separate observations in nested
  subgroups (below vs above cutoff). The length of `prop.cutoff` is the
  number of nested subgroups. An equally-spaced sequence of proportions
  ending with 1 is recommended. Default is `seq(0.5, 1, length = 6)`.

- surv.min:

  Lower truncation limit for the probability of being censored. It must
  be a positive value and should be chosen close to 0. Only applies for
  survival outcomes. Default is `0.025`.

- xvar.smooth.score:

  A vector of characters indicating the name of the variables used as
  the smooth terms if `score.method = 'gam'`. The variables must be
  selected from the variables listed in `outcome_model`.

- xvar.smooth.init:

  A vector of characters indicating the name of the variables used as
  the smooth terms if `initial.predictor.method = 'gam'`. The variables
  must be selected from the variables listed in `init.model`. Default is
  `NULL`, which uses all variables in `init.model`.

- tree.depth:

  A positive integer specifying the depth of individual trees in
  boosting (usually 2-3). Used only if `score.method = 'boosting'` or if
  `initial.predictor.method = 'boosting'` with `score.method = 'twoReg'`
  or `'contrastReg'`. Default is 2.

- n.trees.rf:

  A positive integer specifying the number of trees in random forest.
  Used if `score.method = 'randomForest'` or if
  `initial.predictor.method = 'randomForest'` with
  `score.method = 'twoReg'` or `'contrastReg'`. Only applies for
  survival outcomes. Default is `1000`.

- n.trees.boosting:

  A positive integer specifying the maximum number of trees in boosting
  (usually 100-1000). Used if `score.method = 'boosting'` or if
  `initial.predictor.method = 'boosting'` with `score.method = 'twoReg'`
  or `'contrastReg'`. Default is `200`.

- B:

  A positive integer specifying the number of time cross-fitting is
  repeated in `score.method = 'twoReg'` and `'contrastReg'`. Default is
  `3`.

- Kfold:

  A positive integer specifying the number of folds used in
  cross-fitting to partition the data in `score.method = 'twoReg'` and
  `'contrastReg'`. Default is `5`.

- error.maxNR:

  A numerical value \> 0 indicating the minimum value of the mean
  absolute error in Newton Raphson algorithm. Used only if
  `score.method = 'contrastReg'`. Default is `0.001`.

- max.iterNR:

  A positive integer indicating the maximum number of iterations in the
  Newton Raphson algorithm. Used only if `score.method = 'contrastReg'`.
  Default is `150`.

- tune:

  A vector of 2 numerical values \> 0 specifying tuning parameters for
  the Newton Raphson algorithm. `tune[1]` is the step size, `tune[2]`
  specifies a quantity to be added to diagonal of the slope matrix to
  prevent singularity. Used only if `score.method = 'contrastReg'`.
  Default is `c(0.5, 2)`.

- seed:

  An optional integer specifying an initial randomization seed for
  reproducibility. Default is `NULL`, corresponding to no seed.

- plot.gbmperf:

  A logical value indicating whether to plot the performance measures in
  boosting. Used only if `score.method = 'boosting'` or if
  `score.method = 'twoReg'` or `'contrastReg'` and
  `initial.predictor.method = 'boosting'`. Default is `FALSE`.

- verbose:

  An integer value indicating whether intermediate progress messages and
  histograms should be printed. `1` indicates messages are printed and
  `0` otherwise. Default is `0`.

- contrast.control:

  A named list of options for count or survival contrast regression.
  Supported names are `solver` (`"C"`, the default, or `"R"`),
  `error.maxNR`, `max.iterNR`, `tune`, and `optim.control` (a named list
  passed to [`stats::optim`](https://rdrr.io/r/stats/optim.html) on
  fallback). List entries override the corresponding top-level NR
  arguments. Options are kept separate from `...`, which is forwarded to
  boosting. Count and survival outcomes support this list; continuous
  outcomes retain their existing solver.

- rmst.control:

  Named options for Cox restricted mean survival integration. Supported
  entry: `solver = "C"` (default) or `solver = "R"`. The C backend
  accumulates each patient's integral without allocating a
  patient-by-event-time matrix. Cox model fitting is unchanged.
  Available only for survival outcomes; independent of
  `contrast.control`.

- ps.control:

  Named options for LASSO propensity-score interaction screening.
  Supported entry: `solver = "C"` (default) or `solver = "R"`. Only
  affects screening when `ps.method = "lasso"`; the LASSO fit and
  interaction selection rule are unchanged. Kept separate from
  boosting's `...` and other numerical controls.

- rf.control:

  Named list of randomForestSRC options: `mtry`, `nodesize`, `nsplit`,
  `nodedepth`, and `splitrule`. Applies to both treatment arms and
  survival initial-predictor forests. Use `n.trees.rf` for the number of
  trees.

- ...:

  Additional boosting arguments; see
  [`boosting-control`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/boosting-control.md).
  The old named arguments `score.method` and `cate.model` are accepted
  with deprecation warnings. Legacy score names and positional calls
  remain supported when `outcome_method` is omitted.

## Value

Fitted objects retain `treatment.labels` in treatment 0/1 order for
score plot labels. For binary response, see
[`catefitbinary()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitbinary.md).
For count response, see description of other outputs in
[`catefitcount()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitcount.md).
For survival response, see description of outputs in
[`catefitsurv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitsurv.md).

## Details

The table summarizes the available score-fitting methods. "Yes" in the
final column indicates a doubly robust score estimator, subject to the
assumptions described below; "–" means that the method is not available
for that outcome.

|                  |            |           |              |                   |
|------------------|------------|-----------|--------------|-------------------|
| **Score method** | **Binary** | **Count** | **Survival** | **Doubly robust** |
| `twoReg`         | Yes        | Yes       | Yes          | Yes               |
| `contrastReg`    | –          | Yes       | Yes          | Yes               |
| `logistic`       | Yes        | –         | –            | No                |
| `poisson`        | –          | Yes       | Yes          | No                |
| `negBin`         | –          | Yes       | –            | No                |
| `boosting`       | Yes        | Yes       | Yes          | No                |
| `randomForest`   | –          | –         | Yes          | No                |

The doubly robust methods combine outcome and treatment models. Their
robustness to misspecification of one of these models is conditional on
identification, positivity, and method-specific model and regularity
assumptions. For example, binary `twoReg` estimates a linear projection
of the conditional risk difference; this projection equals the
conditional risk difference only if the latter is linear in the
specified design. Truncating propensity scores can affect the robustness
property. Survival methods also require appropriate censoring
assumptions and a correctly specified censoring model.

This table describes *score fitting*. Subgroup effects can be estimated
using doubly robust methods even when the score used to define the
subgroups is not doubly robust. Neither score fitting nor subgroup
estimates calculated on the fitting data constitute independent
validation; see
[`catecv`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md).

For binary response, see the Binary outcomes section below and
[`catefitbinary()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitbinary.md).
For count response, see details in
[`catefitcount()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitcount.md).
For survival response, see details in
[`catefitsurv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefitsurv.md).

## Binary outcomes

For `response = "binary"`, `score.method = "logistic"` fits a separate
logistic outcome model in each treatment arm and subtracts their
predicted risks. `"twoReg"` cross-fits separate logistic outcome models
and a treatment model, then projects the augmented inverse-propensity
outcomes for each arm onto the covariates in `outcome_model` using
linear regression. The score is the difference of these projections,
averaged over `B` repetitions of `Kfold` folds. It estimates the linear
projection of the conditional risk difference if either the treatment
model or both outcome models are correct, under identification and
positivity assumptions. Propensity score truncation can affect this
property. The projection need not lie between -1 and 1 and equals the
conditional risk difference only if that difference is linear in the
specified design. The difference between arm coefficients is returned in
`coefficients$twoReg` and by
[`coef()`](https://rdrr.io/r/stats/coef.html); the separate arm
coefficients are stored in `fit$twoReg$coefficients`. Logistic model
coefficients for each arm are on the log-odds scale and stored in
`fit$logistic`; they are not linear risk-difference score coefficients.
Boosting does not have a linear score coefficient vector. With neither
`twoReg` nor another linear binary score,
[`coef()`](https://rdrr.io/r/stats/coef.html) returns `NULL`.
`"boosting"` fits separate Bernoulli gradient-boosting models in each
treatment arm, with the number of trees selected by five-fold
cross-validation. All three scores are predicted risk differences;
subgroup effects are doubly robust risk differences. Binary boosting
requires at least one outcome-model covariate and both outcome values in
each arm. Its fitted models and selected tree counts are stored in
`fit$boosting`. Scores and subgroup effects are estimated on the same
data; they are not independent validation. For extra boosting controls,
see
[`boosting-control`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/boosting-control.md).
Subgroup fits without both arms or with rank-deficient models return NA.

## References

Yadlowsky, S., Pellegrini, F., Lionetto, F., Braune, S., & Tian, L.
(2020). *Estimation and validation of ratio-based conditional average
treatment effects using observational data. Journal of the American
Statistical Association, 1-18.* DOI: 10.1080/01621459.2020.1772080.

Hoogland, J., Debray, T. P. A., Desachy, G., & Nguyen, T.-L. (2026).
Introduction to methods for personalizing medicine. In T. P. A. Debray,
T.-L. Nguyen, & R. W. Platt (Eds.), *Comparative Effectiveness and
Personalized Medicine Research Using Real-World Data* (Chapter 17, pp.
421–463). Chapman & Hall/CRC.
[doi:10.1201/9781003300809-17](https://doi.org/10.1201/9781003300809-17)

Hoogland, J., IntHout, J., Belias, M., et al. (2021). A tutorial on
individualized treatment effect prediction from randomized trials with a
binary endpoint. *Statistics in Medicine*, 40(26), 5961–5981.
[doi:10.1002/sim.9154](https://doi.org/10.1002/sim.9154)

Jiang, X., Simoneau, G., Pellegrini, F., Shen, C., de Moor, C., & Tian,
L. (2026). Visualization and interpretation of individualized treatment
rule results. In T. P. A. Debray, T.-L. Nguyen, & R. W. Platt (Eds.),
*Comparative Effectiveness and Personalized Medicine Research Using
Real-World Data* (Chapter 20, pp. 513–536). Chapman & Hall/CRC.
[doi:10.1201/9781003300809-20](https://doi.org/10.1201/9781003300809-20)

## See also

[`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)

## Examples

``` r
# Count outcome
fit_1 <- catefit(response = "count",
                 data = countExample,
                 outcome_model = y ~ age + female + previous_treatment +
                                  previous_cost + previous_number_relapses +
                                  offset(log(years)),
                 outcome_method = "poisson",
                 benefit_method = "t_learner",
                 ps.model = trt ~ age + previous_treatment,
                 higher.y = TRUE,
                 seed = 999)
#> Variable trt was recoded to 0/1 with drug0->0 and drug1->1.

coef(fit_1)
#>                              poisson
#> (Intercept)              -0.58185218
#> age                      -0.01117701
#> female                    0.58829733
#> previous_treatmentdrugB   0.73724361
#> previous_treatmentdrugC   0.15545985
#> previous_cost            -0.13035642
#> previous_number_relapses -0.05339112

# \donttest{
# Survival outcome
library(survival)
tau0 <- with(survivalExample,
                 min(quantile(y[trt == "drug1"], 0.95), quantile(y[trt == "drug0"], 0.95)))

fit_2 <- catefit(response = "survival",
                 data = survivalExample,
                 outcome_model = Surv(y, d) ~ age + female + previous_cost +
                                           previous_number_relapses,
                 outcome_method = c("poisson", "boosting", "random_forest"),
                 benefit_method = rep("t_learner", 3),
                 ps.model = trt ~ age + previous_treatment,
                 initial.predictor.method = "logistic",
                 ipcw.model = ~ age + previous_cost + previous_treatment,
                 tau0 = tau0, higher.y = TRUE, seed = 999, n.cores = 1)
#> Variable trt was recoded to 0/1 with drug0->0 and drug1->1.
#>     Warning: Warning(s) occurred when estimating the ATEs in the nested subgroup using "poisson" 
#> Warning: Warning: Warning(s) occurred when estimating the ATEs in the nested subgroup using "poisson"; see 'errors/warnings'.

coef(fit_2)
#>                              poisson
#> (Intercept)               0.46109243
#> age                      -0.24335209
#> female                    0.12532711
#> previous_cost             0.89803249
#> previous_number_relapses -0.02470863

# }
```
