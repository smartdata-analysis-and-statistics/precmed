# Estimate conditional treatment-effect scores for binary outcomes

Fit binary-outcome CATE scores and estimate risk differences in nested
subgroups. This is the binary-outcome interface used by
[`catefit`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md).

## Usage

``` r
catefitbinary(
  data,
  outcome_model = NULL,
  outcome_method = NULL,
  benefit_method = NULL,
  ps_model,
  ps_method = "glm",
  minPS = 0.01,
  maxPS = 0.99,
  higher.y = TRUE,
  prop.cutoff = seq(0.5, 1, length = 6),
  seed = NULL,
  ps.control = list(),
  B = 3,
  Kfold = 5,
  tree.depth = 2,
  n.trees.boosting = 200,
  plot.gbmperf = FALSE,
  init.model = NULL,
  ...
)
```

## Arguments

- data:

  A data frame containing the variables in the outcome, treatment, and
  censoring models (if specified); a data frame with `n` rows (1 row per
  observation).

- outcome_model:

  Outcome and score covariate formula; replaces `cate.model`.

- outcome_method:

  Outcome-model methods paired with `benefit_method`; see
  [`catefit`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md).

- benefit_method:

  Score-construction methods; see
  [`catefit`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md).

- ps_model:

  A formula describing the treatment model (propensity-score model) to
  be fitted. The treatment must appear on the left-hand side. For binary
  outcomes, logical treatment values `FALSE` and `TRUE` are mapped to 0
  and 1, respectively. Numeric 0/1 values or two distinct treatment
  labels are also accepted. If data are from a randomized controlled
  trial, specify `ps_model = trt ~ 1` as an intercept-only model.

- ps_method:

  A character value for the method to estimate the propensity score.
  Allowed values include one of: `'glm'` for logistic regression with
  main effects only (default), or `'lasso'` for a logistic regression
  with main effects and LASSO penalization on two-way interactions
  (added to the model if interactions are not specified in `ps_model`).
  Relevant only when `ps_model` has more than one variable.

- minPS:

  A numerical value (in \`\[0, 1\]\`) below which estimated propensity
  scores should be truncated. Default is `0.01`.

- maxPS:

  A numerical value (in \`(0, 1\]\`) above which estimated propensity
  scores should be truncated. Must be strictly greater than `minPS`.
  Default is `0.99`.

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

- seed:

  An optional integer specifying an initial randomization seed for
  reproducibility. Default is `NULL`, corresponding to no seed.

- ps.control:

  Named options for LASSO propensity-score interaction screening.
  Supported entry: `solver = "C"` (default) or `solver = "R"`. Only
  affects screening when `ps_method = "lasso"`; the LASSO fit and
  interaction selection rule are unchanged. Kept separate from
  boosting's `...` and other numerical controls.

- B:

  A positive integer specifying the number of time cross-fitting is
  repeated in `score.method = 'twoReg'` and `'contrastReg'`. Default is
  `3`.

- Kfold:

  A positive integer specifying the number of folds used in
  cross-fitting to partition the data in `score.method = 'twoReg'` and
  `'contrastReg'`. Default is `5`.

- tree.depth:

  A positive integer specifying the depth of individual trees in
  boosting (usually 2-3). Used only if `score.method = 'boosting'` or if
  `initial.predictor.method = 'boosting'` with `score.method = 'twoReg'`
  or `'contrastReg'`. Default is 2.

- n.trees.boosting:

  A positive integer specifying the maximum number of trees in boosting
  (usually 100-1000). Used if `score.method = 'boosting'` or if
  `initial.predictor.method = 'boosting'` with `score.method = 'twoReg'`
  or `'contrastReg'`. Default is `200`.

- plot.gbmperf:

  A logical value indicating whether to plot the performance measures in
  boosting. Used only if `score.method = 'boosting'` or if
  `score.method = 'twoReg'` or `'contrastReg'` and
  `initial.predictor.method = 'boosting'`. Default is `FALSE`.

- init.model:

  Optional formula for the initial logistic outcome models used by
  `twoReg`; `NULL` uses `outcome_model`.

- ...:

  Legacy arguments `ps.model` and `ps.method` are accepted with
  warnings. Additional boosting arguments; see
  [`boosting-control`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/boosting-control.md).

## Value

A `catefit` object containing treatment labels, fitted models, predicted
risk-difference scores and nested subgroup risk differences for each
requested method. Linear `twoReg` score coefficients are stored in
`coefficients`.

## Details

The legacy arguments `score.method` and `cate.model` are accepted
through `...`; see
[`catefit`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md).

See the Binary outcomes section of
[`catefit`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
for the score methods, their assumptions, and interpretation. Subgroup
effects are estimated on the fitting data and do not provide independent
validation.

## Author

Thomas Debray <tdebray@fromdatatowisdom.com>

## Examples

``` r
# \donttest{
dat <- simulate_binary(n = 400, seed = 2026)
fit <- catefitbinary(data = dat,
                     outcome_model = y ~ age + previous_treatment,
                     outcome_method = "logistic",
                     benefit_method = "t_learner",
                     ps_model = trt ~ age + previous_treatment,
                     seed = 2026)
# }
```
