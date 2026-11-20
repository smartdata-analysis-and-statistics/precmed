# Estimation of the conditional average treatment effect (CATE) score for continuous data

Continuous-outcome CATE fitting is not implemented.

## Usage

``` r
catefitmean(
  data,
  outcome_model = NULL,
  outcome_method = NULL,
  benefit_model = NULL,
  benefit_method = NULL,
  ps_model,
  ps_method = "glm",
  minPS = 0.01,
  maxPS = 0.99,
  higher.y = TRUE,
  prop.cutoff = seq(0.5, 1, length = 6),
  xvar.smooth.score = NULL,
  xvar.smooth.init = NULL,
  tree.depth = 2,
  n.trees.rf = 1000,
  n.trees.boosting = 200,
  B = 3,
  Kfold = 6,
  plot.gbmperf = FALSE,
  error.maxNR = 0.001,
  tune = c(0.5, 2),
  seed = NULL,
  verbose = 0,
  ps.control = list(),
  rf.control = list(),
  ...
)
```

## Arguments

- data:

  A data frame containing the variables in the outcome and treatment
  models; a data frame with `n` rows (1 row per observation).

- outcome_model:

  A formula specifying the outcome and baseline covariates.

- outcome_method:

  Outcome-model methods paired with `benefit_method`; see
  [`catefit`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md).

- benefit_model:

  Optional one-sided formula for the two-regression or
  contrast-regression score covariates. Defaults to the covariates and
  terms in `outcome_model`; T-learner scores do not use it. See
  [`catefit`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md).

- benefit_method:

  Score-construction methods; see
  [`catefit`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md).

- ps_model:

  A formula for the treatment model, which estimates propensity scores.
  Put the treatment variable on the left-hand side; for an
  intercept-only treatment model, use for example `trt ~ 1`. See
  [`catefit`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md).

- ps_method:

  A character value for the method to estimate the propensity score.
  Allowed values include one of: `'glm'` for logistic regression with
  main effects only (default), or `'lasso'` for a logistic regression
  with main effects and LASSO penalization on two-way interactions
  (added to the model if not specified in `ps_model`).

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

- xvar.smooth.score:

  A vector of characters indicating the name of the variables used as
  the smooth terms if fitting a generalized additive model. The
  variables must be selected from the variables listed in
  `outcome_model`. Default is `NULL`, which uses all variables in
  `outcome_model`.

- xvar.smooth.init:

  A vector of characters indicating the name of the variables used as
  the smooth terms if `outcome_method = 'gam'`. The variables must be
  selected from the variables listed in `outcome_model`. Default is
  `NULL`, which uses all variables in `outcome_model`.

- tree.depth:

  A positive integer specifying the maximum depth of boosting trees.
  Used when `outcome_method = 'boosting'`. Default is `2`.

- n.trees.rf:

  A positive integer specifying the number of trees. Used only if
  `outcome_method = 'random_forest'`. Default is `1000`.

- n.trees.boosting:

  The maximum number of boosting trees, used when
  `outcome_method = 'boosting'`. Default is `200`.

- B:

  Number of repetitions of cross-fitting for
  `benefit_method = 'two_regression'` or `'contrast_regression'`.
  Default is `3`.

- Kfold:

  A positive integer specifying the number of folds (parts) used in
  cross-fitting to partition the data in
  `benefit_method = 'two_regression'` and `'contrast_regression'`.
  Default is `6`.

- plot.gbmperf:

  Whether to plot boosting performance measures when
  `outcome_method = 'boosting'`. Default is `FALSE`.

- error.maxNR:

  Convergence tolerance for the Newton-Raphson algorithm used by
  `benefit_method = 'contrast_regression'`. Default is `0.001`.

- tune:

  A vector of 2 numerical values \> 0 specifying tuning parameters for
  the Newton Raphson algorithm. `tune[1]` is the step size, `tune[2]`
  specifies a quantity to be added to diagonal of the slope matrix to
  prevent singularity. Used only if
  `benefit_method = 'contrast_regression'`. Default is `c(0.5, 2)`.

- seed:

  An optional integer specifying an initial randomization seed for
  reproducibility. Default is `NULL`, corresponding to no seed.

- verbose:

  An integer value indicating what kind of intermediate progress
  messages should be printed. `0` means no outputs. `1` means only
  progress and run time. `2` means progress, run time, and all errors
  and warnings. Default is `0`.

- ps.control:

  Named options for LASSO propensity-score interaction screening.
  Supported entry: `solver = "C"` (default) or `solver = "R"`. Only
  affects screening when `ps_method = "lasso"`; the LASSO fit and
  interaction selection rule are unchanged. Kept separate from
  boosting's `...` and other numerical controls.

- rf.control:

  Named list of randomForestSRC options: `mtry`, `nodesize`, `nsplit`,
  `nodedepth`, and `splitrule`. Applies to both treatment arms and
  survival initial-predictor forests. Use `n.trees.rf` for the number of
  trees.

- ...:

  Additional arguments. Deprecated argument names from precmed 1.1.0
  remain supported for backward compatibility and emit warnings.

## Value

Continuous-outcome fitting is not implemented. This function stops with
an error.
