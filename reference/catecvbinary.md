# Cross-validation of conditional risk-difference scores

Repeatedly fits binary treatment-effect scores on treatment-stratified
training splits and applies them to the held-out observations.

## Usage

``` r
catecvbinary(
  data,
  outcome_model,
  outcome_method = "logistic",
  benefit_model = NULL,
  benefit_method = "t_learner",
  ps_model,
  ps_method = "glm",
  minPS = 0.01,
  maxPS = 0.99,
  higher.y = TRUE,
  prop.cutoff = seq(0.5, 1, length = 6),
  prop.multi = c(0, 1/3, 2/3, 1),
  abc = TRUE,
  train.prop = 3/4,
  cv.n = 10,
  tree.depth = 2,
  n.trees.boosting = 200,
  B = 3,
  Kfold = 5,
  seed = NULL,
  plot.gbmperf = FALSE,
  verbose = 0,
  ps.control = list(),
  ...
)
```

## Arguments

- data:

  A data frame containing the variables in the outcome, treatment, and
  censoring models (if specified); a data frame with `n` rows (1 row per
  observation).

- outcome_model:

  A formula specifying the outcome and baseline covariates. Use a `Surv`
  response for survival and include an exposure offset for count
  outcomes.

- outcome_method:

  Outcome-model method for each `benefit_method`. Choose `"logistic"`
  (default) or `"boosting"` for a T-learner. Two-regression uses
  `"logistic"`.

- benefit_model:

  Optional one-sided formula for the covariates and terms used in
  two-regression and contrast-regression scores. The default uses the
  terms in `outcome_model`. T-learners use `outcome_model` directly.

- benefit_method:

  Score-construction method, with one entry per `outcome_method`. The
  default is `'t_learner'`. See
  [`catefit`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
  for supported combinations.

- ps_model:

  A formula for the treatment model, with treatment on the left. It
  estimates propensity scores. For a common randomization probability,
  use an intercept-only formula such as `trt ~ 1`.

- ps_method:

  Treatment-model fitting method. The default `'glm'` uses logistic
  regression. `'lasso'` uses a logistic model with LASSO screening for
  two-way interactions.

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
  of the estimated CATE scores to define nested subgroups. Each element
  represents the cutoff to separate observations in nested subgroups
  (below vs above cutoff). The length of `prop.cutoff` is the number of
  nested subgroups. An equally-spaced sequence of proportions ending
  with 1 is recommended. Default is `seq(0.5, 1, length = 6)`.

- prop.multi:

  A vector of numerical values (in \`\[0, 1\]\`) specifying percentiles
  of the estimated CATE scores to define mutually exclusive subgroups.
  It should start with 0, end with 1, and be of
  `length(prop.multi) > 2`. Each element represents the cutoff to
  separate the observations into `length(prop.multi) - 1` mutually
  exclusive subgroups. Default is `c(0, 1/3, 2/3, 1)`.

- abc:

  A logical value indicating whether the area between curves (ABC)
  should be calculated at each cross-validation iterations, for each
  scoring method. Default is `TRUE`. Binary ABC requires at least two
  distinct `prop.cutoff` values.

- train.prop:

  A numerical value (in \`(0, 1)\`) indicating the proportion of total
  data used for training. Default is `3/4`.

- cv.n:

  A positive integer value indicating the number of cross-validation
  iterations. Default is `10`.

- tree.depth:

  Maximum depth of boosting trees. Defaults to 2.

- n.trees.boosting:

  Maximum number of trees for boosting. Defaults to 200.

- B:

  Number of cross-fitting repetitions used by two-regression and
  contrast regression. Defaults to 3.

- Kfold:

  Number of folds used for cross-fitting within each training set.
  Defaults to 5.

- seed:

  An optional integer specifying an initial randomization seed for
  reproducibility. Default is `NULL`, corresponding to no seed.

- plot.gbmperf:

  Whether to plot boosting performance. Defaults to FALSE.

- verbose:

  An integer value indicating what kind of intermediate progress
  messages should be printed. `0` means no outputs. `1` means only
  progress bar and run time. `2` means progress bar, run time, and all
  errors and warnings. Default is `0`.

- ps.control:

  Named options for LASSO propensity-score interaction screening.
  Supported entry: `solver = "C"` (default) or `solver = "R"`. Only
  affects screening when `ps_method = "lasso"`; the LASSO fit and
  interaction selection rule are unchanged. Kept separate from
  boosting's `...` and other numerical controls.

- ...:

  Additional arguments. Deprecated argument names from precmed 1.1.0
  remain supported for backward compatibility and emit warnings.

## Value

A `"precmed"` object with training and validation subgroup risk
differences, split indices, and diagnostics. See
[`catecv`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md).

## Author

Thomas Debray \<tdebray@fromdatatowisdom.com\>
