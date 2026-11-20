# Cross-validation for continuous outcomes

Continuous-outcome cross-validation is not implemented.

## Usage

``` r
catecvmean(
  data,
  outcome_model,
  outcome_method = NULL,
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
  error.max = 0.1,
  max.iter = 5000,
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

  A formula specifying the outcome and baseline covariates. Use a `Surv`
  response for survival and include an exposure offset for count
  outcomes.

- outcome_method:

  Outcome-model method. Continuous-outcome cross-validation is not
  implemented.

- benefit_model:

  Optional one-sided formula for the covariates and terms used in
  two-regression and contrast-regression scores. The default uses the
  terms in `outcome_model`. T-learners use `outcome_model` directly.

- benefit_method:

  Score-construction method. Continuous-outcome cross-validation is not
  implemented.

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
  scoring method. Default is `TRUE`.

- train.prop:

  A numerical value (in \`(0, 1)\`) indicating the proportion of total
  data used for training. Default is `3/4`.

- cv.n:

  A positive integer value indicating the number of cross-validation
  iterations. Default is `10`.

- error.max:

  A numerical value \> 0 indicating the tolerance (maximum value of
  error) for the largest standardized absolute difference in the
  covariate distributions or in the doubly robust estimated rate ratios
  between the training and validation sets. This is used to define a
  balanced training-validation splitting. Default is `0.1`.

- max.iter:

  A positive integer value indicating the maximum number of iterations
  when searching for a balanced training-validation split. Default is
  `5,000`.

- xvar.smooth.score:

  Reserved for continuous-outcome scores.

- xvar.smooth.init:

  Reserved for continuous-outcome validation. Use `xvar.smooth` for
  count outcomes.

- tree.depth:

  Maximum depth of boosting trees. Defaults to 2.

- n.trees.rf:

  Number of trees for survival random forests. Defaults to 1000.

- n.trees.boosting:

  Maximum number of trees for boosting. Defaults to 200.

- B:

  Number of cross-fitting repetitions used by two-regression and
  contrast regression. Defaults to 3.

- Kfold:

  Number of folds used for cross-fitting within each training set.
  Defaults to 6.

- plot.gbmperf:

  Whether to plot boosting performance. Defaults to FALSE.

- error.maxNR:

  Newton–Raphson convergence tolerance for contrast regression. Defaults
  to 0.001.

- tune:

  Two positive Newton–Raphson tuning values for contrast regression.
  Defaults to `c(0.5, 2)`.

- seed:

  An optional integer specifying an initial randomization seed for
  reproducibility. Default is `NULL`, corresponding to no seed.

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

- rf.control:

  Named list of randomForestSRC options: `mtry`, `nodesize`, `nsplit`,
  `nodedepth`, and `splitrule`. Applies to both treatment arms and
  survival initial-predictor forests. Use `n.trees.rf` for the number of
  trees.

- ...:

  Additional arguments. Deprecated argument names from precmed 1.1.0
  remain supported for backward compatibility and emit warnings.

## Value

This function stops with an error.
