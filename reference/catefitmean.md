# Estimation of the conditional average treatment effect (CATE) score for continuous data

Provides singly robust and doubly robust estimation of CATE score with
up to 6 scoring methods among the following: Linear regression,
boosting, two regressions, contrast regression, random forest and
generalized additive model.

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

  Formula for the outcome and outcome-model covariates; replaces
  `cate.model`.

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
  the smooth terms if `score.method = 'gam'`. The variables must be
  selected from the variables listed in `cate.model`. Default is `NULL`,
  which uses all variables in `cate.model`.

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
  `score.method = 'randomForest'`. Default is `1000`.

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

  The legacy `init.model` argument specifies a separate initial outcome
  formula under the previous interface; it cannot be combined with
  `benefit_model`. The legacy `initial.predictor.method` is also
  accepted with a deprecation warning; do not supply it with
  `outcome_method`. Legacy arguments `ps.model` and `ps.method` are
  accepted with warnings. Additional boosting arguments; see
  [`boosting-control`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/boosting-control.md).

## Value

Returns a list containing the following components:

- `ate.gaussian`: A vector of numerical values of length `prop.cutoff`
  containing the estimated ATE in nested subgroups (defined by
  `prop.cutoff`) constructed based on the estimated CATE scores with
  Poisson regression. Only provided if `score.method` includes
  `'gaussian'`.

- `ate.boosting`: Same as `ate.gaussian`, but with the nested subgroups
  based the estimated CATE scores with boosting. Only provided if
  `score.method` includes `'boosting'`.

- `ate.twoReg`: Same as `ate.gaussian`, but with the nested subgroups
  based the estimated CATE scores with two regressions. Only provided if
  `score.method` includes `'twoReg'`.

- `ate.contrastReg`: Same as `ate.gaussian`, but with the nested
  subgroups based the estimated CATE scores with contrast regression.
  Only provided if `score.method` includes `'contrastReg'`.

- `ate.randomForest`: Same as `ate.gaussian`, but with the nested
  subgroups based the estimated CATE scores with random forest. Only
  provided if `score.method` includes `'gam'`.

- `ate.gam`: Same as `ate.gaussian`, but with the nested subgroups based
  the estimated CATE scores with generalized additive model. Only
  provided if `score.method` includes `'gam'`.

- `score.gaussian`: A vector of numerical values of length n (number of
  observations in `data`) containing the estimated CATE scores according
  to the linear regression. Only provided if `score.method` includes
  `'gaussian'`.

- `score.boosting`: Same as `score.gaussian`, but with estimated CATE
  score according to boosting. Only provided if `score.method` includes
  `'boosting'`.

- `score.twoReg`: Same as `score.gaussian`, but with estimated CATE
  score according to two regressions. Only provided if `score.method`
  includes `'twoReg'`.

- `score.contrastReg`: Same as `score.gaussian`, but with estimated CATE
  score according to contrast regression. Only provided if
  `score.method` includes `'contrastReg'`.

- `score.randomForest`: Same as `score.gaussian`, but with estimated
  CATE score according to random forest. Only provided if `score.method`
  includes `'randomForest'`.

- `score.gam`: Same as `score.gaussian`, but with estimated CATE score
  according to generalized additive model. Only provided if
  `score.method` includes `'gam'`.

- `fit`: Additional details on model fitting if `score.method` includes
  'boosting' or 'contrastReg':

  - `result.boosting`: Details on the boosting model fitted to
    observations with treatment = 0 `($fit0.boosting)` and to
    observations with treatment = 1 `($fit1.boosting)`. Only provided if
    `score.method` includes `'boosting'`.

  - `result.randomForest`: Details on the boosting model fitted to
    observations with treatment = 0 `($fit0.randomForest)` and to
    observations with treatment = 1 `($fit1.randomForest)`. Only
    provided if `score.method` includes `'randomForest'`.

  - `result.gam`: Details on the boosting model fitted to observations
    with treatment = 0 `($fit0.gam)` and to observations with treatment
    = 1 `($fit1.gam)`. Only provided if `score.method` includes `'gam'`.

  - `result.contrastReg$sigma.contrastReg`: Variance-covariance matrix
    of the estimated CATE coefficients in contrast regression. Only
    provided if `score.method` includes `'contrastReg'`.

- `coefficients`: A data frame with the coefficients of the estimated
  CATE score by `score.method`. The data frame has number of rows equal
  to the number of covariates in `cate.model` and number of columns
  equal to `length(score.method)`. If `score.method` includes
  `'contrastReg'`, the data frame has an additional column containing
  the standard errors of the coefficients estimated with contrast
  regression. `'boosting'`, `'randomForest'`, `'gam'` do not have
  coefficient results because these methods do not express the CATE as a
  linear combination of coefficients and covariates.

## Details

The legacy arguments `score.method` and `cate.model` are accepted
through `...`; see
[`catefit`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md).

The CATE score estimates a covariate-conditional average treatment
effect, estimated with either linear regression, boosting, random forest
and generalized additive model applied separately by treatment group or
with two doubly robust estimators, two regressions and contrast
regression (Yadlowsky, 2020) applied to the entire dataset.

`catefitmean()` provides the coefficients of the CATE score for each
scoring method requested through `score.method`. Currently, contrast
regression is the only method which allows for inference of the CATE
coefficients by providing standard errors of the coefficients. The
coefficients can be used to learn the effect size of each variable and
predict the CATE score for a new observation.

`catefitmean()` also provides the predicted CATE score of each
observation in the data set, for each scoring method. The predictions
allow ranking the observations from potentially high responders to the
treatment to potentially low or standard responders.

The estimated ATE among nested subgroups of high responders are also
provided by scoring method. Note that the ATEs in `catefitmean()` are
derived based on the CATE score which is estimated using the same data
sample. Therefore, overfitting may be an issue. `catefitmean()` is more
suitable to inspect the estimated ATEs across scoring methods as it
implements internal cross validation to reduce optimism.

## References

Yadlowsky, S., Pellegrini, F., Lionetto, F., Braune, S., & Tian, L.
(2020). *Estimation and validation of ratio-based conditional average
treatment effects using observational data. Journal of the American
Statistical Association, 1-18.* DOI: 10.1080/01621459.2020.1772080.

## See also

[`catecvmean()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecvmean.md)
function
