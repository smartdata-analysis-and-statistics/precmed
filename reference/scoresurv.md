# Calculate the log CATE score given the baseline covariates and follow-up time for specified scoring method methods for survival outcomes

Based on intxsurv results of the CATE coefficients estimated with random
forest, boosting, naive Poisson, two regression, contrast regression

## Usage

``` r
scoresurv(
  fit,
  x.cate,
  tau0,
  score.method = c("randomForest", "boosting", "poisson", "twoReg", "contrastReg")
)
```

## Arguments

- fit:

  List of objects generated from intxsurv: outputs of random forest,
  boosting, naive Poisson, two regression, contrast regression

- x.cate:

  Matrix of `p.cate` baseline covariates specified in the outcome model;
  dimension `n` by `p.cate`.

- tau0:

  The truncation time for defining restricted mean time lost.

- score.method:

  A vector of one or multiple methods to estimate the CATE score.
  Allowed values are: `'randomForest'`, `'boosting'`, `'poisson'`,
  `'twoReg'`, `'contrastReg'`. Default specifies all 5 methods.

## Value

score.randomForest: Estimated log CATE score for all `n` observations
with the random forest method; vector of size `n` score.boosting:
Estimated log CATE score for all `n` observations with the boosting
method; vector of size `n` score.poisson: Estimated log CATE score for
all `n` observations with the naive Poisson method; vector of size `n`
score.twoReg: Estimated log CATE score for all `n` observations with the
two regression method; vector of size `n` score.contrastReg: Estimated
log CATE score for all `n` observations with the contrast regression
method; vector of size `n` score = NA if the corresponding method is not
called
