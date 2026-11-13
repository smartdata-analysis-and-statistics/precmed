# Calculate the CATE score given the baseline covariates for specified scoring method methods

Based on intxmean results of the CATE coefficients estimated with
boosting, linear regression, two regression, contrast regression, random
forest, generalized additive model

## Usage

``` r
scoremean(
  fit,
  x.cate,
  score.method = c("boosting", "gaussian", "twoReg", "contrastReg", "randomForest",
    "gam")
)
```

## Arguments

- fit:

  List of objects generated from intxmean: outputs of boosting, linear
  regression, two regression, contrast regression, random forest,
  generalized additive model

- x.cate:

  Matrix of `p.cate` baseline covariates; dimension `n` (observations)
  by `p.cate` (covariates in the outcome model)

- score.method:

  A vector of one or multiple methods to estimate the CATE score.
  Allowed values are: `'boosting'`, `'gaussian'`, `'twoReg'`,
  `'contrastReg'`, `'randomForest'`, `'gam'`. Default specifies all 6
  methods.

## Value

score.boosting: Estimated CATE score for all `n` observations with the
boosting method; vector of size `n` score.gaussian: Estimated CATE score
for all `n` observations with the linear regression method; vector of
size `n` score.twoReg: Estimated CATE score for all `n` observations
with the two regression method; vector of size `n` score.contrastReg:
Estimated CATE score for all `n` observations with the contrast
regression method; vector of size `n` score.randomForest: Estimated CATE
score for all `n` observations with the random forest method; vector of
size `n` score.gam: Estimated CATE score for all `n` observations with
the generalized additive model; vector of size `n` score = NA if the
corresponding method is not called
