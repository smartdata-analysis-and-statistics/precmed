# Calculate the log CATE score given the baseline covariates and follow-up time for specified scoring method methods

Based on intxcount results of the CATE coefficients estimated with
boosting, naive Poisson, two regression, contrast regression, negative
binomial

## Usage

``` r
scorecount(
  fit,
  x.cate,
  time,
  score.method = c("boosting", "poisson", "twoReg", "contrastReg", "negBin")
)
```

## Arguments

- fit:

  List of objects generated from intxcount: outputs of boosting, naive
  Poisson, two regression, contrast regression, negative binomial

- x.cate:

  Matrix of `p.cate` baseline covariates; dimension `n` (observations)
  by `p.cate` (covariates in the outcome model)

- time:

  Log-transformed person-years of follow-up; vector of size `n`

- score.method:

  A vector of one or multiple methods to estimate the CATE score.
  Allowed values are: `'boosting'`, `'poisson'`, `'twoReg'`,
  `'contrastReg'`, `'negBin'`. Default specifies all 5 methods.

## Value

score.boosting: Estimated log CATE score for all `n` observations with
the boosting method; vector of size `n` score.poisson: Estimated log
CATE score for all `n` observations with the naive Poisson method;
vector of size `n` score.twoReg: Estimated log CATE score for all `n`
observations with the two regression method; vector of size `n`
score.contrastReg: Estimated log CATE score for all `n` observations
with the contrast regression method; vector of size `n` score.negBin:
Estimated log CATE score for all `n` observations with the naive Poisson
method; vector of size `n` score = NA if the corresponding method is not
called
