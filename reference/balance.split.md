# Split the given dataset into balanced training and validation sets (within a pre-specified tolerance) Balanced means 1) The ratio of treated and controls is maintained in the training and validation sets 2) The covariate distributions are balanced between the training and validation sets

Split the given dataset into balanced training and validation sets
(within a pre-specified tolerance) Balanced means 1) The ratio of
treated and controls is maintained in the training and validation sets
2) The covariate distributions are balanced between the training and
validation sets

## Usage

``` r
balance.split(
  y,
  trt,
  x.cate,
  x.ps,
  time,
  minPS = 0.01,
  maxPS = 0.99,
  train.prop = 3/4,
  error.max = 0.1,
  max.iter = 5000
)
```

## Arguments

- y:

  Observed outcome; vector of size `n` (observations)

- trt:

  Treatment received; vector of size `n` with treatment coded as 0/1

- x.cate:

  Matrix of `p.cate` baseline covariates; dimension `n` by `p.cate`
  (covariates in the outcome model)

- x.ps:

  Matrix of `p.ps` baseline covariates (plus a leading column of 1 for
  the intercept); dimension `n` by `p.ps + 1` (covariates in the
  propensity score model plus intercept)

- time:

  Log-transformed person-years of follow-up; vector of size `n`

- minPS:

  A numerical value (in \`\[0, 1\]\`) below which estimated propensity
  scores should be truncated. Default is `0.01`.

- maxPS:

  A numerical value (in \`(0, 1\]\`) above which estimated propensity
  scores should be truncated. Must be strictly greater than `minPS`.
  Default is `0.99`.

- train.prop:

  A numerical value (in \`(0, 1)\`) indicating the proportion of total
  data used for training. Default is `3/4`.

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

## Value

A list of 10 objects, 5 training and 5 validation of y, trt, x.cate,
x.ps, time: y.train - observed outcome in the training set; vector of
size `m` (observations in the training set) trt.train - treatment
received in the training set; vector of size `m` coded as 0/1
x.cate.train - baseline covariates for the outcome model in the training
set; matrix of dimension `m` by `p.cate` x.ps.train - baseline
covariates (plus intercept) for the propensity score model in the
training set; matrix of dimension `m` by `p.ps + 1` time.train -
log-transformed person-years of follow-up in the training set; vector of
size `m` y.valid - observed outcome in the validation set; vector of
size `n-m` trt.valid - treatment received in the validation set; vector
of size `n-m` coded as 0/1 x.cate.valid - baseline covariates for the
outcome model in the validation set; matrix of dimension `n-m` by
`p.cate` x.ps.valid - baseline covariates (plus intercept) for the
propensity score model in the validation set; matrix of dimension `n-m`
by `p.ps + 1` time.valid - log-transformed person-years of follow-up in
the validation set; vector of size `n-m`
