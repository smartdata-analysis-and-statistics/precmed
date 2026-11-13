# Split the given time-to-event dataset into balanced training and validation sets (within a pre-specified tolerance) Balanced means 1) The ratio of treated and controls is maintained in the training and validation sets 2) The covariate distributions are balanced between the training and validation sets

Split the given time-to-event dataset into balanced training and
validation sets (within a pre-specified tolerance) Balanced means 1) The
ratio of treated and controls is maintained in the training and
validation sets 2) The covariate distributions are balanced between the
training and validation sets

## Usage

``` r
balancesurv.split(
  y,
  d,
  trt,
  x.cate,
  x.ps,
  x.ipcw,
  yf = NULL,
  train.prop = 3/4,
  error.max = 0.1,
  max.iter = 5000
)
```

## Arguments

- y:

  Observed survival or censoring time; vector of size `n`.

- d:

  The event indicator, normally `1 = event, 0 = censored`; vector of
  size `n`.

- trt:

  Treatment received; vector of size `n` with treatment coded as 0/1.

- x.cate:

  Matrix of `p.cate` baseline covariates specified in the outcome model;
  dimension `n` by `p.cate`.

- x.ps:

  Matrix of `p.ps` baseline covariates specified in the propensity score
  model; dimension `n` by `p.ps`.

- x.ipcw:

  Matrix of `p.ipw` baseline covariate specified in inverse probability
  of censoring weighting; dimension `n` by `p.ipw`.

- yf:

  Follow-up time, interpreted as the potential censoring time; vector of
  size `n` if the potential censoring time is known. If unknown, set
  `yf == NULL` and `yf` will be taken as `y` in the function.

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

A list of 14 objects, 7training and 7 validation of y, trt, x.cate,
x.ps, x.ipcw, time, yf: y.train - observed survival or censoring time in
the training set; vector of size `m` (observations in the training set)
d.train - event indicator in the training set; vector of size `m` coded
as 0/1 trt.train - treatment received in the training set; vector of
size `m` coded as 0/1 x.cate.train - baseline covariates for the outcome
model in the training set; matrix of dimension `m` by `p.cate`
x.ps.train - baseline covariates (plus intercept) for the propensity
score model in the training set; matrix of dimension `m` by `p.ps + 1`
x.ipcw.train - baseline covariates for inverse probability of censoring
in the training set; matrix of dimension `m` by `p.ipw` yf.train -
follow-up time in the training set; if known, vector of size `m`; if
unknown, `yf == NULL` y.valid - observed survival or censoring time in
the validation set; vector of size `n-m` d.valid - event indicator in
the validation set; vector of size `n-m` coded as 0/1 trt.valid -
treatment received in the validation set; vector of size `n-m` coded as
0/1 x.cate.valid - baseline covariates for the outcome model in the
validation set; matrix of dimension `n-m` by `p.cate` x.ps.valid -
baseline covariates (plus intercept) for the propensity score model in
the validation set; matrix of dimension `n-m` by `p.ps + 1`
x.ipcw.valid - baseline covariates for inverse probability of censoring
in the validation set; matrix of dimension `n-m` by `p.ipw` yf.valid -
follow-up time in the training set; if known, vector of size `n-m`; if
unknown, `yf == NULL`
