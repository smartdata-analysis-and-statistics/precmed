test_that("Poisson nuisance predictions are rates at unit exposure", {
  n <- 80L
  trt <- rep(0:1, each = n / 2)
  x <- cbind(x1 = rep(c(0, 1), length.out = n))
  time <- log(rep(c(0.5, 2), each = 2, length.out = n))
  set.seed(41)
  y <- rpois(n, exp(0.3 + 0.2 * x[, 1] + 0.1 * trt + time))
  predictions <- list()
  local_mocked_bindings(
    onearmglmcount.dr = function(y, x.cate, time, trt, ps, f.predictor) {
      predictions[[length(predictions) + 1L]] <<- f.predictor
      rep(0, ncol(x.cate) + 1L)
    },
    .package = "precmed"
  )

  set.seed(12)
  index1 <- sample(precmed:::generate_kfold_indices(n / 2, 2, reverse = FALSE))
  index0 <- sample(precmed:::generate_kfold_indices(n / 2, 2, reverse = TRUE))
  fold <- integer(n)
  fold[trt == 1] <- index1
  fold[trt == 0] <- index0
  set.seed(12)
  precmed:::intxcount(y, trt, x, cbind(1, x), time,
    score.method = "twoReg", initial.predictor.method = "poisson",
    B = 1, Kfold = 2)

  expect_length(predictions, 2)
  for (r in 0:1) {
    for (k in 1:2) {
      train <- fold != k & trt == r
      valid <- fold == k
      dat <- data.frame(y = y, x1 = x[, 1], time = time)
      fit <- glm(y ~ x1 + offset(time), data = dat[train, ], family = poisson())
      expected <- predict(fit, newdata = transform(dat[valid, ], time = 0), type = "response")
      expect_equal(unname(predictions[[2L - r]][valid]), unname(expected))
    }
  }
})
