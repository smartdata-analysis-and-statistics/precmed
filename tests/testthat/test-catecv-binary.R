test_that("binary CV scores held-out data with the training model", {
  set.seed(129)
  d <- data.frame(x = rnorm(500), trt = rbinom(500, 1, .5))
  d$y <- rbinom(500, 1, plogis(-.4 + .4*d$trt + .3*d$x + .5*d$trt*d$x))
  args <- list(response = "binary", data = d, score.method = c("logistic", "twoReg"),
               cate.model = y ~ x, ps.model = trt ~ x, cv.n = 2, B = 1, Kfold = 2,
               prop.cutoff = c(.5, .75, 1), abc = FALSE, seed = 27)
  fit <- do.call(catecv, args)
  expect_identical(fit, do.call(catecv, args))
  expect_s3_class(fit, "precmed")
  expect_equal(dim(fit$ate.logistic$ate.est.valid.high.cv), c(3L, 2L))
  for (i in 1:2) {
    train <- fit$splits[[i]]$train; valid <- fit$splits[[i]]$valid
    expect_length(intersect(train, valid), 0)
    outcome <- glm(y ~ trt*x, data = d[train, ], family = binomial())
    validdata <- d[valid, ]
    score <- predict(outcome, transform(validdata, trt = 1), type = "response") -
      predict(outcome, transform(validdata, trt = 0), type = "response")
    selected <- score >= quantile(score, .5)
    sub <- precmed:::binary_data(validdata[selected, ], y ~ x, trt ~ x)
    expected <- unname(precmed:::binary_effect(sub, "glm", .01, .99, list())[3])
    expect_equal(unname(fit$ate.logistic$ate.est.valid.high.cv[1, i]), expected, tolerance = 1e-8)
    expect_equal(unname(fit$ate.twoReg$ate.est.valid.high.cv[3, i]), fit$overall.ate.valid[i])
  }
  expect_s3_class(plot(fit, show.abc = FALSE), "ggplot")
  expect_s3_class(boxplot(fit), "ggplot")
  expect_match(plot(fit, show.abc = FALSE)$labels$y, "Risk difference")
  expect_error(abc(fit), "must have been calculated")
  args$higher.y <- FALSE
  reversed <- do.call(catecv, args)
  train <- reversed$splits[[1]]$train; validdata <- d[reversed$splits[[1]]$valid, ]
  outcome <- glm(y ~ trt*x, data = d[train, ], family = binomial())
  score <- predict(outcome, transform(validdata, trt = 1), type = "response") -
    predict(outcome, transform(validdata, trt = 0), type = "response")
  sub <- precmed:::binary_data(validdata[score <= quantile(score, .5), ], y ~ x, trt ~ x)
  expect_equal(unname(reversed$ate.logistic$ate.est.valid.high.cv[1, 1]),
               unname(precmed:::binary_effect(sub, "glm", .01, .99, list())[3]))
})

test_that("binary CV computes ABC and records unsupported subgroups", {
  set.seed(7)
  d <- data.frame(y = rbinom(240, 1, .4), trt = rep(0:1, 120))
  fit <- suppressWarnings(catecv("binary", d, "logistic", y ~ 1, trt ~ 1,
      cv.n = 1, prop.cutoff = c(.5, 1), seed = 7))
  expect_equal(dim(abc(fit)), c(1L, 1L))
  expect_equal(unname(abc(fit)[1, 1]), 0, tolerance = 1e-8)
  expect_true(any(is.na(fit$ate.logistic$ate.est.valid.group.cv)))
  expect_true(nrow(diagnostics(fit)) > 0)
  expect_true(all(fit$subgroup.support$n0 >= 0))
  expect_s3_class(plot(fit), "ggplot")
})

test_that("binary CV validates settings and reports score failures", {
  d <- data.frame(y = rep(0:1, 20), trt = rep(0:1, each = 20))
  expect_error(catecv("binary", d, "logistic", y ~ 1, trt ~ 1, cv.n = 0), "cv.n")
  expect_error(catecv("binary", d, "logistic", y ~ 1, trt ~ 1, error.max = .1), "Unsupported binary")
  expect_error(catecv("binary", d, "logistic", y ~ 1, trt ~ 1, prop.cutoff = 1), "ABC requires")
  expect_error(catecv("binary", transform(d, trt = 1), "logistic", y ~ 1, trt ~ 1), "two non-missing")
  fit <- suppressWarnings(catecv("binary", d, "twoReg", y ~ 1, trt ~ 1,
    cv.n = 1, Kfold = 100, abc = FALSE, seed = 1))
  expect_true(all(is.na(fit$ate.twoReg$ate.est.valid.high.cv)))
  expect_true(any(diagnostics(fit)$group == "score"))
})

test_that("binary CV learns formula transformations from training data", {
  d <- data.frame(y = rep(0:1, 50), trt = rep(0:1, each = 50), x = 1:100)
  train <- seq_len(100) %% 2 == 0
  design <- precmed:::binary_cv_design(d, train, y ~ scale(x))
  expected <- (d$x[!train] - mean(d$x[train])) / sd(d$x[train])
  expect_equal(unname(design$valid[, 2]), expected)
})

test_that("binary boosting scores validation observations", {
  d <- simulate_binary(n = 400, seed = 10)
  fit <- suppressWarnings(catecv("binary", d, "boosting", y ~ age, trt ~ 1,
    cv.n = 1, abc = FALSE, prop.cutoff = c(.5, 1), seed = 3,
    n.trees.boosting = 15, plot.gbmperf = FALSE))
  expect_true(is.finite(fit$ate.boosting$ate.est.valid.high.cv[2, 1]))
})
