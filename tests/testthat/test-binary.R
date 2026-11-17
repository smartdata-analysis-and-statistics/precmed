test_that("binary ATE matches the observed risk difference in intercept-only models", {
  d <- data.frame(trt = rep(0:1, each = 50), y = c(rep(1, 10), rep(0, 40), rep(1, 30), rep(0, 20)))
  fit <- atefit("binary", d, y ~ 1, trt ~ 1, n.boot = 20, seed = 12)
  expect_equal(unname(coef(fit)), .4, tolerance = 1e-7)
  expect_equal(fit$risk0$estimate, .2, tolerance = 1e-7)
  expect_equal(fit$risk1$estimate, .6, tolerance = 1e-7)
  expect_equal(dim(confint(fit)), c(1L, 2L))
  expect_s3_class(plot(fit), "ggplot")
  expect_equal(fit$trt.boot, atefit("binary", d, y ~ 1, trt ~ 1, n.boot = 20, seed = 12)$trt.boot)
})

test_that("binary CATE models return probabilities on the risk difference scale", {
  set.seed(42)
  d <- data.frame(trt = rep(0:1, 100), x = rnorm(200))
  d$y <- rbinom(200, 1, plogis(-.5 + .5 * d$trt + .3 * d$x))
  fit <- catefit("binary", d, c("logistic", "twoReg"), y ~ x, trt ~ 1, prop.cutoff = 1)
  expect_equal(length(fit$score.logistic), nrow(d))
  expect_true(all(abs(fit$score.logistic) <= 1))
  beta <- fit$fit$twoReg$coefficients
  expect_equal(fit$score.twoReg, drop(model.matrix(y ~ x, d) %*% (beta[, 2] - beta[, 1])))
  expect_equal(coef(fit)$twoReg, unname(beta[, 2] - beta[, 1]))
  expect_identical(rownames(coef(fit)), rownames(beta))
  expect_match(paste(capture.output(print(fit)), collapse = "\n"), "twoReg")
  expect_false(any(grepl("^NULL$", capture.output(print(fit)))))
  expect_equal(unname(fit$ate.logistic), unname(fit$ate.twoReg))
  expect_output(print(fit), "logistic")
})

test_that("public binary CATE interface matches catefit and validates methods", {
  d <- simulate_binary(n = 200, seed = 42)
  direct <- precmed::catefitbinary(d, "logistic", y ~ age, trt ~ 1,
                                  prop.cutoff = 1, seed = 19)
  dispatched <- catefit("binary", d, "logistic", y ~ age, trt ~ 1,
                        prop.cutoff = 1, seed = 19)
  expect_equal(direct, dispatched)
  expect_error(precmed::catefitbinary(d, "poisson", y ~ age, trt ~ 1),
               "score.method")
  expect_error(precmed::catefitbinary(d[d$trt == levels(d$trt)[1], ],
                                      "logistic", y ~ age, trt ~ 1))
})

test_that("binary logistic CATE uses separate arm-specific outcome models", {
  d <- simulate_binary(n = 500, seed = 442)
  fit <- catefit("binary", d, "logistic", y ~ age, trt ~ age, prop.cutoff = 1)
  control <- glm(y ~ age, data = d, subset = trt == levels(trt)[1], family = binomial())
  treated <- glm(y ~ age, data = d, subset = trt == levels(trt)[2], family = binomial())
  expected <- predict(treated, d, type = "response") - predict(control, d, type = "response")
  expect_equal(unname(fit$score.logistic), unname(expected), tolerance = 1e-8)
  expect_equal(unname(fit$fit$logistic$control$coefficients), unname(coef(control)), tolerance = 1e-8)
  expect_equal(unname(fit$fit$logistic$treated$coefficients), unname(coef(treated)), tolerance = 1e-8)
  bad <- d
  bad$age[bad$trt == levels(bad$trt)[1]] <- 1
  expect_error(catefit("binary", bad, "logistic", y ~ age, trt ~ 1), "rank deficient")
})

test_that("binary fits without linear score coefficients omit NULL when printed", {
  d <- simulate_binary(n = 200, seed = 24)
  fit <- catefit("binary", d, "logistic", y ~ age, trt ~ 1, prop.cutoff = 1)
  expect_null(coef(fit))
  expect_false(any(grepl("^NULL$", capture.output(print(fit)))))
  intercept <- catefit("binary", d, "twoReg", y ~ 1, trt ~ 1,
                       prop.cutoff = 1, B = 1, Kfold = 2, seed = 24)
  expect_identical(rownames(coef(intercept)), "(Intercept)")
  expect_equal(unname(rep(coef(intercept)$twoReg, nrow(d))), unname(intercept$score.twoReg))
})

test_that("binary twoReg cross-fits augmented outcomes and projects both arms", {
  set.seed(62)
  d <- data.frame(trt = rep(0:1, each = 100), x = rnorm(200))
  d$y <- rbinom(200, 1, plogis(-.5 + .4 * d$trt + .3 * d$x))
  prepared <- precmed:::binary_data(d, y ~ x, trt ~ x)
  set.seed(9)
  actual <- precmed:::binary_tworeg_dr(prepared, "glm", .01, .99, list(), 1, 2)
  set.seed(9)
  fold <- integer(nrow(d))
  for (a in 0:1) {
    idx <- which(d$trt == a)
    fold[idx] <- sample(rep(1:2, length.out = length(idx)))
  }
  m0 <- m1 <- p <- numeric(nrow(d))
  for (k in 1:2) {
    train <- fold != k; valid <- !train
    f0 <- glm(y ~ x, data = d, subset = train & trt == 0, family = binomial())
    f1 <- glm(y ~ x, data = d, subset = train & trt == 1, family = binomial())
    fp <- glm(trt ~ x, data = d, subset = train, family = binomial())
    m0[valid] <- predict(f0, d[valid, ], type = "response")
    m1[valid] <- predict(f1, d[valid, ], type = "response")
    p[valid] <- pmin(.99, pmax(.01, predict(fp, d[valid, ], type = "response")))
  }
  z0 <- m0 + (1 - d$trt) * (d$y - m0) / (1 - p)
  z1 <- m1 + d$trt * (d$y - m1) / p
  expected <- cbind(control = coef(lm(z0 ~ x, d)), treated = coef(lm(z1 ~ x, d)))
  expect_equal(actual$fit$coefficients, expected, tolerance = 1e-8)
  expect_error(precmed:::binary_tworeg_dr(prepared, "glm", .01, .99, list(), 1, 101), "at least Kfold")
  expect_error(precmed:::binary_tworeg_dr(prepared, "glm", .01, .99, list(), 0, 2), "B must")
})

test_that("binary twoReg recovers a risk-difference projection with either nuisance model correct", {
  set.seed(124)
  n <- 12000
  x <- rnorm(n)
  trt <- rbinom(n, 1, plogis(-.2 + .6 * x))
  for (outcome_correct in c(FALSE, TRUE)) {
    eta <- -.6 + .2 * x + if (outcome_correct) 0 else .35 * x^2
    mu0 <- plogis(eta); mu1 <- plogis(eta + .5)
    d <- data.frame(x = x, trt = trt, y = rbinom(n, 1, ifelse(trt == 1, mu1, mu0)))
    prepared <- precmed:::binary_data(d, y ~ x, if (outcome_correct) trt ~ 1 else trt ~ x)
    fit <- precmed:::binary_tworeg_dr(prepared, "glm", .001, .999, list(), 1, 2)
    target <- coef(lm(I(mu1 - mu0) ~ x))
    estimate <- fit$fit$coefficients[, 2] - fit$fit$coefficients[, 1]
    expect_equal(unname(estimate), unname(target), tolerance = .04)
  }
})

test_that("binary fitting checks input and handles invalid subgroups", {
  d <- data.frame(y = c(0, 1, 0, 1), trt = c(0, 0, 1, 1), x = 1:4)
  expect_error(atefit("binary", transform(d, y = 2), y ~ 1, trt ~ 1, n.boot = 2), "coded 0/1")
  expect_error(atefit("binary", transform(d, trt = 1), y ~ 1, trt ~ 1, n.boot = 2), "two non-missing")
  expect_warning(fit <- catefit("binary", d, "logistic", y ~ 1, trt ~ 1, prop.cutoff = .25), NA)
  expect_true(is.finite(fit$ate.logistic))
})

test_that("binary effect uses the augmented inverse propensity formula", {
  set.seed(23)
  d <- data.frame(trt = rep(0:1, 80), x = rnorm(160))
  d$y <- rbinom(nrow(d), 1, plogis(-.4 + .5*d$trt + .4*d$x))
  prepared <- precmed:::binary_data(d, y ~ x, trt ~ x)
  pred <- precmed:::binary_predictions(prepared$y, prepared$trt, prepared$x, "logistic")
  p <- fitted(glm(trt ~ x, family = binomial(), data = d))
  expected <- mean(pred$mu1 - pred$mu0 + d$trt * (d$y - pred$mu1) / p - (1-d$trt) * (d$y - pred$mu0) / (1-p))
  actual <- precmed:::binary_effect(prepared, "glm", .01, .99, list())
  expect_equal(unname(actual[3]), expected, tolerance = 1e-8)
  expect_error(precmed:::binary_effect(precmed:::binary_subset(prepared, d$trt == 0), "glm", .01, .99, list()), "Both treatment")
})

test_that("logical treatments match numeric treatments in binary fits", {
  set.seed(41)
  d <- data.frame(trt = rep(0:1, 100), x = rnorm(200))
  d$y <- rbinom(nrow(d), 1, plogis(-0.5 + 0.4 * d$trt + 0.2 * d$x))
  logical_data <- transform(d, trt = as.logical(trt))
  numeric_fit <- atefit("binary", d, y ~ x, trt ~ 1, n.boot = 10, seed = 10)
  logical_fit <- atefit("binary", logical_data, y ~ x, trt ~ 1, n.boot = 10, seed = 10)
  expect_equal(coef(logical_fit), coef(numeric_fit))
  expect_equal(logical_fit$trt.boot, numeric_fit$trt.boot)
  expect_identical(logical_fit$treatment.labels, c("FALSE", "TRUE"))
  numeric_cate <- catefit("binary", d, c("logistic", "twoReg"), y ~ x, trt ~ 1, prop.cutoff = 1, seed = 19)
  logical_cate <- catefit("binary", logical_data, c("logistic", "twoReg"), y ~ x, trt ~ 1, prop.cutoff = 1, seed = 19)
  expect_equal(logical_cate$score.logistic, numeric_cate$score.logistic)
  expect_equal(logical_cate$score.twoReg, numeric_cate$score.twoReg)
  expect_equal(logical_cate$ate.logistic, numeric_cate$ate.logistic)
  expect_error(atefit("binary", transform(d, trt = TRUE), y ~ x, trt ~ 1, n.boot = 2), "two non-missing")
  expect_error(catefit("binary", transform(d, trt = NA), "twoReg", y ~ x, trt ~ 1), "two non-missing")
})

test_that("exported atefitbinary matches binary atefit dispatch", {
  d <- simulate_binary(n = 300, seed = 11)
  direct <- atefitbinary(d, y ~ age, trt ~ age, n.boot = 10, seed = 42)
  dispatched <- atefit("binary", d, y ~ age, trt ~ age, n.boot = 10, seed = 42)
  expect_identical(direct, dispatched)
  expect_identical(getExportedValue("precmed", "atefitbinary"), atefitbinary)
})
