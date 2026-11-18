test_that("catefit scores require an available method and reject new data", {
  fit <- structure(list(response = "survival", score.contrastReg = c(.1, .2),
                        score.randomForest = c(.3, .4),
                        coefficients = data.frame(contrastReg = c(.1, .2),
                                                  row.names = c("(Intercept)", "age"))), class = "catefit")
  expect_equal(predict(fit, method = "contrastReg"), fit$score.contrastReg)
  expect_equal(predict(fit, method = "contrast_regression"), fit$score.contrastReg)
  expect_equal(predict(fit, method = "contrastreg"), fit$score.contrastReg)
  expect_equal(predict(fit, method = "random_forest"), fit$score.randomForest)
  expect_equal(predict(fit, method = "t_learner"), fit$score.randomForest)
  expect_error(predict(fit), "Specify one")
  expect_error(predict(fit, method = "boosting"), "available")
  expect_error(predict(fit, newdata = data.frame(age = 1), method = "contrastReg"), "new data")
  expect_equal(coef(fit), fit$coefficients)
  expect_equal(coef(fit, method = "contrastReg"), c(`(Intercept)` = .1, age = .2))
  expect_equal(coef(fit, method = "contrast_regression"), coef(fit, method = "contrastReg"))
  expect_equal(coef(fit, method = "contrastreg"), coef(fit, method = "contrastReg"))
  expect_error(coef(fit, method = "randomForest"), "No score coefficients")
  expect_error(confint(fit), "standard errors are not available")
  fit$score.randomForest <- NULL
  expect_equal(predict(fit), c(.1, .2))
})

test_that("CATE prediction identifies each score in a multi-method fit", {
  fit <- structure(list(response = "binary", score.logistic = c(.1, .2),
                        score.twoReg = c(.3, .4), score.boosting = c(.5, .6)),
                   class = "catefit")
  expect_equal(predict(fit, method = "two_regression"), fit$score.twoReg)
  expect_equal(predict(fit, method = "twoReg"), fit$score.twoReg)
  expect_equal(predict(fit, method = "logistic"), fit$score.logistic)
  expect_equal(predict(fit, method = "boosting"), fit$score.boosting)
  expect_error(predict(fit, method = "t_learner"), "logistic, boosting")
  fit$score.logistic <- NULL
  expect_equal(predict(fit, method = "t_learner"), fit$score.boosting)
  fit$score.boosting <- NULL
  expect_error(predict(fit, method = "t_learner"), "none")
  count <- structure(list(score.negBin = 1:2), class = "catefit")
  expect_equal(predict(count, method = "negative_binomial"), 1:2)
})

test_that("contrast regression intervals use available coefficient standard errors", {
  fit <- structure(list(response = "count", coefficients = data.frame(
    contrastReg = c(.1, .2), SE_contrastReg = c(.02, .04),
    row.names = c("(Intercept)", "age"))), class = "catefit")
  interval <- confint(fit, method = "contrastReg", level = .90)
  expect_equal(confint(fit, method = "contrast_regression", level = .90), interval)
  expect_equal(confint(fit, method = "contrastreg", level = .90), interval)
  expect_equal(unname(interval[2, ]), .2 + c(-1, 1) * qnorm(.95) * .04)
  expect_equal(rownames(confint(fit, parm = "age")), "age")
  expect_error(confint(fit, parm = Inf), "parm")
  expect_error(confint(fit, level = 1), "level")
  fit$response <- "survival"
  expect_error(confint(fit), "standard errors are not available")
})

test_that("coefficient selectors respect T-learner ambiguity and availability", {
  fit <- structure(list(response = "count", score.poisson = 1:2, score.boosting = 1:2,
                        score.twoReg = 1:2, coefficients = data.frame(poisson = c(1, 2),
                        twoReg = c(3, 4), row.names = c("(Intercept)", "age"))), class = "catefit")
  expect_equal(coef(fit, method = "two_regression"), coef(fit, method = "twoReg"))
  expect_error(coef(fit, method = "t_learner"), "poisson, boosting")
  fit$score.boosting <- NULL
  expect_equal(coef(fit, method = "t_learner"), coef(fit, method = "poisson"))
  fit$score.poisson <- NULL
  fit$score.boosting <- 1:2
  expect_error(coef(fit, method = "t_learner"), "No score coefficients")
})
