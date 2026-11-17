test_that("catefit scores require an available method and reject new data", {
  fit <- structure(list(response = "survival", score.contrastReg = c(.1, .2),
                        score.randomForest = c(.3, .4),
                        coefficients = data.frame(contrastReg = c(.1, .2),
                                                  row.names = c("(Intercept)", "age"))), class = "catefit")
  expect_equal(predict(fit, method = "contrastReg"), fit$score.contrastReg)
  expect_error(predict(fit), "Specify one")
  expect_error(predict(fit, method = "boosting"), "available")
  expect_error(predict(fit, newdata = data.frame(age = 1), method = "contrastReg"), "new data")
  expect_equal(coef(fit), fit$coefficients)
  expect_equal(coef(fit, method = "contrastReg"), c(`(Intercept)` = .1, age = .2))
  expect_error(coef(fit, method = "randomForest"), "No score coefficients")
  expect_error(confint(fit), "standard errors are not available")
  fit$score.randomForest <- NULL
  expect_equal(predict(fit), c(.1, .2))
})

test_that("contrast regression intervals use available coefficient standard errors", {
  fit <- structure(list(response = "count", coefficients = data.frame(
    contrastReg = c(.1, .2), SE_contrastReg = c(.02, .04),
    row.names = c("(Intercept)", "age"))), class = "catefit")
  interval <- confint(fit, method = "contrastReg", level = .90)
  expect_equal(unname(interval[2, ]), .2 + c(-1, 1) * qnorm(.95) * .04)
  expect_equal(rownames(confint(fit, parm = "age")), "age")
  expect_error(confint(fit, parm = Inf), "parm")
  expect_error(confint(fit, level = 1), "level")
  fit$response <- "survival"
  expect_error(confint(fit), "standard errors are not available")
})
