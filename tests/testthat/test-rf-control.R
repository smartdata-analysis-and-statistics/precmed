test_that("forest controls reject invalid and reserved options", {
  expect_equal(precmed:::resolve_rf_control(list(mtry = 2, nodesize = 5, nsplit = 0)),
               list(mtry = 2, nodesize = 5, nsplit = 0))
  expect_error(precmed:::resolve_rf_control(list(ntree = 10)), "n.trees.rf")
  expect_error(precmed:::resolve_rf_control(list(forest = FALSE)), "rf.control")
  expect_error(precmed:::resolve_rf_control(list(mtry = 1.5)), "mtry")
  expect_error(precmed:::resolve_rf_control(list(nodesize = 0)), "nodesize")
  expect_error(precmed:::resolve_rf_control(list(splitrule = "")), "splitrule")
  expect_error(precmed:::resolve_rf_control(list(mtry = 1, mtry = 2)), "rf.control")
})

test_that("forest settings reach randomForestSRC and predictions work", {
  skip_if_not_installed("randomForestSRC")
  set.seed(53)
  dat <- data.frame(y = rexp(100), d = rbinom(100, 1, .8),
                    x = rnorm(100), z = rnorm(100))
  fit <- precmed:::fit_rf(survival::Surv(y, d) ~ ., dat, 20,
                          list(mtry = 1, nodesize = 7))
  expect_equal(fit$ntree, 20)
  expect_equal(fit$mtry, 1)
  expect_equal(fit$nodesize, 7)
  expect_equal(nrow(predict(fit, newdata = dat[1:4, ])$survival), 4)
})
