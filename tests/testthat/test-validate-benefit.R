# Author: Thomas Debray.
run_benefit <- function(...) suppressWarnings(validate_benefit(...))

test_that("group calibration and policy values match hand calculations", {
  trt <- rep(c(0, 0, 1, 1), 2)
  benefit <- rep(c(-1, 1), each = 4)
  y <- c(2, 4, 1, 3, 1, 3, 4, 6)
  x <- run_benefit(y, trt, benefit, groups = 2)
  expect_equal(x$calibration$observed, c(-1, 3))
  expect_equal(x$calibration$predicted, c(-1, 1))
  expect_equal(x$mean.bias, 1)
  expect_equal(unname(x$calibration.fit[c("intercept", "slope", "rmse")]), c(1, 2, sqrt(2)))
  expect_equal(unname(x$policy), c(.5, 2, 1.5, .5))
  low <- run_benefit(y, trt, benefit, groups = 2, higher.y = FALSE)
  expect_equal(unname(low$policy["versus.opposite"]), 2)
})

test_that("ties and missing treatment arms are handled explicitly", {
  x <- run_benefit(c(0, 1, 1, 0), c(0, 1, 0, 1), rep(0, 4), groups = 4)
  expect_equal(nrow(x$calibration), 1)
  expect_true(is.na(x$calibration.fit["slope"]))
  expect_equal(unname(x$policy["versus.control"]), 0)
  x <- run_benefit(c(0, 1, 0, 1), c(0, 0, 1, 1), c(-1, -1, 1, 1), groups = 2)
  expect_true(all(is.na(x$calibration$observed)))
  expect_true(is.na(x$policy["versus.opposite"]))
  expect_error(run_benefit(1:3, c(0, 0, 0), 1:3), "both treatment arms")
  expect_error(run_benefit(1:3, c(0, 1, 0), c(1, NA, 2)), "finite")
})

test_that("policy standardization handles unequal randomization", {
  x <- run_benefit(c(0, 0, 2, 0, 0, 4), c(0, 0, 1, 0, 0, 1),
                        c(-1, -1, -1, 1, 1, 1), groups = 2)
  expect_equal(unname(x$policy["versus.opposite"]), 1)
  expect_equal(unname(x$policy["versus.control"]), 2)
})


test_that("propensity weighting matches normalized arm means and policy values", {
  y <- c(0, 1, 1, 0, 1, 0, 1, 1)
  trt <- rep(c(0, 0, 1, 1), 2)
  benefit <- rep(c(-1, 1), each = 4)
  ps <- rep(c(.2, .8, .2, .8), 2)
  x <- run_benefit(y, trt, benefit, groups = 2, propensity = ps)
  expect_equal(x$calibration$observed, c(0, .8))
  expect_equal(x$mean.bias, .4)
  expect_equal(unname(x$policy[c("versus.opposite", "versus.control", "versus.treatment")]),
               c(.4, .4, 0))
  expect_true(all(is.na(x$calibration$se)))
  expect_error(run_benefit(y, trt, benefit, propensity = rep(0, 8)), "strictly")
  equal <- run_benefit(y, trt, benefit, groups = 2, propensity = rep(.3, 8))
  plain <- run_benefit(y, trt, benefit, groups = 2)
  expect_equal(equal$calibration$observed, plain$calibration$observed)
  expect_equal(equal$policy, plain$policy)
})


test_that("benefit validation helper is internal", {
  expect_false("validate_benefit" %in% getNamespaceExports("precmed"))
})


test_that("benefit validation warns that it is experimental", {
  expect_warning(validate_benefit(c(0, 1), c(0, 1), c(0, 0), groups = 1),
                 "experimental and under development")
})
