library(testthat)

# Basic tests for AUC function

# Test for correct linear interpolation AUC
test_that("AUC with linear interpolation is calculated correctly", {
  x <- c(0, 1, 2, 3, 4)
  y <- c(0, 1, 1, 2, 0)
  expect_equal(auc(x, y, type = "linear"), 4)  # Value obtained via MESS
})

# Test for correct spline interpolation AUC
test_that("AUC with spline interpolation is calculated correctly", {
  x <- c(0, 1, 2, 3, 4)
  y <- c(0, 1, 1, 2, 0)
  expect_equal(auc(x, y, type = "spline", subdivisions = 100),
               integrate(splinefun(x, y, method = "natural"), 0, 4)$value)
})

# Test when x and y lengths differ
test_that("AUC throws error when x and y lengths differ", {
  x <- c(0, 1, 2)
  y <- c(0, 1)
  expect_error(auc(x, y, type = "linear"), "length\\(x\\) == length\\(y\\)")
})

# Test for handling of NA values
test_that("AUC handles NA values in from or to", {
  x <- c(0, 1, 2, 3, 4)
  y <- c(0, 1, 1, 2, 0)
  expect_error(auc(x, y, from = NA), "is.na\\(from\\) is not TRUE")
})

# Test for minimal input, when x has fewer than two unique values
test_that("AUC returns NA for less than two unique x values", {
  x <- c(1, 1, 1)
  y <- c(0, 1, 2)
  expect_true(is.na(auc(x, y, type = "linear")))
})

# Test AUC for a constant function (flat line)
test_that("AUC is calculated correctly for a constant function", {
  x <- c(0, 1, 2, 3, 4)
  y <- c(1, 1, 1, 1, 1)
  expect_equal(auc(x, y, type = "linear"), 4)  # Area under constant y=1 from 0 to 4 is 4
})

# Test AUC for decreasing y values
test_that("AUC is calculated correctly for decreasing y values", {
  x <- c(0, 1, 2, 3, 4)
  y <- c(4, 3, 2, 1, 0)
  expect_equal(auc(x, y, type = "linear"), 8)  # Value obtained via MESS
})

# Test AUC for non-default from and to
test_that("AUC is calculated correctly for custom from and to", {
  x <- c(0, 1, 2, 3, 4)
  y <- c(0, 1, 1, 2, 0)
  expect_equal(auc(x, y, from = 1, to = 3, type = "linear"), 2.5) # Value obtained via MESS
})
