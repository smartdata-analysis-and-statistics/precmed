test_that("ABC summaries retain method names and preserve missing values", {
  fit <- structure(list(abc = TRUE, cv.n = 3,
    ate.first = list(abc.valid = c(1, 2, 3)),
    ate.second = list(abc.valid = c(2, NA, 6))), class = "precmed")
  values <- abc(fit)
  expect_equal(dim(values), c(2L, 3L))
  expect_equal(abc(fit, combine = "mean"), c(first = 2, second = NA_real_))
  expect_equal(abc(fit, combine = mean), apply(values, 1, mean))
  expect_equal(abc(fit, combine = "median"), abc(fit, combine = median))
  expect_equal(abc(fit, combine = function(z) mean(z, na.rm = TRUE)),
               c(first = 2, second = 4))
  expect_error(abc(fit, combine = "sum"), "combine")
  expect_error(abc(fit, combine = function(z) z), "one numeric")
  expect_error(abc(fit, combine = function(z) "text"), "one numeric")
  fit$cv.n <- 1
  fit$ate.first$abc.valid <- 3
  fit$ate.second <- NULL
  expect_equal(abc(fit, combine = "mean"), c(first = 3))
})
