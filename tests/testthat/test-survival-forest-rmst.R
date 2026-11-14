test_that("forest survival integration includes the interval after the last event", {
  integrate <- precmed:::integrate_forest_survival
  curves <- rbind(c(0.8, 0.5), c(0.9, 0.7))
  expect_equal(integrate(c(1, 3), curves, 5), c(3.6, 4.2))
  expect_equal(integrate(c(1, 3), curves, 2), c(1.8, 1.9))
  expect_equal(integrate(c(1, 3), curves, 3), c(2.6, 2.8))
  expect_equal(integrate(1, matrix(0.8, nrow = 1), 3), 2.6)
})

test_that("horizon event status distinguishes censoring from an event", {
  status <- precmed:::survival_event_at_horizon
  expect_equal(status(c(1, 2, 2, 3, 3), c(1, 1, 0, 1, 0), 2),
               c(1, 1, 0, 0, 0))
})

test_that("out-of-range RMST fallback warns and valid estimates are preserved", {
  fallback <- precmed:::rmst_outcome_fallback
  expect_equal(fallback(2, c(1, 3), 5, 1), 2)
  expect_warning(value <- fallback(-1, c(1, 3), 5, 1), "not doubly robust")
  expect_equal(value, 2)
  expect_warning(value <- fallback(6, c(2, 4), 5, 0), "treatment 0")
  expect_equal(value, 3)
})
