test_that("count simulator creates reproducible analysis-ready data", {
  d <- simulate_count(n = 500, seed = 12)
  expect_identical(d, simulate_count(n = 500, seed = 12))
  expect_false(identical(d, simulate_count(n = 500, seed = 13)))
  expect_named(d, c("age", "female", "previous_treatment", "previous_cost",
                    "previous_number_symptoms", "previous_number_relapses", "trt", "y", "years"))
  expect_true(all(d$y >= 0 & d$y == floor(d$y)))
  expect_true(all(d$years >= 0.5 & d$years <= 3))
  expect_equal(attr(d$age, "scaled:center"), 48)
  expect_equal(attr(d$previous_cost, "scaled:scale"), 20172)
  expect_true(all(simulate_count(n = 100, followup_range = c(2, 2))$years == 2))
  expect_true(all(is.finite(simulate_count(n = 100, dispersion = 2)$y)))
  fit <- atefitcount(data = d, outcome_model = y ~ age + offset(log(years)),
                     ps_model = trt ~ age, n.boot = 5, seed = 1)
  expect_true(is.finite(fit$log.rate.ratio$estimate))
})

test_that("count simulator rejects invalid settings", {
  expect_error(simulate_count(n = 0), "positive integer")
  expect_error(simulate_count(seed = 1.5), "seed")
  expect_error(simulate_count(ps_weights = c(other = 1)), "ps_weights")
  expect_error(simulate_count(lp_weights = c(other = 1)), "lp_weights")
  expect_error(simulate_count(conditional_rr = 0), "conditional_rr")
  expect_error(simulate_count(followup_range = c(2, 1)), "followup_range")
  expect_error(simulate_count(dispersion = 0), "dispersion")
})
