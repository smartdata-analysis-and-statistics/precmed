test_that("survival simulator is reproducible and assigns treatment by propensity", {
  dat <- simulate_survival(n = 600, seed = 18)
  expect_identical(dat, simulate_survival(n = 600, seed = 18))
  expect_false(identical(dat, simulate_survival(n = 600, seed = 19)))
  expect_equal(nrow(dat), 600L)
  expect_gt(sum(dat$trt == "drug0"), 0)
  expect_gt(sum(dat$trt == "drug1"), 0)
  expect_true(all(dat$y > 0 & dat$d %in% c(0, 1)))
  expect_equal(attr(dat$age, "scaled:center"), 48)
  expect_equal(attr(dat$previous_cost, "scaled:scale"), 24266)
})

test_that("survival simulator responds to propensity and outcome controls", {
  base <- simulate_survival(seed = 19)
  ps <- stats::predict(stats::glm(trt ~ age + previous_treatment, data = base,
                                  family = stats::binomial()), type = "response")
  expect_equal(mean(base$trt == "drug1"), 2151 / 4000, tolerance = 0.03)
  expect_gt(min(ps), 0.05)
  expect_lt(max(ps), 0.95)
  tilted <- simulate_survival(seed = 19,
    ps_weights = c(age = 1))
  age_gap <- function(x) diff(tapply(x$age, x$trt, mean))
  expect_gt(unname(age_gap(tilted)), unname(age_gap(base)))
  low <- simulate_survival(seed = 19, conditional_hr = 0.5)
  high <- simulate_survival(seed = 19, conditional_hr = 2)
  expect_gt(sum(high$d[high$trt == "drug1"]), sum(low$d[low$trt == "drug1"]))
  expect_identical(high$y[high$trt == "drug0"], low$y[low$trt == "drug0"])
  changed <- simulate_survival(seed = 19,
    lp_weights = c(treatment_age = 0.8))
  expect_false(identical(base$y, changed$y))
  no_interaction <- simulate_survival(seed = 19,
    lp_weights = c(treatment_age = 0, treatment_relapses = 0))
  older_treated_default <- base$trt == "drug1" & base$age > 0
  expect_true(all(base$y[older_treated_default] >= no_interaction$y[older_treated_default]))
  expect_true(any(base$y[older_treated_default] > no_interaction$y[older_treated_default]))
  relapse_interaction <- simulate_survival(seed = 19,
    lp_weights = c(treatment_age = 0, treatment_relapses = -0.4))
  zero_relapses <- no_interaction$previous_number_relapses == 0
  expect_identical(no_interaction$y[zero_relapses], relapse_interaction$y[zero_relapses])
  treated_with_relapses <- no_interaction$trt == "drug1" & !zero_relapses
  expect_lt(sum(relapse_interaction$d[treated_with_relapses]),
            sum(no_interaction$d[treated_with_relapses]))
  older_treated <- no_interaction$trt == "drug1" & no_interaction$age > 0
  age_interaction <- simulate_survival(seed = 19,
    lp_weights = c(treatment_age = -0.4, treatment_relapses = 0))
  expect_lt(sum(age_interaction$d[older_treated]), sum(no_interaction$d[older_treated]))
  expect_identical(age_interaction$y[age_interaction$trt == "drug0"],
                   no_interaction$y[no_interaction$trt == "drug0"])
})

test_that("survival simulator rejects invalid inputs and allows extreme allocation", {
  expect_error(simulate_survival(n = 0), "positive integer")
  expect_error(simulate_survival(n = 1.5), "positive integer")
  expect_error(simulate_survival(ps_weights = c(unknown = 1)), "ps_weights")
  expect_error(simulate_survival(lp_weights = c(unknown = 1)), "lp_weights")
  expect_error(simulate_survival(conditional_hr = 0), "conditional_hr")
  expect_error(simulate_survival(seed = 1.5), "seed")
  extreme <- simulate_survival(n = 100, seed = 1,
    ps_weights = c(intercept = 1000, age = 0, drugB = 0, drugC = 0))
  expect_true(all(extreme$trt == "drug1"))
  lower <- simulate_survival(n = 1000, seed = 1,
    ps_weights = c(intercept = -2, age = 0, drugB = 0, drugC = 0))
  higher <- simulate_survival(n = 1000, seed = 1,
    ps_weights = c(intercept = 1, age = 0, drugB = 0, drugC = 0))
  expect_gt(sum(higher$trt == "drug1"), sum(lower$trt == "drug1"))
})
