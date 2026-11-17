test_that("ATE outcome_model replaces the deprecated alias without changing results", {
  cases <- list(
    list(response = "binary", data = simulate_binary(n = 400, seed = 10), formula = y ~ age),
    list(response = "count", data = simulate_count(n = 400, seed = 11), formula = y ~ age + offset(log(years))),
    list(response = "survival", data = simulate_survival(n = 400, seed = 12), formula = survival::Surv(y, d) ~ age)
  )
  for (case in cases) {
    args <- list(response = case$response, data = case$data, ps.model = trt ~ age,
                 n.boot = 3, seed = 13)
    if (case$response == "survival") args$tau0 <- median(case$data$y) * 1.1
    new <- suppressWarnings(do.call(atefit, c(args, list(outcome_model = case$formula))))
    expect_warning(old <- do.call(atefit, c(args, list(cate.model = case$formula))),
                   "cate.model is deprecated", info = case$response)
    expect_equal(new, old, info = case$response)
    expect_error(do.call(atefit, c(args, list(outcome_model = case$formula, cate.model = case$formula))),
                 "Specify only outcome_model", info = case$response)
    expect_error(do.call(atefit, args), "outcome_model must be supplied", info = case$response)
  }
})

test_that("direct ATE interfaces accept the alias and positional formulas", {
  d <- simulate_binary(n = 400, seed = 20)
  args <- list(data = d, ps.model = trt ~ age, n.boot = 3, seed = 21)
  new <- do.call(atefitbinary, c(args, list(outcome_model = y ~ age)))
  expect_warning(old <- do.call(atefitbinary, c(args, list(cate.model = y ~ age))), "deprecated")
  expect_equal(new, old)
  expect_equal(new, atefitbinary(d, y ~ age, trt ~ age, n.boot = 3, seed = 21))
  expect_error(atefitbinary(d, y ~ age, trt ~ age, cate.model = y ~ age), "Specify only")
})
