test_that("CATE default T-learners match explicit configurations", {
  cases <- list(
    list(response = "binary", fun = catefitbinary, data = simulate_binary(n = 400, seed = 951),
         model = y ~ age, method = "logistic", extra = list()),
    list(response = "count", fun = catefitcount, data = simulate_count(n = 400, seed = 952),
         model = y ~ age + offset(log(years)), method = "poisson", extra = list()),
    list(response = "survival", fun = catefitsurv, data = simulate_survival(n = 400, seed = 953),
         model = survival::Surv(y, d) ~ age, method = "coxph", extra = list()))
  for (case in cases) {
    extra <- case$extra
    if (case$response == "survival") extra$tau0 <- median(case$data$y) * 1.1
    args <- c(list(data = case$data, outcome_model = case$model, ps_model = trt ~ age,
                   prop.cutoff = c(.8, 1), seed = 954), extra)
    explicit <- suppressWarnings(suppressMessages(do.call(case$fun, c(args, list(outcome_method = case$method, benefit_method = "t_learner")))))
    direct <- suppressWarnings(suppressMessages(do.call(case$fun, args)))
    generic <- suppressWarnings(suppressMessages(do.call(catefit, c(list(response = case$response), args))))
    partial <- suppressWarnings(suppressMessages(do.call(case$fun, c(args, list(benefit_method = "t_learner")))))
    selected <- suppressWarnings(suppressMessages(do.call(case$fun, c(args, list(outcome_method = case$method)))))
    nulls <- suppressWarnings(suppressMessages(do.call(case$fun, c(args, list(outcome_method = NULL, benefit_method = NULL)))))
    expect_equal(direct, explicit, info = case$response)
    expect_equal(generic, explicit, info = case$response)
    expect_equal(partial, explicit, info = case$response)
    expect_equal(selected, explicit, info = case$response)
    expect_equal(nulls, explicit, info = case$response)
    expect_error(do.call(case$fun, c(args, list(benefit_method = "two_regression"))), "outcome_method must be supplied")
  }
})
