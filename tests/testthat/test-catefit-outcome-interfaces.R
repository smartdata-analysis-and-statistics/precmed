test_that("outcome-specific CATE interfaces agree with the generic and legacy calls", {
  count <- simulate_count(n = 400, seed = 421)
  binary <- simulate_binary(n = 400, seed = 422)
  surv <- simulate_survival(n = 400, seed = 423)
  cases <- list(
    list(response = "count", fun = catefitcount, data = count,
         formula = y ~ age + offset(log(years)), outcome = "poisson", old = "poisson", extra = list()),
    list(response = "binary", fun = catefitbinary, data = binary,
         formula = y ~ age, outcome = "logistic", old = "logistic", extra = list()),
    list(response = "survival", fun = catefitsurv, data = surv,
         formula = survival::Surv(y, d) ~ age, outcome = "poisson", old = "poisson",
         extra = list(tau0 = median(surv$y))))
  for (case in cases) {
    expect_false(any(c("score.method", "cate.model") %in% names(formals(case$fun))))
    expect_identical(names(formals(case$fun))[1:7],
                     c("data", "outcome_model", "outcome_method", "benefit_model", "benefit_method", "ps_model", "ps_method"))
    common <- c(list(data = case$data, ps_model = trt ~ age,
                     prop.cutoff = c(.8, 1), seed = 424), case$extra)
    modern <- c(common, list(outcome_model = case$formula,
                             benefit_method = "t_learner", outcome_method = case$outcome))
    legacy <- c(common, list(cate.model = case$formula, score.method = case$old))
    direct <- suppressWarnings(suppressMessages(do.call(case$fun, modern)))
    generic <- suppressWarnings(suppressMessages(do.call(catefit, c(list(response = case$response), modern))))
    previous <- suppressWarnings(suppressMessages(do.call(case$fun, legacy)))
    positional <- suppressWarnings(suppressMessages(do.call(case$fun, c(
      list(case$data, case$formula, case$outcome, "t_learner", trt ~ age),
      list(prop.cutoff = c(.8, 1), seed = 424), case$extra))))
    expect_equal(positional, direct)
    expect_equal(direct, generic)
    aliases <- modern
    names(aliases)[names(aliases) == "ps_model"] <- "ps.model"
    aliases$ps.method <- "glm"
    expect_warning(aliased <- suppressMessages(do.call(case$fun, aliases)), "deprecated")
    expect_equal(aliased, direct)
    expect_error(do.call(case$fun, c(modern, list(ps.model = trt ~ 1))), "Specify only ps_model")
    expect_equal(direct, previous)
    expect_error(do.call(case$fun, c(modern, list(cate.model = case$formula))), "only outcome_model")
    expect_error(do.call(case$fun, c(modern, list(score.method = case$old))), "only benefit_method")
    expect_error(do.call(case$fun, c(common, list(outcome_model = case$formula,
      benefit_method = "two_regression"))), "outcome_method must be supplied")
  }
})

test_that("outcome-specific methods resolve initial outcome models", {
  d <- simulate_count(n = 400, seed = 425)
  args <- list(data = d, outcome_model = y ~ age + offset(log(years)),
               ps_model = trt ~ age, benefit_method = "two_regression",
               outcome_method = "poisson", B = 1, Kfold = 2, prop.cutoff = 1, seed = 426)
  direct <- suppressWarnings(do.call(catefitcount, args))
  generic <- suppressWarnings(do.call(catefit, c(list(response = "count"), args)))
  expect_equal(direct, generic)
  expect_error(suppressWarnings(do.call(catefitcount, c(args, list(initial.predictor.method = "boosting")))), "only outcome_method")
})

test_that("continuous CATE interface uses the same argument order", {
  expect_false(any(c("score.method", "cate.model") %in% names(formals(catefitmean))))
  expect_identical(names(formals(catefitmean))[1:7],
                   c("data", "outcome_model", "outcome_method", "benefit_model", "benefit_method", "ps_model", "ps_method"))
})

test_that("CATE treatment model aliases agree and reject conflicts", {
  d <- simulate_binary(n = 250, seed = 517)
  args <- list(data = d, outcome_model = y ~ age, outcome_method = "logistic",
               benefit_method = "t_learner", prop.cutoff = 1, seed = 518)
  for (fun in list(catefitbinary, function(...) catefit("binary", ...))) {
    modern <- do.call(fun, c(args, list(ps_model = trt ~ age, ps_method = "glm")))
    expect_warning(old <- do.call(fun, c(args, list(ps.model = trt ~ age, ps.method = "glm"))), "deprecated")
    expect_equal(old, modern)
    expect_error(do.call(fun, c(args, list(ps_model = trt ~ age, ps.model = trt ~ 1))), "Specify only ps_model")
    expect_error(do.call(fun, c(args, list(ps_model = trt ~ age, ps_method = "glm", ps.method = "lasso"))), "Specify only ps_method")
  }
})
