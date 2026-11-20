test_that("benefit model changes doubly robust scores but not T-learners", {
  cases <- list(
    list(response = "binary", data = simulate_binary(n = 500, seed = 811),
         outcome = y ~ age + female, methods = c("logistic", "logistic"), score = "logistic", extra = list()),
    list(response = "count", data = simulate_count(n = 500, seed = 812),
         outcome = y ~ age + female + offset(log(years)), methods = c("poisson", "poisson"), score = "poisson", extra = list()),
    list(response = "survival", data = simulate_survival(n = 500, seed = 813),
         outcome = survival::Surv(y, d) ~ age + female, methods = c("poisson", "fractional_logit"), score = "poisson", extra = list()))
  for (case in cases) {
    extra <- case$extra
    if (case$response == "survival") extra$tau0 <- median(case$data$y) * 1.1
    args <- c(list(response = case$response, data = case$data, outcome_model = case$outcome,
                   outcome_method = case$methods, benefit_method = c("t_learner", "two_regression"),
                   ps_model = trt ~ age, prop.cutoff = c(.5, 1), B = 1, Kfold = 2, seed = 814), extra)
    full <- suppressWarnings(suppressMessages(do.call(catefit, args)))
    same <- suppressWarnings(suppressMessages(do.call(catefit, c(args, list(benefit_model = ~ age + female)))))
    expect_equal(same$score.twoReg, full$score.twoReg, info = case$response)
    reduced <- suppressWarnings(suppressMessages(do.call(catefit, c(args, list(benefit_model = ~ age)))))
    expect_equal(predict(full, method = case$score), predict(reduced, method = case$score), info = case$response)
    expect_false(isTRUE(all.equal(full$score.twoReg, reduced$score.twoReg)), info = case$response)
    expect_identical(names(coef(reduced, method = "two_regression")), c("(Intercept)", "age"), info = case$response)
    if (case$response != "binary") expect_true("female" %in% names(coef(reduced, method = "poisson")))
    intercept <- suppressWarnings(suppressMessages(do.call(catefit, c(args, list(benefit_model = ~ 1)))))
    expect_identical(names(coef(intercept, method = "two_regression")), "(Intercept)", info = case$response)
    expect_error(do.call(catefit, c(args, list(benefit_model = age ~ female))), "one-sided")
    expect_error(do.call(catefit, c(args, list(benefit_model = ~ age - 1))), "intercept")
    expect_warning(
      expect_error(do.call(catefit, c(args, list(benefit_model = ~ age, init.model = case$outcome))),
                   "only benefit_model"),
      "init.model is deprecated")
  }
})

test_that("benefit model matches the legacy separate initial model", {
  d <- simulate_count(n = 500, seed = 821)
  args <- list(response = "count", data = d, benefit_method = "contrastReg", ps_model = trt ~ age,
               initial.predictor.method = "poisson", prop.cutoff = c(.5, 1), B = 1, Kfold = 2, seed = 822)
  legacy <- suppressWarnings(do.call(catefit, c(args, list(outcome_model = y ~ age + offset(log(years)),
                                                           init.model = y ~ age + female))))
  modern <- suppressWarnings(do.call(catefit, c(args, list(outcome_model = y ~ age + female + offset(log(years)),
                                                           benefit_model = ~ age))))
  expect_equal(legacy$score.contrastReg, modern$score.contrastReg)
  expect_equal(coef(legacy, method = "contrast_regression"), coef(modern, method = "contrast_regression"))
  expect_equal(confint(legacy), confint(modern))
})

test_that("intercept-only contrast scores work with both numerical backends", {
  for (response in c("count", "survival")) {
    d <- if (response == "count") simulate_count(n = 450, seed = 831) else simulate_survival(n = 450, seed = 831)
    model <- if (response == "count") y ~ age + offset(log(years)) else survival::Surv(y, d) ~ age
    args <- list(response = response, data = d, outcome_model = model,
                 outcome_method = if (response == "count") "poisson" else "fractional_logit",
                 benefit_model = ~ 1, benefit_method = "contrast_regression", ps_model = trt ~ age,
                 prop.cutoff = c(.5, 1), B = 1, Kfold = 2, seed = 832)
    if (response == "survival") args$tau0 <- median(d$y) * 1.1
    fit_c <- suppressWarnings(suppressMessages(do.call(catefit, c(args, list(contrast.control = list(solver = "C"))))))
    fit_r <- suppressWarnings(suppressMessages(do.call(catefit, c(args, list(contrast.control = list(solver = "R"))))))
    expect_identical(names(coef(fit_c, method = "contrast_regression")), "(Intercept)")
    expect_equal(fit_c$score.contrastReg, fit_r$score.contrastReg, tolerance = 1e-6)
  }
})

test_that("previous positional calls retain the legacy initial model", {
  d <- simulate_count(n = 450, seed = 841)
  named <- suppressWarnings(catefit("count", d, outcome_model = y ~ age + offset(log(years)),
                                   outcome_method = "poisson", benefit_method = "two_regression",
                                   ps_model = trt ~ age, init.model = y ~ age + female,
                                   prop.cutoff = c(.5, 1), B = 1, Kfold = 2, seed = 842))
  positional <- suppressWarnings(catefit("count", d, y ~ age + offset(log(years)), "poisson",
                                        "two_regression", trt ~ age, "glm", y ~ age + female,
                                        prop.cutoff = c(.5, 1), B = 1, Kfold = 2, seed = 842))
  expect_equal(positional$score.twoReg, named$score.twoReg)
})
