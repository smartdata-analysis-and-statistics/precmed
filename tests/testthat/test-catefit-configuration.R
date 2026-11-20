test_that("paired CATE configurations retain existing score calculations", {
  d <- simulate_count(n = 500, seed = 202)
  args <- list(response = "count", data = d, ps_model = trt ~ age,
               prop.cutoff = 1, B = 1, Kfold = 2, seed = 203)
  outcome <- y ~ age + offset(log(years))
  modern <- suppressWarnings(do.call(catefit, c(args, list(
    outcome_model = outcome, benefit_method = c("t_learner", "two_regression", "contrast_regression"),
    outcome_method = c("poisson", "poisson", "poisson")))))
  old <- suppressWarnings(do.call(catefit, c(args, list(
    cate.model = outcome, score.method = c("poisson", "twoReg", "contrastReg"),
    initial.predictor.method = "poisson"))))
  expect_equal(modern, old)
  aliases <- suppressWarnings(do.call(catefit, c(args, list(
    outcome_model = outcome, benefit_method = c("t_learner", "twoReg", "contrastReg"),
    outcome_method = c("poisson", "poisson", "poisson")))))
  expect_equal(aliases, modern)
  mixed <- suppressWarnings(do.call(catefit, c(args, list(
    outcome_model = outcome, benefit_method = c("t_learner", "twoReg", "contrast_regression"),
    outcome_method = c("poisson", "poisson", "poisson")))))
  expect_equal(mixed, modern)
  expect_error(do.call(catefit, c(args, list(outcome_model = outcome,
    benefit_method = c("twoReg", "two_regression"),
    outcome_method = c("poisson", "poisson")))), "same score")
  expect_s3_class(plot(modern), "ggplot")
  expect_equal(predict(modern, method = "twoReg"), old$score.twoReg)
  expect_equal(predict(modern, method = "two_regression"), old$score.twoReg)
  expect_equal(predict(modern, method = "contrast_regression"), old$score.contrastReg)
  expect_equal(coef(modern), coef(old))
  expect_output(print(modern), "twoReg")
  local_fit <- function() {
    local_data <- d
    local_formula <- outcome
    catefit("count", local_data, benefit_method = "t_learner",
            outcome_method = "poisson", outcome_model = local_formula,
            ps_model = trt ~ age, prop.cutoff = 1, seed = 203)
  }
  expect_equal(local_fit()$score.poisson, modern$score.poisson)
})

test_that("CATE configuration checks supported pairs and compatibility aliases", {
  expect_silent(config <- precmed:::resolve_benefit_config(
    "count", c("twoReg", "contrastReg"), c("poisson", "poisson")))
  expect_identical(config$scores, c("twoReg", "contrastReg"))
  d <- simulate_binary(n = 300, seed = 204)
  expect_identical(names(formals(catefit))[1:7],
                   c("response", "data", "outcome_model", "outcome_method", "benefit_method", "ps_model", "ps_method"))
  args <- list(response = "binary", data = d, outcome_model = y ~ age,
               ps_model = trt ~ age, prop.cutoff = 1, B = 1, Kfold = 2, seed = 205)
  fit <- do.call(catefit, c(args, list(benefit_method = c("t_learner", "two_regression"),
                                      outcome_method = c("logistic", "logistic"))))
  legacy <- suppressWarnings(catefit("binary", d, score.method = c("logistic", "twoReg"),
                                      cate.model = y ~ age, ps_model = trt ~ age,
                                      prop.cutoff = 1, B = 1, Kfold = 2, seed = 205))
  expect_equal(fit, legacy)
  positional <- catefit("binary", d, y ~ age, c("logistic", "logistic"),
                        c("t_learner", "two_regression"), trt ~ age,
                        prop.cutoff = 1, B = 1, Kfold = 2, seed = 205)
  expect_equal(positional, fit)
  previous_positional <- catefit("binary", d, c("logistic", "twoReg"), y ~ age,
                                 trt ~ age, "glm", prop.cutoff = 1,
                                 B = 1, Kfold = 2, seed = 205)
  expect_equal(previous_positional, fit)
  previous_new_order <- catefit("binary", d, c("t_learner", "two_regression"),
                                 y ~ age, trt ~ age,
                                 outcome_method = c("logistic", "logistic"),
                                 prop.cutoff = 1, B = 1, Kfold = 2, seed = 205)
  expect_equal(previous_new_order, fit)
  expect_warning(catefit("binary", d, score.method = "logistic", y ~ age, trt ~ age,
                         prop.cutoff = 1), "score.method is deprecated")
  expect_error(do.call(catefit, c(args, list(benefit_method = "two_regression", outcome_method = "boosting"))), "not supported")
  expect_error(do.call(catefit, c(args, list(benefit_method = "two_regression", outcome_method = "logistic",
                                            initial.predictor.method = "boosting"))), "must agree")
  expect_equal(fit, do.call(catefit, c(args, list(benefit_method = c("t_learner", "two_regression"),
                                                    outcome_method = c("logistic", "logistic"),
                                                    initial.predictor.method = "logistic"))))
  expect_error(do.call(catefit, c(args, list(benefit_method = "contrast_regression", outcome_method = "logistic"))), "not supported")
  expect_error(do.call(catefit, c(args, list(benefit_method = "contrastReg", outcome_method = "logistic"))), "not supported")
  expect_error(do.call(catefit, c(args, list(benefit_method = c("t_learner", "two_regression"), outcome_method = "logistic"))), "one non-missing")
  expect_error(catefit("count", data.frame(), benefit_method = c("two_regression", "contrast_regression"),
                       outcome_method = c("poisson", "boosting")), "same outcome_method")
  expect_false(any(c("score.method", "cate.model") %in% names(formals(catefit))))
})

test_that("survival configurations preserve cross-fitted outcome models", {
  d <- simulate_survival(n = 500, seed = 206)
  args <- list(response = "survival", data = d, ps_model = trt ~ age,
               tau0 = median(d$y) * 1.1, prop.cutoff = c(.8, 1), B = 1, Kfold = 2, seed = 207)
  formula <- survival::Surv(y, d) ~ age
  modern <- suppressWarnings(suppressMessages(do.call(catefit, c(args, list(
    outcome_model = formula, benefit_method = c("t_learner", "two_regression", "contrast_regression"),
    outcome_method = c("poisson", "logistic", "logistic"))))))
  old <- suppressWarnings(suppressMessages(do.call(catefit, c(args, list(
    cate.model = formula, score.method = c("poisson", "twoReg", "contrastReg"),
    initial.predictor.method = "logistic")))))
  expect_equal(modern, old)
  aliases <- suppressWarnings(suppressMessages(do.call(catefit, c(args, list(
    outcome_model = formula, benefit_method = c("t_learner", "twoReg", "contrastReg"),
    outcome_method = c("poisson", "logistic", "logistic"))))))
  expect_equal(aliases, modern)
  expect_error(do.call(catefit, c(args, list(outcome_model = formula,
    benefit_method = "two_regression", outcome_method = "poisson"))), "not supported")
})
