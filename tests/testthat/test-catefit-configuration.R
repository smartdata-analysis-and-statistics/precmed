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
  expect_identical(names(formals(catefit))[1:8],
                   c("response", "data", "outcome_model", "outcome_method", "benefit_model", "benefit_method", "ps_model", "ps_method"))
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
  expect_error(suppressWarnings(do.call(catefit, c(args, list(benefit_method = "two_regression", outcome_method = "logistic",
                                            initial.predictor.method = "boosting")))), "only outcome_method")
  expect_error(suppressWarnings(do.call(catefit, c(args, list(benefit_method = c("t_learner", "two_regression"),
                                                    outcome_method = c("logistic", "logistic"),
                                                    initial.predictor.method = "logistic")))), "only outcome_method")
  expect_error(do.call(catefit, c(args, list(benefit_method = "contrast_regression", outcome_method = "logistic"))), "not supported")
  expect_error(do.call(catefit, c(args, list(benefit_method = "contrastReg", outcome_method = "logistic"))), "not supported")
  expect_error(do.call(catefit, c(args, list(benefit_method = c("t_learner", "two_regression"), outcome_method = "logistic"))), "one non-missing")
  config <- precmed:::resolve_benefit_config("count", c("two_regression", "contrast_regression"), c("poisson", "boosting"))
  expect_identical(unname(config$initial_methods), c("poisson", "boosting"))
  expect_false(any(c("score.method", "cate.model") %in% names(formals(catefit))))
})

test_that("survival configurations preserve cross-fitted outcome models", {
  d <- simulate_survival(n = 500, seed = 206)
  args <- list(response = "survival", data = d, ps_model = trt ~ age,
               tau0 = median(d$y) * 1.1, prop.cutoff = c(.8, 1), B = 1, Kfold = 2, seed = 207)
  formula <- survival::Surv(y, d) ~ age
  modern <- suppressWarnings(suppressMessages(do.call(catefit, c(args, list(
    outcome_model = formula, benefit_method = c("t_learner", "two_regression", "contrast_regression"),
    outcome_method = c("poisson", "fractional_logit", "fractional_logit"))))))
  old <- suppressWarnings(suppressMessages(do.call(catefit, c(args, list(
    cate.model = formula, score.method = c("poisson", "twoReg", "contrastReg"),
    initial.predictor.method = "fractional_logit")))))
  expect_equal(modern, old)
  aliases <- suppressWarnings(suppressMessages(do.call(catefit, c(args, list(
    outcome_model = formula, benefit_method = c("t_learner", "twoReg", "contrastReg"),
    outcome_method = c("poisson", "fractional_logit", "fractional_logit"))))))
  expect_equal(aliases, modern)
  expect_error(do.call(catefit, c(args, list(outcome_model = formula,
    benefit_method = "two_regression", outcome_method = "poisson"))), "not supported")
})

test_that("survival fractional logit accepts the former label as an alias", {
  d <- simulate_survival(n = 350, seed = 618)
  args <- list(data = d, outcome_model = survival::Surv(y, d) ~ age,
               ps_model = trt ~ age, benefit_method = c("two_regression", "contrast_regression"),
               outcome_method = c("fractional_logit", "fractional_logit"),
               tau0 = median(d$y), prop.cutoff = c(.8, 1), B = 1, Kfold = 2, seed = 619)
  direct <- suppressWarnings(suppressMessages(do.call(catefitsurv, args)))
  generic <- suppressWarnings(suppressMessages(do.call(catefit, c(list(response = "survival"), args))))
  expect_equal(generic, direct)
  expect_true(all(is.finite(direct$score.twoReg)))
  expect_true(all(is.finite(direct$score.contrastReg)))
  args$outcome_method <- c("logistic", "logistic")
  expect_warning(old_direct <- suppressMessages(do.call(catefitsurv, args)), "deprecated")
  expect_warning(old_generic <- suppressMessages(do.call(catefit, c(list(response = "survival"), args))), "deprecated")
  expect_equal(old_direct, direct)
  expect_equal(old_generic, generic)
  old_args <- args
  old_args$outcome_method <- NULL
  old_args$benefit_method <- c("twoReg", "contrastReg")
  expect_warning(old_initial <- suppressMessages(do.call(catefitsurv, c(
    old_args, list(initial.predictor.method = "logistic")))), "deprecated")
  expect_equal(old_initial, direct)
})
