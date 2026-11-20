test_that("survival T-learners predict arm-specific restricted mean time lost", {
  d <- simulate_survival(n = 500, seed = 941)
  tau <- as.numeric(quantile(d$y, .75))
  args <- list(response = "survival", data = d, outcome_model = survival::Surv(y, d) ~ age + female,
               outcome_method = c("coxph", "weibull", "log_normal"),
               benefit_method = rep("t_learner", 3), ps_model = trt ~ age,
               tau0 = tau, prop.cutoff = c(.8, 1), seed = 942)
  fit <- suppressWarnings(suppressMessages(do.call(catefit, args)))
  x <- stats::model.matrix(~ age + female, d)[, -1, drop = FALSE]
  arm <- as.integer(d$trt == "drug1")
  for (method in args$outcome_method) {
    calc <- function(k) {
      if (method == "coxph") precmed:::cox.rmst(d$y[arm == k], d$d[arm == k], x[arm == k, , drop = FALSE], x, tau)
      else precmed:::parametric_survival_outcome(d$y[arm == k], d$d[arm == k], x[arm == k, , drop = FALSE], x, tau, method)
    }
    expected <- log(tau - calc(1)) - log(tau - calc(0))
    expect_equal(predict(fit, method = method), expected, tolerance = 1e-9)
    expect_length(subgroup_estimates(fit, method = method, measure = "rmtl"), 2L)
  }
  expect_s3_class(plot(fit, type = "density"), "ggplot")
  expect_output(print(fit), "log_normal")
  expect_error(predict(fit, method = "t_learner"), "Select one T-learner")
  args$outcome_method <- "lognormal"; args$benefit_method <- "t_learner"
  alias <- suppressWarnings(suppressMessages(do.call(catefit, args)))
  expect_equal(alias$score.log_normal, fit$score.log_normal)
  expect_equal(predict(alias, method = "t_learner"), alias$score.log_normal)
})

test_that("survival T-learners handle intercept-only and reject unsupported combinations", {
  d <- simulate_survival(n = 400, seed = 943)
  tau <- median(d$y) * 1.1
  args <- list(data = d, outcome_model = survival::Surv(y, d) ~ 1, outcome_method = "coxph",
               benefit_method = "t_learner", ps_model = trt ~ 1,
               tau0 = tau, prop.cutoff = c(.8, 1), seed = 944)
  direct <- suppressWarnings(suppressMessages(do.call(catefitsurv, args)))
  generic <- suppressWarnings(suppressMessages(do.call(catefit, c(list(response = "survival"), args))))
  expect_equal(generic$score.coxph, direct$score.coxph)
  expect_true(all(is.finite(direct$score.coxph)))
  expect_true(all(is.finite(direct$ate.coxph)))
  expect_equal(length(unique(direct$score.coxph)), 1L)
  args$outcome_method <- "weibull"
  expect_true(all(is.finite(suppressWarnings(suppressMessages(do.call(catefitsurv, args)))$score.weibull)))
  args$benefit_method <- "contrast_regression"
  expect_error(do.call(catefitsurv, args), "not supported")
  expect_error(do.call(catecvsurv, args), "not supported")
  args$benefit_method <- "t_learner"; args$outcome_method <- "coxph"
  one_arm <- d[d$trt == d$trt[1], , drop = FALSE]
  args$data <- one_arm
  expect_error(suppressMessages(do.call(catefitsurv, args)), "treatment|Treatment|arm|group|trt")
  bad <- transform(d, d = 0)
  args$data <- bad; args$benefit_method <- "t_learner"; args$outcome_method <- "coxph"
  expect_error(suppressMessages(do.call(catefitsurv, args)), "at least one event")
})

test_that("intercept-only Cox integration agrees with a survival reference", {
  y <- c(1, 2, 2, 3, 4, 5); d <- c(1, 1, 1, 0, 1, 0); tau <- 3.5
  x <- matrix(numeric(), nrow = length(y), ncol = 0)
  fit <- survival::survfit(survival::Surv(y, d) ~ 1, stype = 2)
  events <- fit$time[fit$time <= tau]
  survival <- fit$surv[fit$time <= tau]
  expected <- sum(c(1, survival) * diff(c(0, events, tau)))
  expect_equal(precmed:::cox.rmst(y, d, x, x, tau), rep(expected, length(y)), tolerance = 1e-10)
})

test_that("new survival T-learners can be combined with doubly robust scores", {
  d <- simulate_survival(n = 450, seed = 945)
  args <- list(response = "survival", data = d, outcome_model = survival::Surv(y, d) ~ age + female,
               outcome_method = c("coxph", "fractional_logit"),
               benefit_method = c("t_learner", "two_regression"), benefit_model = ~ age,
               ps_model = trt ~ age, tau0 = median(d$y) * 1.1,
               prop.cutoff = c(.8, 1), B = 1, Kfold = 2, seed = 946)
  fit <- suppressWarnings(suppressMessages(do.call(catefit, args)))
  expect_true(all(is.finite(fit$score.coxph)))
  expect_true(all(is.finite(fit$score.twoReg)))
  expect_identical(names(coef(fit, method = "two_regression")), c("(Intercept)", "age"))
  expect_equal(predict(fit, method = "t_learner"), fit$score.coxph)
})
