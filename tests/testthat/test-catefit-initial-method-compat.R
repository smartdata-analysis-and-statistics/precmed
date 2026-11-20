test_that("initial predictor method is hidden from CATE fitting signatures", {
  for (fun in list(catefit, catefitbinary, catefitcount, catefitsurv, catefitmean))
    expect_false("initial.predictor.method" %in% names(formals(fun)))
  expect_true("initial.predictor.method" %in% names(formals(catecv)))
})

test_that("legacy count initial method is preserved in named and positional calls", {
  d <- simulate_count(n = 400, seed = 971)
  outcome <- y ~ age + offset(log(years))
  common <- list(prop.cutoff = c(.8, 1), B = 1, Kfold = 2, seed = 972)
  modern <- suppressWarnings(do.call(catefit, c(list(response = "count", data = d,
    outcome_model = outcome, outcome_method = "poisson", benefit_method = "two_regression",
    ps_model = trt ~ age), common)))
  old_named <- suppressWarnings(do.call(catefit, c(list(response = "count", data = d,
    score.method = "twoReg", cate.model = outcome, ps.model = trt ~ age,
    initial.predictor.method = "poisson"), common)))
  old_positional <- suppressWarnings(do.call(catefit, c(list("count", d, "twoReg", outcome,
    trt ~ age, "glm", NULL, "poisson"), common)))
  old_direct <- suppressWarnings(do.call(catefitcount, c(list(d, "twoReg", outcome,
    trt ~ age, "glm", "poisson"), common)))
  expect_equal(old_named$score.twoReg, modern$score.twoReg)
  expect_equal(old_positional$score.twoReg, modern$score.twoReg)
  expect_equal(old_direct$score.twoReg, modern$score.twoReg)
  expect_error(suppressWarnings(do.call(catefitcount, c(list(data = d, outcome_model = outcome,
    outcome_method = "poisson", benefit_method = "two_regression", ps_model = trt ~ age,
    initial.predictor.method = "gam"), common))), "only outcome_method")
})

test_that("legacy survival initial method remains supported in outcome-specific calls", {
  d <- simulate_survival(n = 400, seed = 973)
  outcome <- survival::Surv(y, d) ~ age
  common <- list(tau0 = median(d$y) * 1.1, prop.cutoff = c(.8, 1), B = 1, Kfold = 2, seed = 974)
  modern <- suppressWarnings(suppressMessages(do.call(catefitsurv, c(list(data = d,
    outcome_model = outcome, outcome_method = "fractional_logit", benefit_method = "two_regression",
    ps_model = trt ~ age), common))))
  old <- suppressWarnings(suppressMessages(do.call(catefitsurv, c(list(d, "twoReg", outcome,
    trt ~ age, "glm", "fractional_logit"), common))))
  expect_equal(old$score.twoReg, modern$score.twoReg)
})

test_that("positional initial methods conflict with a supplied outcome method", {
  d <- simulate_count(n = 400, seed = 975)
  outcome <- y ~ age + offset(log(years))
  common <- list(prop.cutoff = c(.8, 1), B = 1, Kfold = 2, seed = 976)
  expect_error(suppressWarnings(do.call(catefitcount, c(list(d, outcome, "poisson", NULL,
    "two_regression", trt ~ age, "glm", "poisson"), common))), "only outcome_method")
  expect_error(suppressWarnings(do.call(catefit, c(list("count", d, outcome, "poisson", NULL,
    "two_regression", trt ~ age, "glm", "poisson"), common))), "only outcome_method")
})

test_that("the initial predictor alias warns and conflicts with any explicit outcome method", {
  d <- simulate_count(n = 400, seed = 981)
  args <- list(data = d, outcome_model = y ~ age + offset(log(years)),
               benefit_method = "twoReg", ps_model = trt ~ age, initial.predictor.method = "poisson",
               prop.cutoff = c(.8, 1), B = 1, Kfold = 2, seed = 982)
  expect_warning(direct <- do.call(catefitcount, args), "initial.predictor.method is deprecated")
  expect_warning(generic <- do.call(catefit, c(list(response = "count"), args)), "initial.predictor.method is deprecated")
  expect_equal(direct$score.twoReg, generic$score.twoReg)
  for (method in list("poisson", "boosting", NULL)) {
    supplied <- c(args, list(outcome_method = method))
    for (fun in list(catefitcount, function(...) catefit("count", ...))) {
      expect_warning(expect_error(do.call(fun, supplied), "Specify only outcome_method"),
                     "initial.predictor.method is deprecated")
    }
  }
  args$initial.predictor.method <- NULL
  args["initial.predictor.method"] <- list(NULL)
  # Even an explicitly supplied NULL is identified as a legacy argument.
  args$benefit_method <- "poisson"
  expect_warning(do.call(catefitcount, args), "initial.predictor.method is deprecated")
})



test_that("deprecated initial method works with modern benefit-method names", {
  d <- simulate_count(n = 450, seed = 984)
  common <- list(data = d, outcome_model = y ~ age + offset(log(years)), ps_model = trt ~ age,
                 benefit_method = c("t_learner", "two_regression", "contrast_regression"),
                 prop.cutoff = c(.8, 1), B = 1, Kfold = 2, seed = 985)
  modern <- suppressWarnings(do.call(catefitcount, c(common, list(outcome_method = rep("poisson", 3)))))
  expect_warning(legacy <- do.call(catefitcount, c(common, list(initial.predictor.method = "poisson"))),
                 "initial.predictor.method is deprecated")
  expect_warning(generic <- do.call(catefit, c(list(response = "count"), common,
                                               list(initial.predictor.method = "poisson"))),
                 "initial.predictor.method is deprecated")
  expect_equal(legacy$score.twoReg, modern$score.twoReg)
  expect_equal(legacy$score.contrastReg, modern$score.contrastReg)
  expect_equal(generic$score.twoReg, modern$score.twoReg)
  expect_equal(generic$score.contrastReg, modern$score.contrastReg)
})

test_that("binary legacy initial method also resolves modern two-regression", {
  d <- simulate_binary(n = 350, seed = 986)
  common <- list(data = d, outcome_model = y ~ age, benefit_method = "two_regression",
                 ps_model = trt ~ age, prop.cutoff = 1, B = 1, Kfold = 2, seed = 987)
  modern <- do.call(catefitbinary, c(common, list(outcome_method = "logistic")))
  expect_warning(old <- do.call(catefitbinary, c(common, list(initial.predictor.method = "logistic"))),
                 "initial.predictor.method is deprecated")
  expect_warning(generic <- do.call(catefit, c(list(response = "binary"), common,
                                               list(initial.predictor.method = "logistic"))),
                 "initial.predictor.method is deprecated")
  expect_equal(old$score.twoReg, modern$score.twoReg)
  expect_equal(generic$score.twoReg, modern$score.twoReg)
})
