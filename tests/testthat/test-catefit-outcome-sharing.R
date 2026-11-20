test_that("count scores with different outcome methods agree with separate fits", {
  d <- simulate_count(n = 500, seed = 991)
  args <- list(data = d, outcome_model = y ~ age + offset(log(years)), ps_model = trt ~ age,
               benefit_model = ~ age, prop.cutoff = c(.8, 1), B = 2, Kfold = 2, seed = 992)
  combined <- suppressWarnings(do.call(catefitcount, c(args, list(
    benefit_method = c("two_regression", "contrast_regression"), outcome_method = c("poisson", "gam")))))
  two <- suppressWarnings(do.call(catefitcount, c(args, list(benefit_method = "two_regression", outcome_method = "poisson"))))
  contrast <- suppressWarnings(do.call(catefitcount, c(args, list(benefit_method = "contrast_regression", outcome_method = "gam"))))
  expect_equal(combined$score.twoReg, two$score.twoReg)
  expect_equal(combined$score.contrastReg, contrast$score.contrastReg, tolerance = 1e-8)
  expect_identical(combined$outcome_methods, c(two_regression = "poisson", contrast_regression = "gam"))
  generic <- suppressWarnings(do.call(catefit, c(list(response = "count"), args, list(
    benefit_method = c("two_regression", "contrast_regression"), outcome_method = c("poisson", "gam")))))
  expect_equal(combined, generic)
})

test_that("cross-fitting shares treatment models and matching outcome predictions", {
  d <- simulate_count(n = 400, seed = 993)
  x <- model.matrix(~ age, d)[, -1, drop = FALSE]
  z <- model.matrix(~ age, d)
  trt <- as.integer(d$trt == "drug1")
  original_ps <- precmed:::glm.simplereg.ps
  count <- 0L
  local_mocked_bindings(glm.simplereg.ps = function(...) { count <<- count + 1L; original_ps(...) }, .package = "precmed")
  args <- list(y = d$y, trt = trt, x.cate = x, x.ps = z, time = log(d$years),
               score.method = c("twoReg", "contrastReg"), B = 1, Kfold = 2,
               initial.predictor.method = "poisson", contrast.control = list(solver = "R"))
  set.seed(994)
  shared <- suppressWarnings(do.call(precmed:::intxcount, c(args, list(initial.methods = c(twoReg = "poisson", contrastReg = "poisson")))))
  expect_equal(count, 2L)
  count <- 0L
  set.seed(994)
  separate <- suppressWarnings(do.call(precmed:::intxcount, c(args, list(initial.methods = c(twoReg = "poisson", contrastReg = "gam")))))
  expect_equal(count, 2L)
  expect_equal(shared$result.twoReg, separate$result.twoReg)
  expect_false(isTRUE(all.equal(shared$result.contrastReg$delta.contrastReg, separate$result.contrastReg$delta.contrastReg)))
})

test_that("survival scores can use different initial outcome methods", {
  d <- simulate_survival(n = 400, seed = 995)
  args <- list(data = d, outcome_model = survival::Surv(y, d) ~ age + female,
               ps_model = trt ~ age, tau0 = median(d$y) * 1.1, prop.cutoff = c(.8, 1),
               B = 1, Kfold = 2, n.trees.rf = 40, seed = 996)
  combined <- suppressWarnings(suppressMessages(do.call(catefitsurv, c(args, list(
    benefit_method = c("two_regression", "contrast_regression"), outcome_method = c("fractional_logit", "random_forest"))))))
  two <- suppressWarnings(suppressMessages(do.call(catefitsurv, c(args, list(
    benefit_method = "two_regression", outcome_method = "fractional_logit")))))
  expect_equal(combined$score.twoReg, two$score.twoReg)
  expect_true(all(is.finite(combined$score.contrastReg)))
  expect_identical(combined$outcome_methods, c(two_regression = "fractional_logit", contrast_regression = "random_forest"))
})

test_that("matching boosting outcome models are fitted once per arm and fold", {
  d <- simulate_count(n = 400, seed = 997)
  x <- model.matrix(~ age + female, d)[, -1, drop = FALSE]
  z <- model.matrix(~ age, d)
  original_fit <- precmed:::fit_boosting
  fit_count <- 0L
  local_mocked_bindings(fit_boosting = function(...) { fit_count <<- fit_count + 1L; original_fit(...) }, .package = "precmed")
  set.seed(998)
  suppressWarnings(precmed:::intxcount(y = d$y, trt = as.integer(d$trt == "drug1"), x.cate = x, x.ps = z,
    time = log(d$years), score.method = c("twoReg", "contrastReg"), B = 1, Kfold = 2,
    initial.predictor.method = "boosting", initial.methods = c(twoReg = "boosting", contrastReg = "boosting"),
    n.trees.boosting = 20, contrast.control = list(solver = "R"), n.cores = 1))
  expect_equal(fit_count, 4L)
})
