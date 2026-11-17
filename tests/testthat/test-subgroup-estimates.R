test_that("subgroup extraction validates methods and supports older fits", {
  fit <- structure(list(response = "survival", score.contrastReg = c(1, 2),
    ate.contrastReg = c(prop0.5 = .8, prop1 = .9),
    chr.contrastReg = c(prop0.5 = .7, prop1 = .8),
    rd.contrastReg = c(prop0.5 = -.2, prop1 = -.1)), class = "catefit")
  expect_equal(subgroup_estimates(fit), fit$ate.contrastReg)
  expect_equal(subgroup_estimates(fit, measure = "cumulative.hr"), fit$chr.contrastReg)
  expect_equal(subgroup_estimates(fit, measure = "cumulative_hr"), fit$chr.contrastReg)
  expect_equal(subgroup_estimates(fit, measure = "absolute_risk_reduction"), -fit$rd.contrastReg)
  expect_error(subgroup_estimates(fit, measure = "risk0"), "refit")
  expect_error(subgroup_estimates(fit, method = "unknown"), "available scoring")
  expect_error(subgroup_estimates(fit, measure = "unknown"), "measure")
  fit$score.poisson <- c(1, 2)
  expect_error(subgroup_estimates(fit), "Specify one")
  fit$response <- "continuous"
  expect_error(subgroup_estimates(fit), "binary, count or survival")
})

test_that("survival subgroup risks and differences come from the same fits", {
  d <- simulate_survival(n = 600, seed = 81)
  expect_message(fit <- catefitsurv(data = d, score.method = "poisson",
    cate.model = survival::Surv(y, d) ~ age, ps.model = trt ~ age,
    tau0 = 300, prop.cutoff = c(.5, 1), seed = 12, verbose = 0), "recoded")
  expect_equal(predict(fit), fit$score.poisson)
  risk0 <- subgroup_estimates(fit, measure = "risk0")
  risk1 <- subgroup_estimates(fit, measure = "risk1")
  expect_named(risk0, c("prop0.5", "prop1"))
  expect_equal(subgroup_estimates(fit, measure = "risk.difference"), risk1 - risk0)
  expect_equal(subgroup_estimates(fit, measure = "absolute_risk_reduction"), risk0 - risk1)
  expect_true(all(is.finite(subgroup_estimates(fit, measure = "rmst0"))))
  expect_true(all(is.finite(subgroup_estimates(fit, measure = "cumulative_hr"))))
})


test_that("count subgroup estimates expose stored rate ratios", {
  fit <- structure(list(response = "count", score.poisson = c(.1, .2),
    ate.poisson = c(prop0.5 = .8, prop1 = 1.1)), class = "catefit")
  expect_identical(subgroup_estimates(fit), fit$ate.poisson)
  expect_identical(subgroup_estimates(fit, method = "poisson", measure = "rate_ratio"), fit$ate.poisson)
  expect_identical(subgroup_estimates(fit, measure = "rate.ratio"), fit$ate.poisson)
  expect_error(subgroup_estimates(fit, measure = "risk0"), "measure")
})


test_that("binary subgroup estimates preserve signed risk differences", {
  set.seed(123)
  d <- data.frame(x = rnorm(300), trt = rep(0:1, each = 150))
  d$y <- rbinom(300, 1, plogis(-.4 + .3 * d$x - .5 * d$trt))
  fit <- catefitbinary(d, "logistic", y ~ x, trt ~ 1,
                       prop.cutoff = c(.5, 1), seed = 19)
  expect_identical(subgroup_estimates(fit), fit$ate.logistic)
  expect_identical(subgroup_estimates(fit, method = "logistic", measure = "risk_difference"), fit$ate.logistic)
  expect_error(subgroup_estimates(fit, measure = "risk0"), "measure")
  expect_equal(names(subgroup_estimates(fit)), c("prop0.5", "prop1"))
})
