test_that("count ATE fits separate outcome models", {
  d <- simulate_count(n = 800, seed = 104)
  args <- list(response = "count", data = d,
               outcome_model = y ~ age + female + offset(log(years)),
               ps.model = trt ~ age, n.boot = 5, seed = 105)
  default <- suppressMessages(do.call(atefit, args))
  arm <- as.numeric(d$trt == "drug1")
  design <- model.matrix(~ age + female, d)
  fits <- lapply(0:1, function(a) glm(y ~ age + female + offset(log(years)),
                                     family = poisson(), data = d, subset = arm == a))
  m0 <- exp(drop(design %*% coef(fits[[1]])))
  m1 <- exp(drop(design %*% coef(fits[[2]])))
  ps <- precmed:::glm.simplereg.ps(arm, model.matrix(~ age, d), minPS = .01, maxPS = .99)
  rates <- c(mean(m0 + (1 - arm) * (d$y / d$years - m0) / (1 - ps)),
             mean(m1 + arm * (d$y / d$years - m1) / ps))
  expect_equal(c(default$rate0$estimate, default$rate1$estimate), unname(rates), tolerance = 1e-7)
})

test_that("count outcome models handle small and missing treatment arms", {
  d <- data.frame(y = rep(c(1, 2), 11), trt = rep(0:1, each = 11), years = 1)
  expect_no_error(atefit("count", d[-22, ], y ~ 1 + offset(log(years)), trt ~ 1,
                         n.boot = 2, seed = 6))
  one_arm <- transform(d, trt = 0)
  expect_error(atefit("count", one_arm, y ~ 1 + offset(log(years)), trt ~ 1,
                      n.boot = 2), "two|groups|treatment", ignore.case = TRUE)
})
