test_that("binary ATE uses separate outcome models", {
  d <- simulate_binary(n = 800, seed = 101)
  dat <- precmed:::binary_data(d, y ~ age + female, trt ~ age)
  separate <- precmed:::binary_predictions(dat$y, dat$trt, dat$x, "twoReg")
  args <- list(response = "binary", data = d, outcome_model = y ~ age + female,
               ps.model = trt ~ age, n.boot = 5, seed = 102)
  fit <- do.call(atefit, args)
  p <- precmed:::glm.simplereg.ps(dat$trt, dat$z, minPS = .01, maxPS = .99)
  r0 <- mean(separate$mu0 + (1 - dat$trt) * (dat$y - separate$mu0) / (1 - p))
  r1 <- mean(separate$mu1 + dat$trt * (dat$y - separate$mu1) / p)
  expect_equal(fit$risk0$estimate, r0)
  expect_equal(fit$risk1$estimate, r1)
})

test_that("ATE interfaces no longer accept interactions", {
  for (name in c("atefit", "atefitbinary", "atefitcount", "atefitsurv", "atefitmean"))
    expect_false("interactions" %in% names(formals(get(name, envir = asNamespace("precmed")))), info = name)
  d <- simulate_binary(n = 100, seed = 103)
  expect_error(atefit("binary", d, y ~ age, trt ~ age, interactions = TRUE), "unused argument")
})
