test_that("negative binomial outcome model agrees with MASS and preserves exposure", {
  set.seed(221)
  n <- 700
  d <- data.frame(trt = rep(0:1, each = n / 2), x = rnorm(n),
                  years = runif(n, .4, 2))
  mu <- d$years * exp(.3 + .25 * d$x + .3 * d$trt)
  d$y <- rnbinom(n, mu = mu, size = 1.1)
  design <- model.matrix(~ x, d)
  for (arm in 0:1) {
    selected <- d$trt == arm
    reference <- MASS::glm.nb(y ~ x + offset(log(years)), data = d[selected, ])
    actual <- precmed:::fit_ate_count_outcome(design[selected, , drop = FALSE],
                                               d$y[selected], log(d$years[selected]),
                                               "negative_binomial")
    expect_equal(actual, unname(coef(reference)), tolerance = 1e-6)
  }
  args <- list(response = "count", data = d, outcome_model = y ~ x + offset(log(years)),
               ps_model = trt ~ x, n.boot = 5, seed = 222)
  fit <- do.call(atefit, c(args, list(outcome_method = "negative_binomial")))
  alias <- do.call(atefit, c(args, list(outcome_method = "negbin")))
  expect_equal(fit$log.rate.ratio, alias$log.rate.ratio)
  expect_equal(fit$trt.boot, alias$trt.boot)
  expect_true(all(is.finite(fit$trt.boot)))
  f0 <- MASS::glm.nb(y ~ x + offset(log(years)), data = d, subset = trt == 0)
  f1 <- MASS::glm.nb(y ~ x + offset(log(years)), data = d, subset = trt == 1)
  m0 <- exp(drop(design %*% coef(f0))); m1 <- exp(drop(design %*% coef(f1)))
  p <- precmed:::glm.simplereg.ps(d$trt, model.matrix(~ x, d), minPS = .01, maxPS = .99)
  rates <- c(mean(m0 + (1-d$trt) * (d$y / d$years - m0) / (1-p)),
             mean(m1 + d$trt * (d$y / d$years - m1) / p))
  expect_equal(c(fit$rate0$estimate, fit$rate1$estimate), rates, tolerance = 1e-6)
  pois <- do.call(atefit, args)
  expect_equal(pois$trt.boot, do.call(atefit, c(args, list(outcome_method = "poisson")))$trt.boot)
})

test_that("outcome methods are validated for each response", {
  d <- simulate_binary(n = 250, seed = 223)
  args <- list(data = d, outcome_model = y ~ age, ps_model = trt ~ age, n.boot = 2, seed = 1)
  expect_error(do.call(atefit, c(list(response = "binary"), args, list(outcome_method = "negbin"))), "logistic")
  expect_error(do.call(atefitbinary, c(args, list(outcome_method = NA_character_))), "logistic")
  dc <- simulate_count(n = 250, seed = 224)
  expect_error(atefitcount(dc, y ~ age, trt ~ age, outcome_method = "logistic"), "poisson")
  ds <- simulate_survival(n = 250, seed = 225)
  expect_error(atefitsurv(ds, survival::Surv(y, d) ~ age, trt ~ age, outcome_method = "logistic"), "coxph")
  expect_error(atefit("count", dc, y ~ age, trt ~ age, outcome_method = c("poisson", "negbin")), "outcome_method")
})

test_that("negative binomial fitting errors are recorded for bootstrap samples", {
  d <- simulate_count(n = 250, seed = 226)
  original <- precmed:::fit_ate_count_outcome
  calls <- 0L
  local_mocked_bindings(fit_ate_count_outcome = function(design, y, time, outcome_method) {
    expect_identical(outcome_method, "negative_binomial")
    calls <<- calls + 1L
    if (calls == 3L) stop("test negative binomial failure")
    # Keep this test focused on propagation and error handling.
    original(design, y, time, "poisson")
  }, .package = "precmed")
  expect_warning(fit <- atefitcount(d, y ~ age + offset(log(years)), trt ~ age,
                                    outcome_method = "negbin", n.boot = 3, seed = 227), "bootstrap")
  expect_true(any(grepl("test negative binomial failure", fit$bootstrap.diagnostics$reason)))
  expect_true(is.na(fit$log.rate.ratio$SE))
})


test_that("explicit logistic and Cox methods match the defaults", {
  b <- simulate_binary(n = 300, seed = 228)
  args <- list(response = "binary", data = b, outcome_model = y ~ age,
               ps_model = trt ~ age, n.boot = 3, seed = 229)
  expect_equal(do.call(atefit, args), do.call(atefit, c(args, list(outcome_method = "logistic"))))
  d <- simulate_survival(n = 350, seed = 230)
  args <- list(response = "survival", data = d, outcome_model = survival::Surv(y, d) ~ age,
               ps_model = trt ~ age, tau0 = median(d$y) * 1.1, n.boot = 3, seed = 231)
  implicit <- suppressWarnings(do.call(atefit, args))
  explicit <- suppressWarnings(do.call(atefit, c(args, list(outcome_method = "coxph"))))
  expect_equal(implicit, explicit)
})
