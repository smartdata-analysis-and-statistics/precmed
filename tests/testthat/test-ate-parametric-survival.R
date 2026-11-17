test_that("parametric survival predictions agree with numerical integration", {
  set.seed(401)
  n <- 600
  x <- rnorm(n)
  for (method in c("weibull", "log_normal")) {
    event <- if (method == "weibull") rweibull(n, shape = 1.4, scale = exp(1 + .25 * x)) else
      rlnorm(n, meanlog = 1 + .25 * x, sdlog = .7)
    censor <- rexp(n, rate = .12)
    y <- pmin(event, censor); d <- as.numeric(event <= censor)
    tau <- 4
    newx <- matrix(c(-1, 0, 1), ncol = 1)
    result <- precmed:::parametric_survival_outcome(y, d, matrix(x, ncol = 1), newx,
                                                   tau, method, TRUE)
    reference <- survival::survreg(survival::Surv(y, d) ~ x,
                                   dist = if (method == "weibull") "weibull" else "lognormal")
    eta <- unname(coef(reference)[1] + coef(reference)[2] * newx[, 1])
    surv <- lapply(eta, function(e) {
      if (method == "weibull") function(t) pweibull(t, shape = 1 / reference$scale, scale = exp(e), lower.tail = FALSE)
      else function(t) plnorm(t, meanlog = e, sdlog = reference$scale, lower.tail = FALSE)
    })
    expect_equal(result$rmst, vapply(surv, function(f) integrate(f, 0, tau)$value, numeric(1)), tolerance = 1e-6, info = method)
    expect_equal(result$survival, vapply(surv, function(f) f(tau), numeric(1)), tolerance = 1e-7, info = method)
    expect_equal(precmed:::parametric_survival_outcome(y, d, matrix(x, ncol = 1), newx, tau, method), result$rmst)
  }
})

test_that("ATE survival parametric methods fit and bootstrap arm-specific models", {
  d <- simulate_survival(n = 500, seed = 402)
  tau <- unname(quantile(d$y, .7))
  args <- list(response = "survival", data = d,
               outcome_model = survival::Surv(y, d) ~ age + female,
               ps.model = trt ~ age, tau0 = tau, n.boot = 4, seed = 403)
  for (method in c("weibull", "log_normal")) {
    fit <- suppressMessages(do.call(atefit, c(args, list(outcome_method = method))))
    expect_true(all(is.finite(fit$trt.boot)), info = method)
    expect_equal(fit$risk.difference$estimate, fit$risk1$estimate - fit$risk0$estimate, info = method)
    expect_equal(fit$risk.difference$SE, sd(fit$trt.boot[, "risk.difference"]), info = method)
    expect_true(all(c(fit$rmst0$estimate, fit$rmst1$estimate) >= 0 &
                      c(fit$rmst0$estimate, fit$rmst1$estimate) <= tau), info = method)
  }
  a <- suppressMessages(do.call(atefit, c(args, list(outcome_method = "log_normal"))))
  b <- suppressMessages(do.call(atefit, c(args, list(outcome_method = "lognormal"))))
  expect_equal(a$trt.boot, b$trt.boot)
})

test_that("parametric survival models handle intercepts and invalid arms", {
  set.seed(404)
  y <- rweibull(100, 1.5, scale = 4)
  x <- matrix(numeric(), 100, 0)
  for (method in c("weibull", "log_normal")) {
    expect_length(precmed:::parametric_survival_outcome(y, rep(1, 100), x, x[1:3, , drop = FALSE], 3, method), 3)
    expect_error(precmed:::parametric_survival_outcome(y, rep(0, 100), x, x, 3, method), "at least one event")
    expect_error(precmed:::parametric_survival_outcome(replace(y, 1, 0), rep(1, 100), x, x, 3, method), "positive observed times")
  }
})

test_that("parametric model bootstrap failures are recorded", {
  d <- simulate_survival(n = 350, seed = 405)
  original <- precmed:::parametric_survival_outcome
  calls <- 0L
  local_mocked_bindings(parametric_survival_outcome = function(...) {
    calls <<- calls + 1L
    if (calls == 3L) stop("test parametric fit failure")
    original(...)
  }, .package = "precmed")
  expect_warning(fit <- suppressMessages(atefit("survival", d, survival::Surv(y, d) ~ age,
                                 trt ~ age, tau0 = median(d$y) * 1.1,
                                 outcome_method = "weibull", n.boot = 3, seed = 406)), "bootstrap")
  expect_true(any(grepl("test parametric fit failure", fit$bootstrap.diagnostics$reason)))
  expect_true(is.na(fit$rmst0$SE))
})

test_that("separate Weibull shapes can represent non-proportional treatment hazards", {
  set.seed(407)
  n <- 700
  time0 <- rweibull(n, shape = .8, scale = 4)
  time1 <- rweibull(n, shape = 1.8, scale = 4)
  x <- matrix(numeric(), n, 0)
  prediction <- function(time, tau) {
    precmed:::parametric_survival_outcome(time, rep(1, n), x, x[1, , drop = FALSE],
                                         tau, "weibull", TRUE)$survival
  }
  hr_at_horizon <- function(tau) -log(prediction(time1, tau)) / -log(prediction(time0, tau))
  expect_gt(abs(log(hr_at_horizon(2) / hr_at_horizon(6))), .5)
  for (method in c("weibull", "log_normal")) {
    short <- precmed:::parametric_survival_outcome(time0, rep(1, n), x, x[1, , drop = FALSE], .001, method)
    expect_true(is.finite(short) && short >= 0 && short <= .001)
    expect_error(precmed:::parametric_survival_outcome(time0, rep(1, n), x, x[1, , drop = FALSE], 0, method), "tau0")
  }
})
