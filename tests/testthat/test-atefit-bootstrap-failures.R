test_that("binary bootstrap failures retain reasons and invalidate inference", {
  calls <- 0L
  local_mocked_bindings(binary_effect = function(...) {
    calls <<- calls + 1L
    if (calls == 2L) stop("test outcome model failure")
    c(risk0 = .3, risk1 = .4, risk.difference = .1 + calls / 1000)
  }, .package = "precmed")
  d <- data.frame(y = rep(c(0, 1), 50), trt = rep(c(0, 0, 1, 1), 25))
  expect_warning(fit <- atefit("binary", d, y ~ 1, trt ~ 1, n.boot = 3, seed = 123),
                 "1 bootstrap replicate.*test outcome model failure")
  expect_true(is.na(fit$risk.difference$SE))
  expect_true(all(is.na(confint(fit))))
  expect_equal(unique(fit$bootstrap.diagnostics$iteration), 1L)
  expect_true(all(fit$bootstrap.diagnostics$reason == "test outcome model failure"))
  expect_equal(nrow(diagnostics(fit)), 3L)
})

test_that("count undefined and failed bootstraps never yield partial inference", {
  calls <- 0L
  local_mocked_bindings(drcount = function(...) {
    calls <<- calls + 1L
    if (calls == 2L) return(list(log.rate.ratio = NA_real_, rate0 = 0, rate1 = 0))
    if (calls == 3L) stop("test count model failure")
    list(log.rate.ratio = .1, rate0 = 1, rate1 = exp(.1))
  }, .package = "precmed")
  expect_warning(fit <- atefitcount(data = countExample, outcome_model = y ~ age + offset(log(years)),
    ps_model = trt ~ age, n.boot = 3, seed = 42), "2 bootstrap replicate")
  expect_true(is.na(fit$log.rate.ratio$SE))
  expect_true(all(is.na(confint(fit))))
  issues <- fit$bootstrap.diagnostics
  expect_true(any(issues$iteration == 1 & issues$measure == "log.rate.ratio" & grepl("undefined", issues$reason)))
  expect_true(any(issues$iteration == 2 & grepl("test count model failure", issues$reason)))
  expect_equal(sum(diagnostics(fit)$split == "bootstrap", na.rm = TRUE), nrow(issues))
})

test_that("survival failures are tracked per measure and preserve unaffected inference", {
  calls <- 0L
  local_mocked_bindings(drsurv = function(...) {
    calls <<- calls + 1L
    list(rmst1 = 5 + calls / 10, rmst0 = 4 + calls / 10,
         log.rmtl.ratio = .2 + calls / 100,
         log.cumulative.hazard.ratio = if (calls == 2L) NA_real_ else .1 + calls / 100,
         risk0 = .3, risk1 = .4, risk.difference = .1)
  }, .package = "precmed")
  d <- simulate_survival(n = 200, seed = 55)
  expect_message(expect_warning(fit <- atefitsurv(d, survival::Surv(y, d) ~ age,
    trt ~ age, tau0 = 500, n.boot = 3, seed = 1), "1 bootstrap replicate"), "recoded")
  expect_true(is.na(fit$log.cumulative.hazard.ratio$SE))
  expect_true(is.finite(fit$log.rmtl.ratio$SE))
  expect_true(is.finite(fit$rmst0$SE))
  expect_equal(fit$bootstrap.diagnostics$measure, "log.cumulative.hazard.ratio")
  expect_match(fit$bootstrap.diagnostics$reason, "undefined")
})

test_that("non-finite bootstrap estimates are diagnosed and counted", {
  values <- cbind(first = c(1, Inf, 2), second = c(1, 2, NaN))
  issues <- precmed:::ate_bootstrap_diagnostics(values, rep(NA_character_, 3))
  expect_equal(issues$iteration, c(2L, 3L))
  expect_warning(precmed:::ate_bootstrap_warn(issues), "2 bootstrap replicate")
})

test_that("survival fitting errors invalidate each affected measure", {
  calls <- 0L
  local_mocked_bindings(drsurv = function(...) {
    calls <<- calls + 1L
    if (calls == 2L) stop("test survival model failure")
    list(rmst1 = 5, rmst0 = 4, log.rmtl.ratio = .2,
         log.cumulative.hazard.ratio = .1, risk0 = .3, risk1 = .4, risk.difference = .1)
  }, .package = "precmed")
  d <- simulate_survival(n = 200, seed = 55)
  expect_message(expect_warning(fit <- atefitsurv(d, survival::Surv(y, d) ~ age,
    trt ~ age, tau0 = 500, n.boot = 3, seed = 1), "test survival model failure"), "recoded")
  expect_equal(nrow(fit$bootstrap.diagnostics), 7L)
  expect_true(all(fit$bootstrap.diagnostics$reason == "test survival model failure"))
  expect_true(all(is.na(confint(fit))))
  expect_true(is.na(fit$rmst0$SE))
})

test_that("single-arm resamples are reported rather than silently discarded", {
  d <- data.frame(trt = c(0, 0, 1, 1), y = c(0, 1, 0, 1))
  expect_warning(fit <- atefit("binary", d, y ~ 1, trt ~ 1, n.boot = 20, seed = 2),
                 "Both treatment groups are required")
  expect_true(any(grepl("Both treatment groups are required", fit$bootstrap.diagnostics$reason)))
  expect_true(is.na(fit$risk.difference$SE))
})
