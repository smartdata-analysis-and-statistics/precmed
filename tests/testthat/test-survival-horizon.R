test_that("survival horizon preserves events and removes the conditional HR from ATE output", {
  data(survivalExample, package = "precmed")
  dat <- survivalExample
  tau <- unname(quantile(dat$y, 0.85))
  # Exercise both an event and a censoring observation exactly at the horizon.
  dat$y[1:2] <- tau
  dat$d[1:2] <- c(1, 0)
  expect_warning(fit <- atefitsurv(
    data = dat,
    cate.model = survival::Surv(y, d) ~ age + female + previous_cost + previous_number_relapses,
    ps.model = trt ~ age + previous_treatment,
    tau0 = tau, n.boot = 2, seed = 123, rmst.control = list(solver = "R")
  ), "recoded to 0/1")
  expect_null(fit$log.hazard.ratio)
  expect_false("log.hazard.ratio" %in% colnames(fit$trt.boot))
  expect_true("log.cumulative.hazard.ratio" %in% colnames(fit$trt.boot))
  expect_equal(survival_event_at_horizon(c(tau, tau, tau + 1), c(1, 0, 1), tau),
               c(1, 0, 0))
})
