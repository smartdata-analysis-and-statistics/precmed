test_that("survival hazard ratio censors follow-up beyond the horizon", {
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
  dat$time_at_horizon <- pmin(dat$y, tau)
  dat$event_at_horizon <- dat$d * (dat$y <= tau)
  expected <- survival::coxph(
    survival::Surv(time_at_horizon, event_at_horizon) ~
      trt + age + female + previous_cost + previous_number_relapses,
    data = dat
  )
  expect_equal(fit$log.hazard.ratio$estimate, unname(coef(expected)[1]), tolerance = 1e-8)
  set.seed(123)
  expected.boot <- replicate(2, {
    sample.dat <- dat[sample(nrow(dat), replace = TRUE), ]
    unname(coef(update(expected, data = sample.dat))[1])
  })
  expect_equal(unname(fit$trt.boot[, "log.hazard.ratio"]), expected.boot, tolerance = 1e-8)
})
