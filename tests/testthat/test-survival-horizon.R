test_that("survival horizon preserves events and removes the conditional HR from ATE output", {
  data(survivalExample, package = "precmed")
  dat <- survivalExample
  tau <- unname(quantile(dat$y, 0.85))
  # Exercise both an event and a censoring observation exactly at the horizon.
  dat$y[1:2] <- tau
  dat$d[1:2] <- c(1, 0)
  expect_message(fit <- atefitsurv(
    data = dat,
    cate.model = survival::Surv(y, d) ~ age + female + previous_cost + previous_number_relapses,
    ps.model = trt ~ age + previous_treatment,
    tau0 = tau, n.boot = 2, seed = 123, rmst.control = list(solver = "R")
  ), "recoded to 0/1")
  expect_equal(sum(diagnostics(fit, "propensity", "arms")$n), nrow(dat))
  expect_null(fit$log.hazard.ratio)
  expect_false("log.hazard.ratio" %in% colnames(fit$trt.boot))
  expect_true(all(c("log.cumulative.hazard.ratio", "risk0", "risk1", "risk.difference") %in% colnames(fit$trt.boot)))
  expect_equal(fit$risk.difference$estimate, fit$risk1$estimate - fit$risk0$estimate)
  expect_equal(fit$trt.boot[, "risk.difference"], fit$trt.boot[, "risk1"] - fit$trt.boot[, "risk0"], ignore_attr = TRUE)
  expect_equal(fit$risk.difference$SE, sd(fit$trt.boot[, "risk.difference"]))
  expect_true("risk.difference" %in% names(coef(fit)))
  p <- plot(fit, metrics = c("risk0", "risk1", "risk.difference"))
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_equal(nrow(built$data[[3]]), 1L)
  expect_equal(built$data[[3]]$xintercept, 0)
  expect_equal(fit$treatment.labels, as.character(sort(unique(dat$trt))))
  expect_equal(survival_event_at_horizon(c(tau, tau, tau + 1), c(1, 0, 1), tau),
               c(1, 0, 0))
})
