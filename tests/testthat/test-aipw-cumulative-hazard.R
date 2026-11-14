test_that("AIPW cumulative hazard agrees with empirical survival in a balanced uncensored sample", {
  trt <- rep(0:1, each = 100)
  y <- c(rep(0.5, 20), rep(2, 80), rep(0.5, 40), rep(2, 60))
  estimate <- aipw_log_cumulative_hr(y, rep(1, 200), trt,
                                     rep(0.8, 200), rep(0.6, 200),
                                     rep(0.5, 200), rep(1, 200), 1)
  expect_equal(estimate, log(-log(0.6) / -log(0.8)), tolerance = 1e-8)
})

test_that("AIPW cumulative hazard is undefined when neither group has events", {
  expect_true(is.na(aipw_log_cumulative_hr(rep(2, 20), rep(1, 20), rep(0:1, each = 10),
                                           rep(1, 20), rep(1, 20),
                                           rep(0.5, 20), rep(1, 20), 1)))
})
