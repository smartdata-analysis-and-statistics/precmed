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


test_that("horizon risks agree with empirical event risks and allow boundaries", {
  trt <- rep(0:1, each = 100)
  y <- c(rep(0.5, 20), rep(2, 80), rep(0.5, 40), rep(2, 60))
  s <- aipw_horizon_survival(y, rep(1, 200), trt, rep(.8, 200),
                             rep(.6, 200), rep(.5, 200), rep(1, 200), 1)
  expect_equal(unname(1 - s), c(.2, .4))
  expect_equal(unname(s[1] - s[2]), .2)
  s <- aipw_horizon_survival(rep(2, 20), rep(1, 20), rep(0:1, each = 10),
                             rep(1, 20), rep(1, 20), rep(.5, 20), rep(1, 20), 1)
  expect_equal(unname(s[1] - s[2]), 0)
})
