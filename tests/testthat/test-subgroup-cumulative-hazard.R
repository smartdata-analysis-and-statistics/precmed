test_that("subgroup cumulative hazard uses the same estimator as the overall fit", {
  set.seed(18)
  n <- 800
  x <- matrix(rnorm(n), ncol = 1)
  trt <- rbinom(n, 1, plogis(0.2 * x[, 1]))
  event <- rexp(n, exp(-2 + 0.2 * x[, 1] + 0.3 * trt))
  censor <- rexp(n, 0.03)
  y <- pmin(event, censor)
  d <- as.integer(event <= censor)
  x.ps <- cbind(1, x)
  x.ipcw <- cbind(x, trt)
  args <- list(y = y, d = d, trt = trt, x.cate = x, x.ps = x.ps,
               x.ipcw = x.ipcw, yf = NULL, tau0 = 3,
               rmst.control = list(solver = "R"))
  overall <- do.call(drsurv, c(args, list(cumulative.hr = TRUE)))
  sub <- do.call(estsurv.bilevel.subgroups,
                 c(args, list(score = x[, 1], higher.y = TRUE,
                              prop = 1, onlyhigh = TRUE)))
  expect_equal(sub$chr.high, exp(overall$log.cumulative.hazard.ratio))
  expect_equal(sub$hr.high, exp(overall$log.hazard.ratio))
  expect_true(is.finite(sub$chr.high))
})

test_that("subgroup fitting failures preserve all result fields", {
  out <- survCatch(stop("deliberate failure"))
  expect_true(is.na(out$log.cumulative.hazard.ratio))
  expect_true(is.na(out$log.hazard.ratio))
  expect_match(out$errors, "deliberate failure")
})
