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

test_that("catefit and catecv expose cumulative hazards alongside Cox HRs", {
  data(survivalExample, package = "precmed")
  common <- list(data = survivalExample,
                 cate.model = survival::Surv(y, d) ~ age + female + previous_cost + previous_number_relapses,
                 ps.model = trt ~ age + previous_treatment, tau0 = 522.66,
                 rmst.control = list(solver = "R"))
  options <- list(score.method = "contrastReg", initial.predictor.method = "logistic",
                  prop.cutoff = c(0.5, 1), contrast.control = list(solver = "R"), seed = 7)
  ate <- suppressWarnings(do.call(atefitsurv, c(common, list(n.boot = 2, seed = 1))))
  fit <- suppressWarnings(do.call(catefitsurv, c(common, options)))
  cv <- suppressWarnings(do.call(catecvsurv, c(common, options,
                     list(prop.multi = c(0, 0.5, 1), cv.n = 2, abc = FALSE))))
  expect_equal(unname(fit$chr.contrastReg["prop1"]),
               exp(ate$log.cumulative.hazard.ratio$estimate))
  expect_length(fit$hr.contrastReg, 2)
  expect_equal(unname(cv$chr.contrastReg$chr.est.valid.high.cv["prop1", ]),
               cv$overall.chr.valid)
  expect_equal(dim(cv$chr.contrastReg$chr.est.valid.group.cv), c(2L, 2L))
  expect_equal(dim(cv$hr.contrastReg$hr.est.valid.group.cv), c(2L, 2L))
})
