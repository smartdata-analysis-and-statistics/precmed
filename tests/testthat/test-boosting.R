boosting_example <- function() {
  set.seed(81)
  n <- 300L
  x <- cbind(x1 = rnorm(n), x2 = rnorm(n))
  trt <- rep(0:1, each = n / 2)
  time <- log(runif(n, 0.2, 2))
  y <- rpois(n, exp(0.2 + 0.3 * x[, 1] - 0.2 * trt + time))
  list(y = y, trt = trt, x = x, time = time,
       data = data.frame(y, x, time))
}

test_that("GBM3 adapter preserves historical training settings and offsets", {
  d <- boosting_example()$data
  set.seed(51)
  output <- capture.output(fit <- precmed:::fit_boosting(
    y ~ x1 + x2 + offset(time), d, "poisson", n.trees = 20,
    interaction.depth = 2, cv.folds = 3, n.cores = 1))
  expect_length(output, 0)
  expect_s3_class(fit, "GBMFit")
  expect_equal(fit$params$shrinkage, 0.1)
  expect_equal(fit$params$num_train, nrow(d))
  expect_equal(fit$params$num_features, 2)
  best <- precmed:::best_boosting_iteration(fit, "cv", plot.it = FALSE)
  expect_true(is.numeric(best) && length(best) == 1 && best >= 1 && best <= 20)
  link <- suppressWarnings(predict(fit, d, n.trees = best))
  response <- suppressWarnings(predict(fit, d, n.trees = best, type = "response"))
  expect_equal(response, exp(link))
  newdata <- d
  newdata$time <- newdata$time + log(2)
  expect_equal(suppressWarnings(predict(fit, newdata, n.trees = best)), link)
  expect_error(precmed:::fit_boosting(y ~ x1 + x2, d, "gaussian",
    n.cores = 1, par.details = gbm3::gbmParallel(1)), "only one")
})

test_that("count scores and boosting nuisance fits work with GBM3", {
  d <- boosting_example()
  set.seed(51)
  fit <- precmed:::intxcount(y = d$y, trt = d$trt, x.cate = d$x,
    x.ps = cbind(1, d$x), time = d$time,
    score.method = c("boosting", "twoReg"), initial.predictor.method = "boosting",
    B = 1, Kfold = 2, n.trees.boosting = 20, plot.gbmperf = FALSE, n.cores = 1)
  scores <- precmed:::scorecount(fit, d$x, d$time,
                                 score.method = c("boosting", "twoReg"))
  expect_true(all(is.finite(scores$score.boosting)))
  expect_true(all(is.finite(scores$score.twoReg)))
  expect_s3_class(fit$result.boosting$fit0.boosting, "GBMFit")
})

test_that("single-predictor count boosting works in fitting and validation", {
  set.seed(81)
  n <- 500L
  dat <- data.frame(x = rnorm(n), trt = rep(0:1, each = n / 2),
                    years = runif(n, 0.3, 2))
  dat$y <- rpois(n, exp(0.4 + 0.25 * dat$x - 0.2 * dat$trt) * dat$years)
  args <- list(data = dat, score.method = c("boosting", "twoReg", "contrastReg"),
               cate.model = y ~ x + offset(log(years)), ps.model = trt ~ x,
               initial.predictor.method = "boosting", B = 1, Kfold = 2,
               n.trees.boosting = 30, plot.gbmperf = FALSE, n.cores = 1,
               seed = 51)
  messages <- character()
  fit <- withCallingHandlers(do.call(catefitcount, args), warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  for (method in args$score.method) {
    expect_true(all(is.finite(fit[[paste0("score.", method)]])))
  }
  expect_false(any(grepl("OOB|out.of.bag", messages, ignore.case = TRUE)))
  messages <- character()
  cv <- withCallingHandlers(do.call(catecvcount, c(args, list(cv.n = 1,
    abc = FALSE, error.max = 0.5))), warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  for (method in args$score.method) {
    expect_true(all(is.finite(cv[[paste0("ate.", method)]]$ate.est.valid.high.cv)))
  }
  expect_false(any(grepl("OOB|out.of.bag", messages, ignore.case = TRUE)))
  args$score.method <- "boosting"
  args$seed <- 11
  messages <- character()
  tied <- withCallingHandlers(do.call(catecvcount, c(args, list(cv.n = 1,
    abc = TRUE, error.max = 0.5))), warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("score.boosting, validation iteration 1, low subgroup", messages, fixed = TRUE)))
  expect_true(anyNA(tied$ate.boosting$ate.est.valid.low.cv))
  expect_true(any(is.finite(tied$ate.boosting$ate.est.valid.high.cv)))
})

test_that("continuous boosting retains both selected tree counts", {
  d <- boosting_example()
  set.seed(51)
  fit <- precmed:::intxmean(y = log1p(d$y), trt = d$trt,
    x.cate = d$x, x.init = d$x, x.ps = cbind(1, d$x), score.method = "boosting",
    B = 1, Kfold = 2, n.trees.boosting = 20, plot.gbmperf = FALSE,
    shrinkage = 0.05, n.cores = 1)
  expect_true(is.finite(fit$result.boosting$best0.iter))
  expect_true(is.finite(fit$result.boosting$best1.iter))
  expect_equal(fit$result.boosting$fit0.boosting$params$shrinkage, 0.05)
  scores <- precmed:::scoremean(fit, d$x, score.method = "boosting")
  expect_true(all(is.finite(scores$score.boosting)))
})

test_that("survival boosting supports weights and generic prediction", {
  d <- boosting_example()
  set.seed(451)
  event <- rexp(length(d$y), rate = exp(0.2 * d$x[, 1] - 0.3 * d$trt))
  censor <- rexp(length(d$y), rate = 0.4)
  fit <- precmed:::intxsurv(y = pmin(event, censor), d = as.numeric(event <= censor),
    trt = d$trt, x.cate = d$x, x.ps = cbind(1, d$x), x.ipcw = d$x,
    tau0 = 2, score.method = "boosting", B = 1, Kfold = 2,
    n.trees.boosting = 20, plot.gbmperf = FALSE, n.cores = 1)
  expect_s3_class(fit$result.boosting$fit1.boosting, "GBMFit")
  score <- precmed:::scoresurv(fit, d$x, tau0 = 2, score.method = "boosting")
  expect_true(all(is.finite(score$score.boosting)))
})
