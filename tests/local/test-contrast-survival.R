# Opt-in local regression/benchmark tests; excluded from source packages.
testthat::test_that("native survival contrast regression matches the original R solver", {
  testthat::skip_if(Sys.getenv("PRECMED_RUN_LOCAL_TESTS") != "true", "Requires local opt-in")
  testthat::skip_if(tolower(Sys.getenv("CI")) %in% c("true", "1") ||
                      nzchar(Sys.getenv("GITHUB_ACTIONS")) ||
                      nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")), "Local tests only")
  testthat::skip_on_cran()
  make_args <- function(seed, n = 1000L, tau0 = 2) {
    set.seed(seed)
    x <- cbind(x1 = rnorm(n), x2 = runif(n), x3 = rbinom(n, 1, .5))
    ps <- plogis(.2 * x[, 1] - .3 * x[, 3])
    trt <- rbinom(n, 1, ps)
    mu0 <- plogis(-1 + .1 * x[, 1] + .2 * x[, 2])
    mu1 <- pmin(.95, mu0 * exp(-.2 + .1 * x[, 1] - .1 * x[, 3]))
    mu <- ifelse(trt == 1, mu1, mu0)
    weights <- rexp(n)
    weights[seq(1, n, by = 11)] <- 0
    list(ynew = tau0 * (1 - rbeta(n, 20 * mu, 20 * (1 - mu))),
         dnew = as.numeric(weights > 0), trt = trt, x.cate = x, tau0 = tau0,
         weightsurv = weights, ps = ps,
         f1.predictor = tau0 * (1 - mu1), f0.predictor = tau0 * (1 - mu0))
  }
  differences <- numeric()
  compare <- function(args, controls = list()) {
    before <- do.call(original_contrast_survival, c(args, controls))
    r <- do.call(precmed:::twoarmsurv.dr, c(args, controls, list(solver = "R")))
    native <- do.call(precmed:::twoarmsurv.dr, c(args, controls, list(solver = "C")))
    testthat::expect_identical(r, before)
    testthat::expect_equal(native$coef, r$coef, tolerance = 1e-7)
    testthat::expect_identical(native$converge, r$converge)
    testthat::expect_identical(do.call(precmed:::twoarmsurv.dr, c(args, controls)), native)
    differences <<- c(differences, max(abs(native$coef - r$coef)))
    native
  }
  for (seed in c(1L, 42L, 999L)) {
    args <- make_args(seed)
    testthat::expect_true(compare(args)$converge)
    compare(args, list(max.iterNR = 0))
    compare(args, list(tune = c(.7, 1.5), error.maxNR = 1e-5))
    # A small tau forces the historical survival-specific eta cap and fallback.
    compare(make_args(seed, tau0 = .0001))
  }
  args <- make_args(42L)
  args$x.cate <- matrix(numeric(), length(args$ynew), 0)
  compare(args)
  args <- make_args(42L)
  args$weightsurv <- rep(1, length(args$ynew))
  compare(args)
  # Normalize arbitrary weight scale identically in both implementations.
  args$weightsurv <- args$weightsurv * 17
  compare(args)

  # Aggressive steps overflow eta; both implementations must use fallback.
  compare(make_args(42L), list(tune = c(1e6, 2)))
  # Survival must continue with first-iteration regularization even when the
  # initial slope is negative (unlike the count solver's positivity guard).
  args <- make_args(42L)
  args$ynew[] <- args$tau0
  args$f1.predictor[] <- .9 * args$tau0
  args$f0.predictor[] <- .1 * args$tau0
  args$weightsurv <- ifelse(args$trt == 0, 2, .2)
  w <- args$weightsurv / mean(args$weightsurv)
  ps <- args$ps
  slopewt <- (.9 * (args$trt / ps - 1) / 2 +
              .1 * ((1 - args$trt) / (1 - ps) - 1) / 2) * ps * (1 - ps)
  x <- cbind(1, args$x.cate)
  testthat::expect_lt(min(eigen(crossprod(x * w * slopewt, x),
                               only.values = TRUE)$values), 0)
  compare(args)

  for (solver in c("R", "C")) {
    args <- make_args(42L)
    fit <- do.call(precmed:::twoarmsurv.dr, c(args,
      list(solver = solver, max.iterNR = 0, optim.control = list(maxit = 1))))
    testthat::expect_false(fit$converge)
  }
  for (bad in c(0, NA_real_, Inf)) {
    args <- make_args(42L)
    args$tau0 <- bad
    testthat::expect_error(do.call(precmed:::twoarmsurv.dr, args), "tau0")
  }
  args <- make_args(42L)
  args$weightsurv[] <- 0
  testthat::expect_error(do.call(precmed:::twoarmsurv.dr, args), "weights")
  args <- make_args(42L)
  args$ps[1] <- 0
  testthat::expect_error(do.call(precmed:::twoarmsurv.dr, args), "0 < ps < 1")

  # Real public fit/CV pipelines, including cross-fitting and RMST evaluation.
  set.seed(19)
  n <- 500L
  x1 <- rnorm(n); x2 <- runif(n); trt <- rbinom(n, 1, .5)
  event <- rexp(n, exp(.1 * x1 - .2 * trt))
  censor <- rexp(n, .6)
  data <- data.frame(y = pmin(event, censor), d = as.numeric(event <= censor), trt, x1, x2)
  common <- list(response = "survival", data = data,
                 cate.model = survival::Surv(y, d) ~ x1 + x2,
                 ps.model = trt ~ x1 + x2, ipcw.model = ~ x1 + x2,
                 tau0 = 1.5, score.method = "contrastReg", initial.predictor.method = "logistic",
                 B = 1, Kfold = 3, seed = 42, prop.cutoff = c(.5, 1))
  for (fun in list(precmed::catefit, precmed::catecv)) {
    extra <- if (identical(fun, precmed::catecv)) list(cv.n = 2, max.iter = 100) else list()
    r <- do.call(fun, c(common, extra, list(contrast.control = list(solver = "R"))))
    native <- do.call(fun, c(common, extra, list(contrast.control = list(solver = "C"))))
    testthat::expect_equal(native, r, tolerance = 1e-7)
    testthat::expect_identical(do.call(fun, c(common, extra)), native)
  }
  args <- make_args(19L, n = 10000L)
  timings <- sapply(c("R", "C"), function(solver) replicate(3,
    system.time(do.call(precmed:::twoarmsurv.dr, c(args, list(solver = solver))))[["elapsed"]]))
  cat(sprintf("\nSurvival solver maximum coefficient difference: %.3g\n", max(differences)))
  cat(sprintf("Median solver: R %.3fs; C %.3fs (%.2fx)\n",
              median(timings[, "R"]), median(timings[, "C"]),
              median(timings[, "R"]) / median(timings[, "C"])))
})
