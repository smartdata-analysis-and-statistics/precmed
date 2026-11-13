# Local only: opt in using PRECMED_RUN_LOCAL_TESTS=true and NOT_CRAN=true.
testthat::test_that("C contrast regression agrees with the original R implementation", {
  testthat::skip_if(Sys.getenv("PRECMED_RUN_LOCAL_TESTS") != "true", "Requires local opt-in")
  testthat::skip_if(tolower(Sys.getenv("CI")) %in% c("true", "1") ||
                      nzchar(Sys.getenv("GITHUB_ACTIONS")) ||
                      nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")), "Local tests only")
  testthat::skip_on_cran()
  make_data <- function(seed, n = 800L) {
    set.seed(seed)
    x <- cbind(x1 = rnorm(n), x2 = rbinom(n, 1, .5), x3 = runif(n))
    ps <- plogis(.2 * x[, 1] - .3 * x[, 2])
    trt <- rbinom(n, 1, ps)
    f0 <- exp(.2 + .1 * x[, 1] + .2 * x[, 2])
    f1 <- f0 * exp(-.3 + .2 * x[, 1] - .1 * x[, 2])
    time <- log(runif(n, .5, 2))
    list(y = rpois(n, exp(time) * ifelse(trt == 1, f1, f0)),
         x.cate = x, time = time, trt = trt, ps = ps,
         f1.predictor = f1, f0.predictor = f0)
  }
  max.diff <- 0
  for (seed in c(1L, 42L, 999L)) {
    args <- make_data(seed)
    for (controls in list(list(), list(max.iterNR = 0),
                          list(tune = c(.7, 1.5), error.maxNR = 1e-5))) {
      reference <- do.call(original_contrast_count, c(args, controls))
      r <- do.call(precmed:::twoarmglmcount.dr, c(args, controls, list(solver = "R")))
      native <- do.call(precmed:::twoarmglmcount.dr, c(args, controls, list(solver = "C")))
      testthat::expect_identical(r, reference)
      testthat::expect_identical(
        do.call(precmed:::twoarmglmcount.dr, c(args, controls)), native)
      # LAPACK eigen/BLAS rounding can alter Nelder-Mead's path on fallback.
      testthat::expect_equal(native$coef, r$coef, tolerance = 1e-7)
      testthat::expect_equal(native$vcov, r$vcov, tolerance = 1e-7)
      testthat::expect_identical(native$converge, r$converge)
      max.diff <- max(max.diff, abs(native$coef - r$coef))
    }
  }
  # Explicit fallback controls and non-convergence must apply to both solvers.
  args <- make_data(42L)
  for (solver in c("R", "C")) {
    fit <- do.call(precmed:::twoarmglmcount.dr, c(args,
      list(solver = solver, max.iterNR = 0, optim.control = list(maxit = 1))))
    testthat::expect_false(fit$converge)
  }
  # No predictors: keep a one-column intercept-only design in the C backend.
  args$x.cate <- matrix(numeric(), length(args$y), 0)
  testthat::expect_equal(do.call(precmed:::twoarmglmcount.dr, c(args, list(solver = "C"))),
                         do.call(original_contrast_count, args), tolerance = 1e-7)
  # Singular designs should report errors in both implementations.
  args <- make_data(42L)
  args$x.cate <- cbind(args$x.cate, duplicate = args$x.cate[, 1])
  testthat::expect_error(do.call(original_contrast_count, args))
  testthat::expect_error(do.call(precmed:::twoarmglmcount.dr, c(args, list(solver = "C"))))

  # Test the exported fit/CV interfaces with real cross-fitting.
  args <- make_data(7L, 500L)
  data <- data.frame(y = args$y, trt = args$trt, args$x.cate, time = args$time)
  common <- list(response = "count", data = data,
                 cate.model = y ~ x1 + x2 + x3 + offset(time),
                 ps.model = trt ~ x1 + x2,
                 score.method = "contrastReg", initial.predictor.method = "poisson",
                 B = 1, Kfold = 3, seed = 42)
  r <- do.call(precmed::catefit, c(common, list(contrast.control = list(solver = "R"))))
  native <- do.call(precmed::catefit, c(common, list(contrast.control = list(solver = "C"))))
  testthat::expect_equal(native, r, tolerance = 1e-7)
  r <- do.call(precmed::catecv, c(common, list(cv.n = 2, contrast.control = list(solver = "R"))))
  native <- do.call(precmed::catecv, c(common, list(cv.n = 2, contrast.control = list(solver = "C"))))
  testthat::expect_equal(native, r, tolerance = 1e-7)

  testthat::expect_identical(
    precmed:::resolve_contrast_control(list(), .001, 150, c(.5, 2))$solver, "C")
  testthat::expect_equal(do.call(precmed::catefit, common),
                         do.call(precmed::catefit, c(common, list(contrast.control = list(solver = "C")))))

  # Named controls reject typos and override legacy NR arguments explicitly.
  testthat::expect_error(precmed:::resolve_contrast_control(
    list(solvr = "C"), .001, 150, c(.5, 2)), "named list")
  testthat::expect_error(precmed:::resolve_contrast_control(
    list(solver = "Fortran"), .001, 150, c(.5, 2)), "arg")
  testthat::expect_error(precmed:::resolve_contrast_control(
    list(optim.control = list(maxitt = 5)), .001, 150, c(.5, 2)), "optim.control")
  testthat::expect_error(precmed:::resolve_contrast_control(
    list(max.iterNR = NA_real_), .001, 150, c(.5, 2)), "max.iterNR")
  ctrl <- precmed:::resolve_contrast_control(
    list(solver = "C", max.iterNR = 0, optim.control = list(maxit = 1)),
    .001, 150, c(.5, 2))
  testthat::expect_identical(ctrl$max.iterNR, 0)
  testthat::expect_identical(ctrl$optim.control, list(maxit = 1))
  testthat::expect_error(precmed::catecv(response = "continuous", data = data,
    cate.model = y ~ x1, ps.model = trt ~ x1,
    contrast.control = list(solver = "C")), "count or survival")
  bad <- make_data(42L)
  bad$ps[1] <- 0
  testthat::expect_error(do.call(precmed:::twoarmglmcount.dr,
    c(bad, list(solver = "C"))), "0 < ps < 1")
  bad <- make_data(42L)
  bad$x.cate[1, 1] <- NA_real_
  testthat::expect_error(do.call(precmed:::twoarmglmcount.dr,
    c(bad, list(solver = "C"))), "finite numeric")

  # Benchmark complete solver, including covariance, without a timing assertion.
  args <- make_data(19L, 10000L)
  timings <- sapply(c("R", "C"), function(solver) replicate(3,
    system.time(do.call(precmed:::twoarmglmcount.dr, c(args, list(solver = solver))))[["elapsed"]]))
  cat(sprintf("\nMaximum coefficient difference: %.3g\n", max.diff))
  cat(sprintf("Median solver time: R %.3fs, C %.3fs (%.2fx)\n",
              median(timings[, "R"]), median(timings[, "C"]),
              median(timings[, "R"]) / median(timings[, "C"])))
})
