# Run explicitly from the repository root:
# Sys.setenv(PRECMED_RUN_LOCAL_TESTS = "true", NOT_CRAN = "true")
# pkgload::load_all(".")
# testthat::test_dir("tests/local")
# This directory is excluded from source packages and automatic test discovery.

testthat::test_that("GBM3 count cross-fitting is reproducible without redundant fits", {
  testthat::skip_if(Sys.getenv("PRECMED_RUN_LOCAL_TESTS") != "true",
                    "Requires an explicit local opt-in")
  testthat::skip_if(tolower(Sys.getenv("CI")) %in% c("true", "1") ||
                      nzchar(Sys.getenv("GITHUB_ACTIONS")) ||
                      nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")),
                    "Local regression test must not run in CI or R CMD check")
  testthat::skip_on_cran()

  set.seed(2026)
  n <- 600L
  x <- cbind(age = rnorm(n), male = rbinom(n, 1, 0.5))
  trt <- rep(0:1, each = n / 2)
  time <- log(runif(n, 0.5, 2))
  y <- rpois(n, exp(0.3 + 0.2 * x[, 1] - 0.1 * x[, 2] +
                     trt * (-0.3 + 0.1 * x[, 1]) + time))

  run_fit <- function(seed) {
    calls <- 0L
    checked <- precmed:::intxcount
    environment(checked) <- list2env(list(fit_boosting = function(...) {
      calls <<- calls + 1L
      precmed:::fit_boosting(...)
    }), parent = environment(checked))
    set.seed(seed)
    result <- checked(y = y, trt = trt, x.cate = x, x.ps = cbind(1, x), time = time,
      score.method = c("twoReg", "contrastReg"), initial.predictor.method = "boosting",
      B = 2, Kfold = 3, n.trees.boosting = 50, plot.gbmperf = FALSE, n.cores = 1)
    list(result = result, calls = calls)
  }
  for (seed in c(11L, 29L)) {
    first <- run_fit(seed)
    second <- run_fit(seed)
    testthat::expect_equal(first, second)
    testthat::expect_equal(first$calls, 2L * 2L * 3L)
    testthat::expect_true(all(is.finite(first$result$result.twoReg)))
    testthat::expect_true(all(is.finite(first$result$result.contrastReg$delta.contrastReg)))
  }
})
