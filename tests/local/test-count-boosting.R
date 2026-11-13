# Run explicitly from the repository root:
# Sys.setenv(PRECMED_RUN_LOCAL_TESTS = "true", NOT_CRAN = "true")
# pkgload::load_all(".")
# testthat::test_dir("tests/local")
# This directory is excluded from source packages and automatic test discovery.

testthat::test_that("count cross-fitting matches the original without discarded GBM fits", {
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

  run_fit <- function(fun, seed, initial.method) {
    calls <- 0L
    # Control GBM randomness without changing the outer fold-generation RNG.
    # Thus the original's discarded fit and retained fit are identical, and
    # subsequent folds/arms are comparable despite the removed random draws.
    controlled_gbm <- function(...) {
      calls <<- calls + 1L
      rng <- .Random.seed
      on.exit(assign(".Random.seed", rng, envir = .GlobalEnv))
      set.seed(415L)
      gbm::gbm(...)
    }
    environment(fun) <- list2env(list(gbm = controlled_gbm),
                                 parent = asNamespace("precmed"))
    set.seed(seed)
    invisible(capture.output(timing <- system.time(result <- fun(
      y = y, trt = trt, x.cate = x, x.ps = cbind(1, x), time = time,
      score.method = c("twoReg", "contrastReg"),
      initial.predictor.method = initial.method,
      B = 2, Kfold = 3, n.trees.boosting = 50,
      plot.gbmperf = FALSE, n.cores = 1
    ))))
    list(result = result, calls = calls, elapsed = unname(timing["elapsed"]))
  }

  for (seed in c(11L, 29L)) {
    before <- run_fit(original_intxcount, seed, "boosting")
    after <- run_fit(precmed:::intxcount, seed, "boosting")
    testthat::expect_equal(after$result, before$result, tolerance = 1e-10)
    testthat::expect_true(all(is.finite(after$result$result.twoReg)))
    testthat::expect_true(all(is.finite(
      after$result$result.contrastReg$delta.contrastReg)))
    testthat::expect_true(after$result$result.contrastReg$converge.contrastReg)
    testthat::expect_equal(before$calls, 3L * 2L * 3L)
    testthat::expect_equal(after$calls, 2L * 2L * 3L)
    message(sprintf("Seed %d: original %.2fs (%d GBM calls); revised %.2fs (%d calls)",
                    seed, before$elapsed, before$calls, after$elapsed, after$calls))
  }

  # The unrelated Poisson-initialization branch remains exactly equivalent.
  before <- run_fit(original_intxcount, 11L, "poisson")
  after <- run_fit(precmed:::intxcount, 11L, "poisson")
  testthat::expect_equal(after$result, before$result, tolerance = 1e-10)
  testthat::expect_equal(after$calls, 0L)
})
