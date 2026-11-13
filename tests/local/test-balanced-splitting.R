# Opt-in only; same local runner and CI/CRAN guards as test-count-boosting.R.
testthat::test_that("matrix GLMs preserve count estimates and balanced splits", {
  testthat::skip_if(Sys.getenv("PRECMED_RUN_LOCAL_TESTS") != "true",
                    "Requires an explicit local opt-in")
  testthat::skip_if(tolower(Sys.getenv("CI")) %in% c("true", "1") ||
                      nzchar(Sys.getenv("GITHUB_ACTIONS")) ||
                      nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")),
                    "Local regression test must not run in CI or R CMD check")
  testthat::skip_on_cran()

  set.seed(48)
  n <- 600L
  x <- cbind(x1 = rnorm(n), x2 = rbinom(n, 1, 0.4))
  trt <- rbinom(n, 1, plogis(0.2 * x[, 1]))
  time <- log(runif(n, 0.5, 2))
  y <- rpois(n, exp(0.3 + 0.2 * x[, 1] - 0.3 * trt + time))
  families <- list(outcome = poisson(), propensity = binomial())

  compare_estimates <- function(ids, cate = x, ps = cbind(1, x),
                                outcome = y, minPS = 0.01, maxPS = 0.99) {
    before <- precmed:::drcount(outcome[ids], trt[ids], cate[ids, , drop = FALSE],
                                ps[ids, , drop = FALSE], time[ids],
                                minPS = minPS, maxPS = maxPS)
    after <- precmed:::drcount.matrix(outcome[ids], trt[ids],
                                      cbind(1, cate)[ids, , drop = FALSE],
                                      ps[ids, , drop = FALSE], time[ids], families,
                                      minPS = minPS, maxPS = maxPS)
    testthat::expect_equal(after, before, tolerance = 1e-12)
  }
  compare_estimates(seq_len(n))
  compare_estimates(seq_len(n), ps = matrix(1, n, 1))
  compare_estimates(seq_len(n), cate = cbind(x, duplicate = x[, 1]),
                    ps = cbind(1, x, duplicate = x[, 1]))
  compare_estimates(c(which(trt == 1)[1:8], which(trt == 0)[1:50]))
  compare_estimates(seq_len(n), minPS = 0.4, maxPS = 0.6)
  compare_estimates(seq_len(n), outcome = rep(0, n))

  args <- list(y = y, trt = trt, x.cate = x, x.ps = cbind(1, x), time = time)
  run_split <- function(fun, seed, error.max, max.iter) {
    set.seed(seed)
    warnings <- character()
    elapsed <- system.time(value <- withCallingHandlers(
      do.call(fun, c(args, list(error.max = error.max, max.iter = max.iter))),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ))[["elapsed"]]
    list(value = value, rng = .Random.seed, warnings = warnings, elapsed = elapsed)
  }
  for (seed in c(1L, 42L, 999L)) {
    for (cutoff in c(0.1, 0)) {
      before <- run_split(original_balance_split, seed, cutoff, 75)
      after <- run_split(precmed:::balance.split, seed, cutoff, 75)
      testthat::expect_identical(after$value, before$value)
      testthat::expect_identical(after$rng, before$rng)
      testthat::expect_identical(after$warnings, before$warnings)
    }
  }

  # Missing values retain the original formula-based path (including na.action).
  args$x.cate[1, 1] <- NA_real_
  set.seed(42L)
  testthat::expect_error(do.call(original_balance_split,
    c(args, list(error.max = 0, max.iter = 5))), "No balanced data split found")
  before.rng <- .Random.seed
  set.seed(42L)
  testthat::expect_error(do.call(precmed:::balance.split,
    c(args, list(error.max = 0, max.iter = 5))), "No balanced data split found")
  testthat::expect_identical(.Random.seed, before.rng)
  args$x.cate <- x

  # Compare the complete five-split workflow, isolating only balance.split.
  cv <- precmed::catecvcount
  original.env <- new.env(parent = environment(cv))
  original.env$balance.split <- original_balance_split
  environment(cv) <- original.env
  data <- data.frame(y = y, trt = trt, x1 = x[, 1], x2 = x[, 2], time = time)
  run_cv <- function(fun) {
    invisible(capture.output(value <- fun(
      data = data, cate.model = y ~ x1 + x2 + offset(time), ps.model = trt ~ x1 + x2,
      score.method = c("poisson", "contrastReg"), initial.predictor.method = "poisson",
      B = 1, Kfold = 3, cv.n = 5, seed = 999, verbose = 0, n.cores = 1
    )))
    value
  }
  before <- run_cv(cv)
  after <- run_cv(precmed::catecvcount)
  testthat::expect_equal(after, before, tolerance = 1e-12)

  # Timing is informational; do not assert a wall-clock threshold.
  timings <- sapply(c(11L, 29L, 53L), function(seed) {
    c(original = run_split(original_balance_split, seed, 0, 100)$elapsed,
      revised = run_split(precmed:::balance.split, seed, 0, 100)$elapsed)
  })
  message(sprintf("Balanced splitting median: original %.3fs; revised %.3fs (%.2fx)",
                  median(timings[1, ]), median(timings[2, ]),
                  median(timings[1, ]) / median(timings[2, ])))
  saveRDS(timings, file.path(tempdir(), "precmed-balanced-splitting-timings.rds"))
})
