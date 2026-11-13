# Local-only comparisons: excluded from source packages and normal test discovery.
testthat::test_that("continuous balanced splitting preserves estimates and selected splits", {
  testthat::skip_if(Sys.getenv("PRECMED_RUN_LOCAL_TESTS") != "true", "Requires local opt-in")
  testthat::skip_if(tolower(Sys.getenv("CI")) %in% c("true", "1") ||
                      nzchar(Sys.getenv("GITHUB_ACTIONS")) ||
                      nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")), "Local tests only")
  testthat::skip_on_cran()
  set.seed(48)
  n <- 600L
  x <- cbind(x1 = rnorm(n), x2 = rbinom(n, 1, .4), x3 = runif(n))
  trt <- rbinom(n, 1, plogis(.2 * x[, 1]))
  y <- .3 + .2*x[, 1] - .1*x[, 2] + trt*(.5 + .1*x[, 1]) + rnorm(n)
  families <- list(outcome = stats::gaussian(), propensity = stats::binomial())
  compare_estimates <- function(ids, cate = x, ps = cbind(1, x), outcome = y,
                                minPS = .01, maxPS = .99) {
    before <- precmed:::drmean(outcome[ids], trt[ids], cate[ids, , drop = FALSE],
                               ps[ids, , drop = FALSE], minPS = minPS, maxPS = maxPS)
    after <- precmed:::drmean.matrix(outcome[ids], trt[ids],
                                     cbind(1, cate)[ids, , drop = FALSE],
                                     ps[ids, , drop = FALSE], families, minPS, maxPS)
    testthat::expect_equal(after, before, tolerance = 1e-12)
  }
  compare_estimates(seq_len(n))
  compare_estimates(seq_len(n), ps = matrix(1, n, 1))
  compare_estimates(seq_len(n), cate = cbind(x, duplicate = x[, 1]),
                    ps = cbind(1, x, duplicate = x[, 1]))
  compare_estimates(c(which(trt == 1)[1:8], which(trt == 0)[1:50]))
  compare_estimates(seq_len(n), minPS = .4, maxPS = .6)
  compare_estimates(seq_len(n), outcome = rep(0, n))

  args <- list(y = y, trt = trt, x.cate = x, x.ps = cbind(1, x))
  run <- function(fun, seed, cutoff, max.iter) {
    set.seed(seed)
    warnings <- character()
    elapsed <- system.time(value <- withCallingHandlers(
      do.call(fun, c(args, list(error.max = cutoff, max.iter = max.iter))),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ))[["elapsed"]]
    list(value = value, rng = .Random.seed, warnings = warnings, elapsed = elapsed)
  }
  for (seed in c(1L, 42L, 999L)) {
    for (cutoff in c(.1, 0, 10)) {
      before <- run(original_balancemean_split, seed, cutoff, 75)
      after <- run(precmed:::balancemean.split, seed, cutoff, 75)
      testthat::expect_identical(after$value, before$value)
      testthat::expect_identical(after$rng, before$rng)
      testthat::expect_identical(after$warnings, before$warnings)
      if (cutoff == 10) testthat::expect_length(after$warnings, 0L)
    }
  }
  # The single-candidate/max-iteration path must preserve the saved best IDs.
  before <- run(original_balancemean_split, 42L, 0, 0)
  after <- run(precmed:::balancemean.split, 42L, 0, 0)
  testthat::expect_identical(after$value, before$value)
  testthat::expect_identical(after$warnings, before$warnings)
  testthat::expect_identical(after$rng, before$rng)

  # Preserve the original missing-data error and consumed RNG on fallback.
  args$x.cate[1, 1] <- NA_real_
  set.seed(42L)
  testthat::expect_error(do.call(original_balancemean_split,
    c(args, list(error.max = 0, max.iter = 5))), "No balanced data split found")
  rng <- .Random.seed
  set.seed(42L)
  testthat::expect_error(do.call(precmed:::balancemean.split,
    c(args, list(error.max = 0, max.iter = 5))), "No balanced data split found")
  testthat::expect_identical(.Random.seed, rng)
  args$x.cate <- x

  # Zero-column outcome designs were unsupported by the formula path; do not
  # silently change that behavior as part of a performance-only optimization.
  args$x.cate <- matrix(numeric(), n, 0)
  set.seed(42L)
  testthat::expect_error(do.call(original_balancemean_split,
    c(args, list(error.max = 10, max.iter = 1))), "subscript out of bounds")
  rng <- .Random.seed
  set.seed(42L)
  testthat::expect_error(do.call(precmed:::balancemean.split,
    c(args, list(error.max = 10, max.iter = 1))), "subscript out of bounds")
  testthat::expect_identical(.Random.seed, rng)
  args$x.cate <- x

  # This optimization must not enable the intentionally disabled public APIs.
  common <- list(response = "continuous", data = data.frame(y, trt, x),
                  cate.model = y ~ x1 + x2, ps.model = trt ~ x1 + x2)
  testthat::expect_error(do.call(precmed::atefit, common), "not implemented yet")
  testthat::expect_error(do.call(precmed::catefit, c(common, list(score.method = "gaussian"))), "not implemented yet")
  testthat::expect_error(do.call(precmed::catecv, c(common, list(score.method = "gaussian"))), "not implemented yet")

  # Informational timings, no wall-clock assertion; same candidate budget/seeds.
  timings <- sapply(c(11L, 29L, 53L), function(seed) c(
    original = run(original_balancemean_split, seed, 0, 100)$elapsed,
    revised = run(precmed:::balancemean.split, seed, 0, 100)$elapsed))
  cat(sprintf("\nContinuous split median: original %.3fs; revised %.3fs (%.2fx)\n",
              median(timings[1, ]), median(timings[2, ]),
              median(timings[1, ]) / median(timings[2, ])))
})
