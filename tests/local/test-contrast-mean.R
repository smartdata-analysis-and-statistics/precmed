# Local opt-in comparisons/benchmarks; excluded from CI and CRAN source packages.
testthat::test_that("matrix continuous contrast regression matches the original", {
  testthat::skip_if(Sys.getenv("PRECMED_RUN_LOCAL_TESTS") != "true", "Requires local opt-in")
  testthat::skip_if(tolower(Sys.getenv("CI")) %in% c("true", "1") ||
                      nzchar(Sys.getenv("GITHUB_ACTIONS")) ||
                      nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")), "Local tests only")
  testthat::skip_on_cran()
  make_args <- function(seed, n = 500L, p = 5L) {
    set.seed(seed)
    x <- matrix(rnorm(n*p), n, p, dimnames = list(paste0("row", seq_len(n)), paste0("x", seq_len(p))))
    ps <- plogis(.2*x[, 1])
    trt <- rbinom(n, 1, ps)
    f0 <- .5 + .3*x[, 1]
    f1 <- f0 + .5 + .1*x[, 2]
    list(y = ifelse(trt == 1, f1, f0) + rnorm(n), x.cate = x, trt = trt,
         ps = ps, f1.predictor = f1, f0.predictor = f0)
  }
  compare <- function(args) {
    before <- do.call(original_contrast_mean, args)
    after <- do.call(precmed:::twoarmglmmean.dr, args)
    testthat::expect_identical(after$coef, before$coef)
    testthat::expect_identical(after$vcov, before$vcov)
  }
  for (seed in c(1L, 42L, 999L)) {
    args <- make_args(seed)
    compare(args)
    args$x.cate <- unname(args$x.cate)
    compare(args)
    args$x.cate <- cbind(args$x.cate, duplicate = args$x.cate[, 1])
    compare(args)
  }
  args <- make_args(42L)
  args$x.cate <- args$x.cate[, 1, drop = FALSE]
  compare(args)
  args$x.cate <- matrix(numeric(), length(args$y), 0)
  compare(args)
  args <- make_args(42L)
  args$ps <- ifelse(args$ps < .5, .01, .99)
  compare(args)
  args$ps[] <- .5
  compare(args)
  args <- make_args(29L, n = 30L, p = 40L)
  compare(args)

  # Missing responses/predictions preserve formula-based omission and warnings.
  run <- function(fun, args) {
    warnings <- character()
    result <- withCallingHandlers(tryCatch(do.call(fun, args),
      error = function(e) list(error = conditionMessage(e))),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    list(result = result, warnings = warnings)
  }
  for (field in c("y", "f1.predictor", "x.cate")) {
    args <- make_args(42L)
    args[[field]][1] <- NA_real_
    testthat::expect_identical(run(precmed:::twoarmglmmean.dr, args),
                               run(original_contrast_mean, args))
  }

  # Compare the real cross-fitting caller without enabling disabled public APIs.
  set.seed(42)
  n <- 600L
  x <- cbind(x1 = rnorm(n), x2 = runif(n))
  trt <- rbinom(n, 1, .5)
  y <- .5*x[,1] + trt*(.5 + .1*x[,2]) + rnorm(n)
  ref <- precmed:::intxmean
  env <- new.env(parent = environment(ref))
  env$twoarmglmmean.dr <- original_contrast_mean
  environment(ref) <- env
  args <- list(y = y, trt = trt, x.cate = x, x.init = x, x.ps = cbind(1, x),
                score.method = "contrastReg", initial.predictor.method = "gaussian",
                B = 2, Kfold = 3, plot.gbmperf = FALSE)
  set.seed(19); before <- do.call(ref, args); rng <- .Random.seed
  set.seed(19); after <- do.call(precmed:::intxmean, args)
  testthat::expect_identical(after, before)
  testthat::expect_identical(.Random.seed, rng)

  # Informational whole-function timing, including sandwich covariance.
  for (p in c(5L, 40L)) {
    args <- make_args(19L, n = 10000L, p = p)
    # Warm up before timing to reduce compilation/cache effects.
    invisible(do.call(original_contrast_mean, args))
    invisible(do.call(precmed:::twoarmglmmean.dr, args))
    batch <- function(fun) system.time(
      for (i in seq_len(20L)) invisible(do.call(fun, args))
    )[["elapsed"]] / 20
    timing <- replicate(3, c(original = batch(original_contrast_mean),
                              revised = batch(precmed:::twoarmglmmean.dr)))
    cat(sprintf("\nContinuous contrast p=%d: original %.4fs; revised %.4fs (%.2fx)\n",
                p, median(timing[1, ]), median(timing[2, ]),
                median(timing[1, ]) / median(timing[2, ])))
  }
})
