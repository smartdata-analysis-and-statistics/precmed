# Expensive CV comparisons are opt-in and excluded from CI/CRAN source packages.
testthat::test_that("PS LASSO reuses the CV fit without changing predictions or RNG", {
  testthat::skip_if(Sys.getenv("PRECMED_RUN_LOCAL_TESTS") != "true", "Requires local opt-in")
  testthat::skip_if(tolower(Sys.getenv("CI")) %in% c("true", "1") ||
                      nzchar(Sys.getenv("GITHUB_ACTIONS")) ||
                      nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")), "Local tests only")
  testthat::skip_on_cran()
  for (seed in c(1L, 42L, 999L)) {
    set.seed(seed)
    n <- 600L
    x <- matrix(rnorm(n*6), n, 6)
    trt <- rbinom(n, 1, plogis(.8*x[,1]*x[,2] + .2*x[,3]))
    testthat::expect_gt(ncol(precmed:::ps_interactions(x, trt)$pairs), 0L)
    for (solver in c("R", "C")) {
      for (ids in list(NULL, 1:17, 1L)) {
        args <- list(trt = trt, x.ps = cbind(1, x),
                      xnew = if (!is.null(ids)) cbind(1, x[ids, , drop = FALSE]),
                      solver = solver, minPS = .1, maxPS = .9)
        set.seed(19); before <- do.call(original_ps_lasso_fit, args); rng <- .Random.seed
        set.seed(19); after <- do.call(precmed:::glm.ps, args)
        testthat::expect_identical(after, before)
        testthat::expect_identical(.Random.seed, rng)
      }
    }
  }
  # Count direct full-data fits separately from fitting inside cv.glmnet's own
  # namespace. Its fold fits and its single full-data fit must remain intact.
  run <- function(fun) {
    direct <- cv <- 0L
    env <- new.env(parent = asNamespace("precmed"))
    env$glmnet <- function(...) {
      direct <<- direct + 1L
      glmnet::glmnet(...)
    }
    env$cv.glmnet <- function(...) {
      cv <<- cv + 1L
      glmnet::cv.glmnet(...)
    }
    environment(fun) <- env
    set.seed(19)
    value <- fun(trt, cbind(1, x))
    list(value = value, direct = direct, cv = cv)
  }
  before <- run(original_ps_lasso_fit)
  after <- run(precmed:::glm.ps)
  testthat::expect_identical(after$value, before$value)
  testthat::expect_equal(before$direct, 1L)
  testthat::expect_equal(after$direct, 0L)
  testthat::expect_equal(before$cv, 1L)
  testthat::expect_equal(after$cv, 1L)

  # No selected interactions still use the unchanged ordinary-logistic branch.
  trt <- rep(0:1, each = 200)
  x <- matrix(rep(seq(-2, 2, length.out = 200), 2), ncol = 1)
  testthat::expect_identical(precmed:::glm.ps(trt, cbind(1, x)),
                             original_ps_lasso_fit(trt, cbind(1, x)))
})
