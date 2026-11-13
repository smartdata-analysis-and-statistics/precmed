# Real boosting/CV/bootstrap comparisons are local only and excluded from builds.
testthat::test_that("survival workflows preserve intercept-only PS matrix dimensions", {
  testthat::skip_if(Sys.getenv("PRECMED_RUN_LOCAL_TESTS") != "true", "Requires local opt-in")
  testthat::skip_if(tolower(Sys.getenv("CI")) %in% c("true", "1") ||
                      nzchar(Sys.getenv("GITHUB_ACTIONS")) ||
                      nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")), "Local tests only")
  testthat::skip_on_cran()
  # Recreate the four pre-fix slices in the two affected branches. Random-forest
  # cross-fitting already had drop=FALSE and remains unchanged in this reference.
  original <- precmed:::intxsurv
  code <- deparse(body(original), width.cutoff = 500L)
  for (pattern in c("x_ps_train <- x.ps[index != k, , drop = FALSE]",
                    "x_ps_valid <- x.ps[index == k, , drop = FALSE]")) {
    lines <- grep(pattern, code, fixed = TRUE)
    stopifnot(length(lines) == 3L)
    code[lines[-1L]] <- sub(", , drop = FALSE]", ", ]", code[lines[-1L]], fixed = TRUE)
  }
  body(original) <- parse(text = paste(code, collapse = "\n"))[[1L]]
  set.seed(19)
  n <- 600L
  x <- cbind(x1 = rnorm(n), x2 = runif(n))
  trt <- rbinom(n, 1, .5)
  event <- rexp(n, exp(.1*x[,1] - .2*trt)); censor <- rexp(n, .6)
  y <- pmin(event, censor); d <- as.numeric(event <= censor)
  args <- list(y = y, d = d, trt = trt, x.cate = x,
                x.ps = matrix(1, n, 1), x.ipcw = x, tau0 = 1.5,
                score.method = "contrastReg", B = 1, Kfold = 3,
                n.trees.boosting = 30, plot.gbmperf = FALSE, n.cores = 1)
  for (method in c("logistic", "boosting")) {
    args$initial.predictor.method <- method
    set.seed(42)
    testthat::expect_error(do.call(original, args), "incorrect number of dimensions")
    calls <- 0L
    checked <- precmed:::intxsurv
    env <- new.env(parent = environment(checked))
    env$glm.simplereg.ps <- function(trt, x.ps, xnew = NULL, minPS = .01, maxPS = .99) {
      calls <<- calls + 1L
      testthat::expect_true(is.matrix(x.ps))
      testthat::expect_true(is.matrix(xnew))
      testthat::expect_equal(ncol(x.ps), 1L)
      testthat::expect_equal(ncol(xnew), 1L)
      prob <- precmed:::glm.simplereg.ps(trt, x.ps, xnew, minPS, maxPS)
      testthat::expect_equal(prob, rep(pmax(minPS, pmin(maxPS, mean(trt))), nrow(xnew)),
                             tolerance = 1e-12)
      prob
    }
    environment(checked) <- env
    set.seed(42); fit <- do.call(checked, args)
    testthat::expect_equal(calls, 3L)
    testthat::expect_length(fit$result.contrastReg$delta.contrastReg, 3L)
    testthat::expect_true(all(is.finite(fit$result.contrastReg$delta.contrastReg)))
    testthat::expect_true(fit$result.contrastReg$converge.contrastReg)

    # Models with multiple PS columns must retain their original output and RNG.
    args$x.ps <- cbind(1, x)
    set.seed(42); before <- do.call(original, args); rng <- .Random.seed
    set.seed(42); after <- do.call(precmed:::intxsurv, args)
    testthat::expect_identical(after, before)
    testthat::expect_identical(.Random.seed, rng)
    args$x.ps <- matrix(1, n, 1)
  }

  # Exercise subgroup evaluation and bootstrap paths beyond the repaired folds.
  data <- data.frame(y = y, d = d, trt = trt, x)
  common <- list(response = "survival", data = data,
                  cate.model = survival::Surv(y, d) ~ x1 + x2,
                  ps.model = trt ~ 1, ipcw.model = ~ x1 + x2,
                  tau0 = 1.5, seed = 42)
  score <- list(score.method = "contrastReg", initial.predictor.method = "logistic",
                 B = 1, Kfold = 3, prop.cutoff = c(.5, 1))
  fit <- do.call(precmed::catefit, c(common, score))
  testthat::expect_s3_class(fit, "catefit")
  testthat::expect_length(fit$ate.contrastReg, 2L)
  testthat::expect_length(fit$hr.contrastReg, 2L)
  testthat::expect_true(all(is.finite(fit$ate.contrastReg)))
  testthat::expect_true(all(is.finite(fit$hr.contrastReg)))
  testthat::expect_null(fit[["errors/warnings"]]$contrastReg$est.high$errors)
  cv <- do.call(precmed::catecv, c(common, score, list(cv.n = 2, max.iter = 1000)))
  testthat::expect_s3_class(cv, "precmed")
  testthat::expect_identical(dim(precmed::abc(cv)), c(1L, 2L))
  testthat::expect_true(all(is.finite(precmed::abc(cv))))
  testthat::expect_identical(dim(cv$ate.contrastReg$ate.est.valid.group.cv), c(3L, 2L))
  testthat::expect_true(all(is.finite(cv$ate.contrastReg$ate.est.valid.group.cv)))
  invisible(capture.output(ate <- do.call(precmed::atefit, c(common, list(n.boot = 3)))))
  testthat::expect_s3_class(ate, "atefit")
  testthat::expect_identical(dim(ate$trt.boot), c(3L, 4L))
  testthat::expect_true(all(is.finite(ate$trt.boot)))
})
