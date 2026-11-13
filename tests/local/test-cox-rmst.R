# Expensive comparisons are opt-in and excluded from CI and CRAN builds.
testthat::test_that("native RMST agrees with R and avoids dense survival matrices", {
  testthat::skip_if(Sys.getenv("PRECMED_RUN_LOCAL_TESTS") != "true", "Requires local opt-in")
  testthat::skip_if(tolower(Sys.getenv("CI")) %in% c("true", "1") ||
                      nzchar(Sys.getenv("GITHUB_ACTIONS")) ||
                      nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")), "Local tests only")
  testthat::skip_on_cran()

  # Deterministic Cox output isolates integration boundaries from model fitting.
  env <- new.env(parent = asNamespace("precmed"))
  env$coxph <- function(...) list(coef = c(0, 0))
  env$coxph.detail <- function(...) list(time = c(1, 2, 4), hazard = c(.1, .2, .3))
  original <- original_cox_rmst
  revised <- precmed:::cox.rmst
  environment(original) <- environment(revised) <- env
  x <- cbind(x1 = 1:6, x2 = (1:6)^2)
  rownames(x) <- letters[1:6]
  args <- list(y = 1:6, d = rep(0:1, 3), x.cate = x, xnew = x, tau0 = .5)
  for (tau in c(.5, 1, 2, 3, 4, 8)) {
    args$tau0 <- tau
    before <- do.call(original, args)
    testthat::expect_identical(do.call(revised, c(args, list(solver = "R"))), before)
    testthat::expect_identical(do.call(revised, args),
                               do.call(revised, c(args, list(solver = "C"))))
    testthat::expect_equal(do.call(revised, c(args, list(solver = "C"))),
                           before, tolerance = 1e-14)
  }
  args$tau0 <- .5
  testthat::expect_identical(do.call(revised, c(args, list(solver = "C"))), rep(.5, 6))
  args$xnew <- args$xnew[1, , drop = FALSE]
  testthat::expect_equal(do.call(revised, c(args, list(solver = "C"))),
                         do.call(original, args), tolerance = 1e-14)

  # Exercise exact zeros, ones, zero-width intervals and extreme hazard ratios.
  native <- function(s, h, g) .Call(precmed:::precmed_rmst, as.double(s), as.double(h), as.double(g))
  dense <- function(s, h, g) {
    surv <- crossprod(t(rep(1, length(h))), s)^h
    colSums(t(cbind(1, surv)) * g)
  }
  testthat::expect_equal(native(c(1, .5, 0), c(0, .01, 1, 1000, Inf), c(1, 0, 1, 2)),
                         dense(c(1, .5, 0), c(0, .01, 1, 1000, Inf), c(1, 0, 1, 2)),
                         tolerance = 1e-14)
  testthat::expect_identical(native(numeric(), c(1, 2), 3), c(3, 3))
  testthat::expect_identical(native(.5, numeric(), c(1, 2)), numeric())
  testthat::expect_error(native(.5, 1, 1), "interval lengths")
  testthat::expect_error(native(.5, NA_real_, c(1, 2)), "hazard ratios")
  testthat::expect_error(precmed:::resolve_rmst_control(list(solvr = "C")), "named list")
  testthat::expect_error(precmed:::resolve_rmst_control(list(solver = "other")), "arg")
  testthat::expect_error(precmed::atefit(response = "count", data = data.frame(),
    cate.model = y ~ x, ps.model = trt ~ 1, rmst.control = list(solver = "C")), "only for survival")

  make_data <- function(seed, n = 400L) {
    set.seed(seed)
    x <- cbind(x1 = rnorm(n), x2 = rnorm(n))
    trt <- rbinom(n, 1, .5)
    event <- rexp(n, exp(.2 * x[, 1] - .3 * trt))
    censor <- rexp(n, .5)
    list(y = pmax(.001, round(pmin(event, censor), 3)),
         d = as.numeric(event <= censor), trt = trt, x = x)
  }
  max.diff <- 0
  for (seed in c(1L, 42L, 999L)) {
    z <- make_data(seed)
    args <- list(y = z$y, d = z$d, x.cate = z$x,
                 xnew = z$x[1:17, , drop = FALSE], tau0 = 1.5)
    # Constant and aliased predictors exercise existing column filtering.
    for (design in list(z$x, cbind(z$x, const = 1), cbind(z$x, duplicate = z$x[, 1]))) {
      args$x.cate <- design
      args$xnew <- design[1:17, , drop = FALSE]
      before <- do.call(original_cox_rmst, args)
      r <- do.call(precmed:::cox.rmst, c(args, list(solver = "R")))
      c <- do.call(precmed:::cox.rmst, c(args, list(solver = "C")))
      testthat::expect_identical(r, before)
      testthat::expect_equal(c, r, tolerance = 1e-14)
      max.diff <- max(max.diff, abs(c - r))
    }
  }

  # Public interfaces: bootstrap, subgroup estimates and repeated validation.
  z <- make_data(7L)
  data <- data.frame(y = z$y, d = z$d, trt = z$trt, z$x)
  common <- list(response = "survival", data = data,
                 cate.model = survival::Surv(y, d) ~ x1 + x2,
                 ps.model = trt ~ x1 + x2, ipcw.model = ~ x1 + x2, tau0 = 1.5,
                 seed = 42)
  fit_args <- c(common, list(score.method = "poisson", prop.cutoff = c(.5, 1)))
  before <- do.call(precmed::catefit, c(fit_args, list(rmst.control = list(solver = "R"))))
  after <- do.call(precmed::catefit, c(fit_args, list(rmst.control = list(solver = "C"))))
  testthat::expect_equal(after, before, tolerance = 1e-12)
  cv_args <- c(fit_args, list(cv.n = 2, max.iter = 100))
  before <- do.call(precmed::catecv, c(cv_args, list(rmst.control = list(solver = "R"))))
  after <- do.call(precmed::catecv, c(cv_args, list(rmst.control = list(solver = "C"))))
  testthat::expect_equal(after, before, tolerance = 1e-12)
  invisible(capture.output(before <- do.call(precmed::atefit, c(common, list(n.boot = 3, rmst.control = list(solver = "R"))))))
  invisible(capture.output(after <- do.call(precmed::atefit,
    c(common, list(n.boot = 3, rmst.control = list(solver = "C"))))))
  testthat::expect_equal(after, before, tolerance = 1e-12)

  testthat::expect_identical(precmed:::resolve_rmst_control(list())$solver, "C")
  testthat::expect_equal(do.call(precmed::catefit, fit_args),
    do.call(precmed::catefit, c(fit_args, list(rmst.control = list(solver = "C")))))

  # Integration-only benchmark, with allocated bytes measured (not peak RSS).
  s <- exp(-seq(.001, 5, length.out = 1000))
  h <- seq(.1, 3, length.out = 5000)
  g <- rep(.001, 1001)
  testthat::expect_equal(native(s, h, g), dense(s, h, g), tolerance = 1e-14)
  measure <- function(fun) {
    file <- tempfile()
    Rprofmem(file)
    timing <- system.time(invisible(fun(s, h, g)))[["elapsed"]]
    Rprofmem(NULL)
    entries <- readLines(file)
    unlink(file)
    bytes <- suppressWarnings(as.numeric(sub(" .*", "", entries)))
    c(seconds = timing, allocated = sum(bytes, na.rm = TRUE))
  }
  r <- apply(replicate(3, measure(dense)), 1, median)
  c <- apply(replicate(3, measure(native)), 1, median)
  cat(sprintf("\nLargest RMST difference: %.3g\n", max.diff))
  cat(sprintf("Integration: R %.3fs / %.1f MB allocated; C %.3fs / %.1f KB allocated\n",
              r[1], r[2]/1e6, c[1], c[2]/1e3))
  testthat::expect_lt(c[2], r[2] / 10)
})
