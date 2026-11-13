# Expensive LASSO fits and allocation benchmarks are opt-in and local only.
testthat::test_that("native PS screening preserves selected interactions and propensity scores", {
  testthat::skip_if(Sys.getenv("PRECMED_RUN_LOCAL_TESTS") != "true", "Requires local opt-in")
  testthat::skip_if(tolower(Sys.getenv("CI")) %in% c("true", "1") ||
                      nzchar(Sys.getenv("GITHUB_ACTIONS")) ||
                      nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")), "Local tests only")
  testthat::skip_on_cran()
  make_data <- function(seed, n = 600L, p = 8L) {
    set.seed(seed)
    x <- matrix(rnorm(n*p), n, p)
    colnames(x) <- paste0("x", seq_len(p))
    trt <- rbinom(n, 1, plogis(.3*x[, 1]))
    list(x = x, trt = trt)
  }
  compare_screen <- function(x, trt, xnew = x[1:17, , drop = FALSE]) {
    r <- precmed:::ps_interactions(x, trt, xnew, solver = "R")
    native <- precmed:::ps_interactions(x, trt, xnew, solver = "C")
    testthat::expect_identical(native, r)
    testthat::expect_identical(precmed:::ps_interactions(x, trt, xnew), native)
    native
  }
  compare_fit <- function(x, trt, reference = original_glm_ps, xnew = x[1:17, , drop = FALSE]) {
    args <- list(trt = trt, x.ps = cbind(1, x),
                  xnew = if (!is.null(xnew)) cbind(1, xnew))
    set.seed(42); before <- do.call(reference, args); before.rng <- .Random.seed
    set.seed(42); r <- do.call(precmed:::glm.ps, c(args, list(solver = "R")))
    testthat::expect_identical(r, before)
    testthat::expect_identical(.Random.seed, before.rng)
    set.seed(42); native <- do.call(precmed:::glm.ps, c(args, list(solver = "C")))
    testthat::expect_identical(native, before)
    testthat::expect_identical(.Random.seed, before.rng)
    testthat::expect_length(native, if (is.null(xnew)) nrow(x) else nrow(xnew))
  }
  corrected <- original_glm_ps
  code <- paste(deparse(body(corrected), width.cutoff = 500L), collapse = "\n")
  code <- gsub("xmat[trt == 1, ]", "xmat[trt == 1, , drop = FALSE]", code, fixed = TRUE)
  code <- gsub("xmat[trt == 0, ]", "xmat[trt == 0, , drop = FALSE]", code, fixed = TRUE)
  code <- gsub("xmatnew[, (diff > 0.1 & sigma1 + sigma0 > 0)]",
               "xmatnew[, (diff > 0.1 & sigma1 + sigma0 > 0), drop = FALSE]", code, fixed = TRUE)
  body(corrected) <- parse(text = code)[[1L]]
  for (seed in c(1L, 42L, 999L)) {
    z <- make_data(seed)
    compare_screen(z$x, z$trt)
    compare_fit(z$x, z$trt)
  }
  z <- make_data(19L)
  compare_screen(cbind(z$x, constant = 1, zero = 0, alias = z$x[, 1]), z$trt)
  compare_screen(z$x * 1e5 + 1e8, z$trt)
  compare_fit(z$x, z$trt, xnew = NULL)
  compare_fit(z$x, z$trt, reference = corrected, xnew = z$x[1, , drop = FALSE])
  # Two predictors previously produced a vector in apply(): change only drop
  # handling in the frozen reference to define the intended old calculation.
  set.seed(7)
  trt <- rep(0:1, each = 200)
  u <- rnorm(400)
  v <- ifelse(trt == 1, u, -u) + rnorm(400, sd = .2)
  two <- cbind(u = u, v = v)
  testthat::expect_error(original_glm_ps(trt, cbind(1, two)), "dim")
  testthat::expect_equal(ncol(compare_screen(two, trt)$pairs), 1L)
  compare_fit(two, trt, reference = corrected)
  # Exactly two predictors and no retained interaction, plus one/zero predictors.
  two <- cbind(rep(seq(-2, 2, length.out = 200), 2), rep(seq(2, -2, length.out = 200), 2))
  testthat::expect_equal(ncol(compare_screen(two, trt)$pairs), 0L)
  compare_fit(two, trt, reference = corrected)
  compare_fit(two[, 1, drop = FALSE], trt)
  zero <- matrix(numeric(), length(trt), 0)
  testthat::expect_equal(precmed:::glm.ps(trt, cbind(1, zero)), rep(mean(trt), length(trt)))

  # Locate the strict SMD threshold and compare decisions on both sides and
  # within native rounding tolerance. These candidates are rechecked in R.
  n <- 100L
  trt <- rep(0:1, each = n)
  values <- seq(-1, 1, length.out = n)
  shift <- .1 * sqrt(2 * stats::var(values) + 1e-6) / sqrt(2)
  for (multiplier in c(1 - 1e-10, 1, 1 + 1e-10)) {
    x <- cbind(rep(1, 2*n), c(values, values + shift*multiplier))
    compare_screen(x, trt)
    compare_screen(cbind(1, x[, 2] + 1e12), trt)
  }
  testthat::expect_error(precmed:::resolve_ps_control(list(solvr = "C")), "named list")
  testthat::expect_error(precmed:::ps_interactions(matrix(1, 3, 2), c(0, 0, 1)), "at least two")
  testthat::expect_error(precmed:::ps_interactions(matrix(NA_real_, 4, 2), c(0, 0, 1, 1)), "finite")

  # Public APIs must forward ps.control, independently of contrast/RMST options.
  z <- make_data(29L, n = 800L, p = 3L)
  set.seed(19)
  count <- data.frame(y = rpois(800, exp(.2*z$x[,1] - .2*z$trt)), trt = z$trt, z$x)
  event <- rexp(800, exp(.2*z$x[,1] - .2*z$trt)); censor <- rexp(800, .5)
  surv <- data.frame(y = pmin(event, censor), d = as.numeric(event <= censor), trt = z$trt, z$x)
  for (response in c("count", "survival")) {
    data <- switch(response, count = count, survival = surv)
    formula <- if (response == "survival") survival::Surv(y, d) ~ x1 + x2 else y ~ x1 + x2 + offset(time)
    if (response == "count") data$time <- 0
    common <- list(response = response, data = data, cate.model = formula,
                   ps.model = trt ~ x1 + x2 + x3, ps.method = "lasso", seed = 42)
    if (response == "survival") common <- c(common, list(ipcw.model = ~ x1 + x2, tau0 = 1.5))
    for (fun in list(precmed::catefit, precmed::catecv, precmed::atefit)) {
      extra <- if (identical(fun, precmed::atefit)) list(n.boot = 2) else
        list(score.method = "poisson", prop.cutoff = c(.5, 1))
      if (identical(fun, precmed::catecv)) extra <- c(extra, list(cv.n = 1, max.iter = 100))
      invisible(capture.output(before <- do.call(fun, c(common, extra, list(ps.control = list(solver = "R"))))))
      invisible(capture.output(after <- do.call(fun, c(common, extra, list(ps.control = list(solver = "C"))))))
      testthat::expect_identical(after, before)
    }
  }
  # Screening-only allocation/time comparison; no wall-clock assertion.
  z <- make_data(7L, n = 2000L, p = 60L)
  xnew <- z$x[1:1000, , drop = FALSE]
  selected <- compare_screen(z$x, z$trt, xnew)
  # Extract only the screening stage from the immutable original function.
  old.screen <- original_glm_ps
  expressions <- as.list(body(old.screen))
  cut <- which(vapply(expressions, identical, logical(1),
                       quote(xmat <- cbind(x, xmat))))[1]
  body(old.screen) <- as.call(c(expressions[seq_len(cut - 1)], list(quote(
    return(list(train = xmat, prediction = if (!is.null(xnew)) xmatnew))
  ))))
  old.x <- cbind(1, z$x); old.new <- cbind(1, xnew)
  testthat::expect_identical(old.screen(z$trt, old.x, old.new),
                             selected[c("train", "prediction")])
  measure <- function(solver) {
    file <- tempfile(); Rprofmem(file)
    time <- system.time(invisible(if (solver == "original")
      old.screen(z$trt, old.x, old.new) else
      precmed:::ps_interactions(z$x, z$trt, xnew, solver)))[["elapsed"]]
    Rprofmem(NULL)
    bytes <- suppressWarnings(as.numeric(sub(" .*", "", readLines(file))))
    unlink(file)
    c(seconds = time, allocated = sum(bytes, na.rm = TRUE))
  }
  r <- apply(replicate(3, measure("original")), 1, median)
  native <- apply(replicate(3, measure("C")), 1, median)
  cat(sprintf("\nPS screening retained %d of %d pairs\n", ncol(selected$pairs), choose(60, 2)))
  cat(sprintf("R %.3fs / %.1f MB allocated; C %.3fs / %.1f MB allocated\n",
              r[1], r[2]/1e6, native[1], native[2]/1e6))
  testthat::expect_lt(native[2], r[2]/10)
})
