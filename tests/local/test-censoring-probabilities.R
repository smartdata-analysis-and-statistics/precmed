# Opt-in only; excluded from source packages and normal test discovery.
testthat::test_that("cumulative hazard lookup preserves censoring probabilities", {
  testthat::skip_if(Sys.getenv("PRECMED_RUN_LOCAL_TESTS") != "true",
                    "Requires an explicit local opt-in")
  testthat::skip_if(tolower(Sys.getenv("CI")) %in% c("true", "1") ||
                      nzchar(Sys.getenv("GITHUB_ACTIONS")) ||
                      nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")),
                    "Local regression test must not run in CI or R CMD check")
  testthat::skip_on_cran()

  differences <- numeric()
  exact <- logical()
  compare <- function(args, original = original_ipcw_surv,
                      revised = precmed:::ipcw.surv) {
    before <- do.call(original, args)
    after <- do.call(revised, args)
    testthat::expect_equal(after, before, tolerance = 1e-14)
    testthat::expect_identical(attributes(after), attributes(before))
    differences <<- c(differences, max(abs(after - before)))
    exact <<- c(exact, identical(after, before))
    after
  }

  # Isolate lookup boundaries independently of the fitted Cox model.
  y <- c(0.5, 1, 1.5, 2, 3, 4, 8, 2)
  x <- cbind(x1 = seq_along(y), x2 = seq_along(y)^2)
  mock_env <- new.env(parent = asNamespace("precmed"))
  mock_env$coxph <- function(...) list(coef = c(0, 0))
  mock_env$coxph.detail <- function(...) list(time = c(1, 2, 4),
                                            hazard = c(0.1, 0.2, 0.3))
  original <- original_ipcw_surv
  revised <- precmed:::ipcw.surv
  environment(original) <- environment(revised) <- mock_env
  args <- list(y = y, d = rep(0:1, 4), x.ipcw = x, tau0 = 3, surv.min = 0)
  actual <- compare(args, original, revised)
  testthat::expect_equal(actual, exp(-c(0, .1, .1, .3, .3, .3, .3, .3)),
                         tolerance = 1e-14)
  args$tau0 <- 10
  compare(args, original, revised)
  args$surv.min <- 0.9
  testthat::expect_true(all(compare(args, original, revised) >= 0.9))
  mock_env$coxph.detail <- function(...) list(time = numeric(), hazard = numeric())
  testthat::expect_identical(compare(args, original, revised), rep(1, length(y)))

  # Real fitted Cox models: tied observations, truncation and external follow-up.
  for (seed in c(11L, 29L, 53L)) {
    set.seed(seed)
    n <- 500L
    x <- cbind(x1 = rnorm(n), x2 = rbinom(n, 1, .4))
    event <- rexp(n, exp(.2 * x[, 1]))
    censor <- rexp(n, exp(-.2 * x[, 1] + .1 * x[, 2]))
    y <- pmax(.01, round(pmin(event, censor), 2))
    args <- list(y = y, d = as.numeric(event <= censor), x.ipcw = x,
                 tau0 = unname(quantile(y, .8)))
    compare(args)
    args$yf <- censor
    compare(args)
    args$ipcw.method <- "aft (weibull)"
    compare(args)
  }

  # Whole-function timing; comparison assertions do not depend on runtime.
  set.seed(19)
  n <- 5000L
  x <- cbind(x1 = rnorm(n), x2 = rnorm(n))
  event <- rexp(n)
  censor <- rexp(n)
  args <- list(y = pmin(event, censor), d = as.numeric(event <= censor),
               x.ipcw = x, tau0 = 2)
  compare(args)
  timings <- replicate(3, c(
    original = system.time(do.call(original_ipcw_surv, args))[["elapsed"]],
    revised = system.time(do.call(precmed:::ipcw.surv, args))[["elapsed"]]
  ))
  cat(sprintf("\nCensoring probabilities: %d/%d comparisons bitwise identical; max difference %.3g\n",
              sum(exact), length(exact), max(differences)))
  cat(sprintf("Median whole-function time: original %.3fs; revised %.3fs (%.2fx)\n",
              median(timings[1, ]), median(timings[2, ]),
              median(timings[1, ]) / median(timings[2, ])))
})
