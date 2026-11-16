# Run from the repository root:
# Sys.setenv(PRECMED_RUN_LOCAL_TESTS = "true", NOT_CRAN = "true")
# pkgload::load_all(".")
# testthat::test_dir("tests/local")

testthat::test_that("catefit and catecv expose cumulative hazards alongside Cox HRs", {
  testthat::skip_if(Sys.getenv("PRECMED_RUN_LOCAL_TESTS") != "true",
                    "Requires an explicit local opt-in")
  testthat::skip_if(tolower(Sys.getenv("CI")) %in% c("true", "1") ||
                      nzchar(Sys.getenv("GITHUB_ACTIONS")) ||
                      nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")),
                    "Local regression test must not run in CI or R CMD check")
  testthat::skip_on_cran()

  data(survivalExample, package = "precmed")
  common <- list(data = survivalExample,
                 cate.model = survival::Surv(y, d) ~ age + female + previous_cost + previous_number_relapses,
                 ps.model = trt ~ age + previous_treatment, tau0 = 522.66,
                 rmst.control = list(solver = "R"))
  options <- list(score.method = "contrastReg", initial.predictor.method = "logistic",
                  prop.cutoff = c(0.5, 1), contrast.control = list(solver = "R"), seed = 7)
  ate <- suppressWarnings(do.call(atefitsurv, c(common, list(n.boot = 2, seed = 1))))
  fit <- suppressWarnings(do.call(catefitsurv, c(common, options)))
  cv <- suppressWarnings(do.call(catecvsurv, c(common, options,
                     list(prop.multi = c(0, 0.5, 1), cv.n = 2, abc = FALSE))))
  testthat::expect_equal(unname(fit$chr.contrastReg["prop1"]),
               exp(ate$log.cumulative.hazard.ratio$estimate))
  testthat::expect_length(fit$hr.contrastReg, 2)
  testthat::expect_equal(unname(cv$chr.contrastReg$chr.est.valid.high.cv["prop1", ]),
               cv$overall.chr.valid)
  testthat::expect_equal(dim(cv$chr.contrastReg$chr.est.valid.group.cv), c(2L, 2L))
  testthat::expect_equal(dim(cv$hr.contrastReg$hr.est.valid.group.cv), c(2L, 2L))
})

testthat::test_that("single-predictor count boosting works in fitting and validation", {
  testthat::skip_if(Sys.getenv("PRECMED_RUN_LOCAL_TESTS") != "true",
                    "Requires an explicit local opt-in")
  testthat::skip_if(tolower(Sys.getenv("CI")) %in% c("true", "1") ||
                      nzchar(Sys.getenv("GITHUB_ACTIONS")) ||
                      nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")),
                    "Local regression test must not run in CI or R CMD check")
  testthat::skip_on_cran()

  set.seed(81)
  n <- 500L
  dat <- data.frame(x = rnorm(n), trt = rep(0:1, each = n / 2),
                    years = runif(n, 0.3, 2))
  dat$y <- rpois(n, exp(0.4 + 0.25 * dat$x - 0.2 * dat$trt) * dat$years)
  args <- list(data = dat, score.method = c("boosting", "twoReg", "contrastReg"),
               cate.model = y ~ x + offset(log(years)), ps.model = trt ~ x,
               initial.predictor.method = "boosting", B = 1, Kfold = 2,
               n.trees.boosting = 30, plot.gbmperf = FALSE, n.cores = 1,
               seed = 51)
  messages <- character()
  fit <- withCallingHandlers(do.call(catefitcount, args), warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  for (method in args$score.method) {
    testthat::expect_true(all(is.finite(fit[[paste0("score.", method)]])))
  }
  testthat::expect_false(any(grepl("OOB|out.of.bag", messages, ignore.case = TRUE)))
  messages <- character()
  cv <- withCallingHandlers(do.call(catecvcount, c(args, list(cv.n = 1,
    abc = FALSE, error.max = 0.5))), warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  for (method in args$score.method) {
    testthat::expect_true(all(is.finite(cv[[paste0("ate.", method)]]$ate.est.valid.high.cv)))
  }
  testthat::expect_false(any(grepl("OOB|out.of.bag", messages, ignore.case = TRUE)))
  args$score.method <- "boosting"
  args$seed <- 11
  messages <- character()
  tied <- withCallingHandlers(do.call(catecvcount, c(args, list(cv.n = 1,
    abc = TRUE, error.max = 0.5))), warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  testthat::expect_true(any(grepl("score.boosting, validation iteration 1, low subgroup", messages, fixed = TRUE)))
  testthat::expect_true(anyNA(tied$ate.boosting$ate.est.valid.low.cv))
  testthat::expect_true(any(is.finite(tied$ate.boosting$ate.est.valid.high.cv)))
})

testthat::test_that("count CV records convergence and support without verbose output", {
  testthat::skip_if(Sys.getenv("PRECMED_RUN_LOCAL_TESTS") != "true",
                    "Requires an explicit local opt-in")
  testthat::skip_if(tolower(Sys.getenv("CI")) %in% c("true", "1") ||
                      nzchar(Sys.getenv("GITHUB_ACTIONS")) ||
                      nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")),
                    "Local regression test must not run in CI or R CMD check")
  testthat::skip_on_cran()

  data(countExample)
  fit <- catecvcount(data = countExample[1:500, ],
    cate.model = y ~ age + previous_treatment + offset(log(years)),
    ps.model = trt ~ age + previous_treatment,
    score.method = "contrastReg", initial.predictor.method = "poisson",
    B = 1, Kfold = 2, cv.n = 1, abc = FALSE, verbose = 0, seed = 45)
  testthat::expect_length(fit$ate.contrastReg$converge.contrastReg.cv, 1)
  testthat::expect_length(fit$ate.contrastReg$converge.repetitions.cv[[1]], 1)
  testthat::expect_true(nrow(fit$subgroup.support) > 0)
  testthat::expect_true(all(c("training", "validation") %in% fit$subgroup.support$split))
  testthat::expect_s3_class(diagnostics(fit), "data.frame")
})

testthat::test_that("simulated rate ratio and exposure match their parameters", {
  testthat::skip_if(Sys.getenv("PRECMED_RUN_LOCAL_TESTS") != "true",
                    "Requires an explicit local opt-in")
  testthat::skip_if(tolower(Sys.getenv("CI")) %in% c("true", "1") ||
                      nzchar(Sys.getenv("GITHUB_ACTIONS")) ||
                      nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")),
                    "Local regression test must not run in CI or R CMD check")
  testthat::skip_on_cran()

  d <- simulate_count(n = 50000, seed = 32, conditional_rr = 0.6,
    ps_weights = c(age = 0, drugB = 0, drugC = 0),
    lp_weights = c(age = 0, female = 0, previous_cost = 0,
                   previous_number_relapses = 0, treatment_age = 0, treatment_relapses = 0))
  rate0 <- with(d, sum(y[trt == "drug0"]) / sum(years[trt == "drug0"]))
  rate1 <- with(d, sum(y[trt == "drug1"]) / sum(years[trt == "drug1"]))
  testthat::expect_equal(rate0, 0.5, tolerance = 0.03)
  testthat::expect_equal(rate1 / rate0, 0.6, tolerance = 0.04)
})
