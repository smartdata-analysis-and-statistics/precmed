test_that("ATE outcome_model replaces the deprecated alias without changing results", {
  cases <- list(
    list(response = "binary", data = simulate_binary(n = 400, seed = 10), formula = y ~ age),
    list(response = "count", data = simulate_count(n = 400, seed = 11), formula = y ~ age + offset(log(years))),
    list(response = "survival", data = simulate_survival(n = 400, seed = 12), formula = survival::Surv(y, d) ~ age)
  )
  for (case in cases) {
    args <- list(response = case$response, data = case$data, ps_model = trt ~ age,
                 n.boot = 3, seed = 13)
    if (case$response == "survival") args$tau0 <- median(case$data$y) * 1.1
    new <- suppressWarnings(do.call(atefit, c(args, list(outcome_model = case$formula))))
    expect_warning(old <- do.call(atefit, c(args, list(cate.model = case$formula))),
                   "cate.model is deprecated", info = case$response)
    expect_equal(new, old, info = case$response)
    expect_error(do.call(atefit, c(args, list(outcome_model = case$formula, cate.model = case$formula))),
                 "Specify only outcome_model", info = case$response)
    expect_error(do.call(atefit, args), "outcome_model must be supplied", info = case$response)
  }
})

test_that("direct ATE interfaces accept the alias and positional formulas", {
  d <- simulate_binary(n = 400, seed = 20)
  args <- list(data = d, ps_model = trt ~ age, n.boot = 3, seed = 21)
  new <- do.call(atefitbinary, c(args, list(outcome_model = y ~ age)))
  expect_warning(old <- do.call(atefitbinary, c(args, list(cate.model = y ~ age))), "deprecated")
  expect_equal(new, old)
  expect_equal(new, atefitbinary(d, y ~ age, trt ~ age, n.boot = 3, seed = 21))
  expect_error(atefitbinary(d, y ~ age, trt ~ age, cate.model = y ~ age), "Specify only")
})
test_that("atefit accepts the deprecated alias through dots", {
  expect_false("cate.model" %in% names(formals(atefit)))
  expect_true("..." %in% names(formals(atefit)))
  d <- simulate_binary(n = 350, seed = 421)
  expect_warning(old <- atefit("binary", d, cate.model = y ~ age,
                               ps_model = trt ~ age, n.boot = 3, seed = 422),
                 "cate.model is deprecated")
  current <- atefit("binary", d, outcome_model = y ~ age,
                    ps_model = trt ~ age, n.boot = 3, seed = 422)
  expect_equal(old, current)
  expect_error(atefit("binary", d, outcome_model = y ~ age,
                      cate.model = y ~ age, ps_model = trt ~ age), "Specify only")
  expect_error(atefit("binary", d, y ~ age, ps_model = trt ~ age,
                      unknown_argument = TRUE), "unused argument")
})

test_that("outcome-specific interfaces accept the alias through dots", {
  cases <- list(
    list(fun = atefitbinary, data = simulate_binary(n = 400, seed = 31), formula = y ~ age),
    list(fun = atefitcount, data = simulate_count(n = 400, seed = 32), formula = y ~ age + offset(log(years))),
    list(fun = atefitsurv, data = simulate_survival(n = 400, seed = 33), formula = survival::Surv(y, d) ~ age)
  )
  for (case in cases) {
    expect_false("cate.model" %in% names(formals(case$fun)))
    expect_true("..." %in% names(formals(case$fun)))
    args <- list(data = case$data, ps_model = trt ~ age, n.boot = 3, seed = 34)
    if (identical(case$fun, atefitsurv)) args$tau0 <- median(case$data$y) * 1.1
    current <- suppressMessages(do.call(case$fun, c(args, list(outcome_model = case$formula))))
    expect_warning(old <- suppressMessages(do.call(case$fun, c(args, list(cate.model = case$formula)))), "cate.model is deprecated")
    expect_equal(old, current)
    expect_error(do.call(case$fun, c(args, list(outcome_model = case$formula, cate.model = case$formula))), "Specify only")
    expect_error(do.call(case$fun, c(args, list(outcome_model = case$formula, unknown = TRUE))), "unused argument")
  }
  mean_fun <- get("atefitmean", envir = asNamespace("precmed"))
  expect_false("cate.model" %in% names(formals(mean_fun)))
  expect_true("..." %in% names(formals(mean_fun)))
  expect_warning(expect_error(mean_fun(data.frame(), cate.model = y ~ x), "not implemented"), "cate.model is deprecated")
  expect_error(mean_fun(data.frame(), outcome_model = y ~ x, cate.model = y ~ x), "Specify only")
  expect_error(mean_fun(data.frame(), outcome_model = y ~ x, unknown = TRUE), "unused argument")
})

test_that("ATE treatment-model aliases warn and preserve results", {
  d <- simulate_binary(n = 400, seed = 70)
  for (fun in list(atefit, atefitbinary, atefitcount, atefitsurv,
                   get("atefitmean", envir = asNamespace("precmed")))) {
    f <- names(formals(fun))
    expect_false(any(c("ps.model", "ps.method") %in% f))
    expect_identical(match("ps_method", f), match("ps_model", f) + 1L)
    base <- list(data = d, outcome_model = y ~ age)
    if (identical(fun, atefit)) base$response <- "binary"
    expect_error(do.call(fun, c(base, list(ps_model = trt ~ 1, ps.model = trt ~ 1))), "Specify only ps_model")
    expect_error(do.call(fun, c(base, list(ps_method = "glm", ps.method = "glm"))), "Specify only ps_method")
  }
  args <- list(response = "binary", data = d, outcome_model = y ~ age,
               n.boot = 3, seed = 71)
  current <- do.call(atefit, c(args, list(ps_model = trt ~ age, ps_method = "glm")))
  expect_warning(old_model <- do.call(atefit, c(args, list(ps.model = trt ~ age))), "ps.model is deprecated")
  expect_equal(old_model, current)
  expect_warning(old_method <- do.call(atefit, c(args, list(ps_model = trt ~ age, ps.method = "glm"))), "ps.method is deprecated")
  expect_equal(old_method, current)
  direct <- list(data = d, outcome_model = y ~ age, n.boot = 3, seed = 71)
  expect_warning(old_direct <- do.call(atefitbinary, c(direct, list(ps.model = trt ~ age))), "ps.model is deprecated")
  expect_equal(old_direct, current)
})

test_that("survival ATE censoring aliases preserve named and positional calls", {
  d <- simulate_survival(n = 400, seed = 121)
  args <- list(data = d, outcome_model = survival::Surv(y, d) ~ age,
               ps_model = trt ~ age, tau0 = median(d$y), n.boot = 3, seed = 122)
  for (fun in list(atefit, atefitsurv)) {
    f <- names(formals(fun))
    expect_false(any(c("ipcw.model", "ipcw.method") %in% f))
    expect_identical(match("censoring_method", f), match("censoring_model", f) + 1L)
    base <- args
    if (identical(fun, atefit)) base$response <- "survival"
    current <- suppressMessages(do.call(fun, c(base, list(censoring_model = ~ age,
                                                         censoring_method = "breslow"))))
    expect_warning(old_model <- suppressMessages(do.call(fun, c(base, list(ipcw.model = ~ age,
                                                                         censoring_method = "breslow")))),
                   "ipcw.model is deprecated")
    expect_equal(old_model, current)
    expect_warning(old_method <- suppressMessages(do.call(fun, c(base, list(censoring_model = ~ age,
                                                                          ipcw.method = "breslow")))),
                   "ipcw.method is deprecated")
    expect_equal(old_method, current)
    expect_error(do.call(fun, c(base, list(censoring_model = ~ age, ipcw.model = ~ age))), "Specify only censoring_model")
    expect_error(do.call(fun, c(base, list(censoring_method = "breslow", ipcw.method = "breslow"))), "Specify only censoring_method")
  }
  positional <- suppressMessages(atefitsurv(d, survival::Surv(y, d) ~ age, "coxph",
    trt ~ age, "glm", ~ age, "breslow", tau0 = median(d$y), n.boot = 3, seed = 122))
  expect_equal(positional, suppressMessages(do.call(atefitsurv, c(args, list(censoring_model = ~ age)))))
  aft <- suppressMessages(do.call(atefitsurv, c(args, list(censoring_method = "aft (weibull)"))))
  expect_warning(old_aft <- suppressMessages(do.call(atefitsurv, c(args, list(ipcw.method = "aft (weibull)")))),
                 "ipcw.method is deprecated")
  expect_equal(old_aft, aft)
})
