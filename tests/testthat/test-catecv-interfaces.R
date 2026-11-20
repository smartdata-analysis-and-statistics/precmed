test_that("validation exposes the same methods and defaults as score fitting", {
  cases <- list(
    list(response = "binary", fun = catecvbinary, method = "logistic",
         data = simulate_binary(n = 440, seed = 1211), outcome = y ~ age),
    list(response = "count", fun = catecvcount, method = "poisson",
         data = simulate_count(n = 440, seed = 1212), outcome = y ~ age + offset(log(years))),
    list(response = "survival", fun = catecvsurv, method = "coxph",
         data = simulate_survival(n = 440, seed = 1213), outcome = survival::Surv(y, d) ~ age)
  )
  expect_identical(formals(catecv)$benefit_method, "t_learner")
  for (case in cases) {
    expect_identical(names(formals(case$fun))[1:7],
                     c("data", "outcome_model", "outcome_method", "benefit_model",
                       "benefit_method", "ps_model", "ps_method"))
    expect_identical(formals(case$fun)$outcome_method, case$method)
    expect_identical(formals(case$fun)$benefit_method, "t_learner")
    controls <- list(cv.n = 1, abc = FALSE, prop.cutoff = c(.7, 1), seed = 1214)
    if (case$response != "binary") controls$error.max <- .8
    if (case$response == "survival") controls$tau0 <- unname(quantile(case$data$y, .75))
    args <- c(list(data = case$data, outcome_model = case$outcome, ps_model = trt ~ age), controls)
    expected <- suppressWarnings(suppressMessages(do.call(case$fun, args)))
    explicit <- c(args, list(outcome_method = case$method, benefit_method = "t_learner"))
    generic <- suppressWarnings(suppressMessages(do.call(function(...) catecv(...),
      c(list(response = case$response), explicit))))
    expect_equal(generic, expected)
    key <- paste0("ate.", case$method)
    expect_true(all(is.finite(expected[[key]]$ate.est.valid.high.cv)))
    full_only <- args
    full_only$prop.cutoff <- 1
    full_fit <- suppressWarnings(suppressMessages(do.call(case$fun, full_only)))
    expect_equal(unname(full_fit[[key]]$ate.est.valid.high.cv[1L, ]),
                 full_fit$overall.ate.valid)
    expect_equal(nrow(full_fit[[key]]$ate.est.valid.low.cv), 0L)
    legacy <- c(list(case$response, case$data, case$method, case$outcome, trt ~ age), controls)
    suppressWarnings(expect_warning(old <- suppressMessages(do.call(catecv, legacy)), "score.method is deprecated"))
    expect_equal(old, expected)
    legacy <- c(list(case$data, case$method, case$outcome, trt ~ age), controls)
    suppressWarnings(expect_warning(old <- suppressMessages(do.call(case$fun, legacy)), "cate.model is deprecated"))
    expect_equal(old, expected)
    expect_error(do.call(case$fun, c(args, list(cate.model = case$outcome))),
                 "Specify only outcome_model")
  }
})

test_that("survival validation supports multiple scores and legacy censoring aliases", {
  d <- simulate_survival(n = 500, seed = 1215)
  args <- list(data = d, outcome_model = survival::Surv(y, d) ~ age + female,
               benefit_model = ~ age, ps_model = trt ~ age,
               outcome_method = c("coxph", "weibull", "log_normal", "fractional_logit"),
               benefit_method = c(rep("t_learner", 3), "two_regression"),
               tau0 = unname(quantile(d$y, .75)), cv.n = 1, error.max = .8, abc = FALSE,
               prop.cutoff = c(.7, 1), B = 1, Kfold = 2, seed = 1216)
  current <- suppressWarnings(suppressMessages(do.call(catecvsurv, c(args,
                list(censoring_model = ~ age, censoring_method = "aft (weibull)")))))
  expect_identical(rownames(current$ate.twoReg$coefficients), c("(Intercept)", "age"))
  for (method in c("coxph", "weibull", "log_normal")) {
    expect_true(all(is.finite(current[[paste0("ate.", method)]]$ate.est.valid.high.cv)))
  }
  suppressWarnings(expect_warning(old <- suppressMessages(do.call(catecv, c(list(response = "survival"), args,
    list(ipcw.model = ~ age, ipcw.method = "aft (weibull)")))), "ipcw.method is deprecated"))
  expect_equal(old, current)
  expect_s3_class(plot(current, show.abc = FALSE), "ggplot")
  expect_s3_class(boxplot(current), "ggplot")
})

test_that("survival T-learners evaluate held-out covariates with training-arm fits", {
  d <- simulate_survival(n = 400, seed = 1217)
  train <- seq_len(300)
  valid <- setdiff(seq_len(nrow(d)), train)
  tau <- unname(quantile(d$y, .75))
  x <- model.matrix(~ age, d)[, -1L, drop = FALSE]
  trt <- as.integer(d$trt == "drug1")
  fit <- precmed:::intxsurv(d$y[train], d$d[train], trt[train],
    x.cate = x[train, , drop = FALSE], x.ps = cbind(1, x[train, , drop = FALSE]),
    x.ipcw = x[train, , drop = FALSE], tau0 = tau, score.method = "coxph", x.predict = x)
  score <- precmed:::scoresurv(fit, x[valid, , drop = FALSE], tau,
                               score.method = "coxph", prediction_rows = valid)$score.coxph
  arm0 <- train[trt[train] == 0]
  arm1 <- train[trt[train] == 1]
  rmst0 <- precmed:::cox.rmst(d$y[arm0], d$d[arm0], x[arm0, , drop = FALSE], x[valid, , drop = FALSE], tau)
  rmst1 <- precmed:::cox.rmst(d$y[arm1], d$d[arm1], x[arm1, , drop = FALSE], x[valid, , drop = FALSE], tau)
  expect_equal(score, log(tau - rmst1) - log(tau - rmst0), tolerance = 1e-9)
})

test_that("count validation fits different initial outcome models within each split", {
  d <- simulate_count(n = 500, seed = 1218)
  args <- list(data = d, outcome_model = y ~ age + female + offset(log(years)),
               ps_model = trt ~ age, benefit_model = ~ scale(age),
               benefit_method = c("two_regression", "contrast_regression"),
               outcome_method = c("poisson", "gam"), xvar.smooth = "age",
               cv.n = 1, error.max = .8, abc = FALSE, prop.cutoff = c(.7, 1),
               Kfold = 2, B = 1, seed = 1219)
  direct <- suppressWarnings(suppressMessages(do.call(catecvcount, args)))
  generic <- suppressWarnings(suppressMessages(do.call(catecv, c(list(response = "count"), args))))
  expect_equal(direct, generic)
  expect_identical(generic$outcome_methods, c(twoReg = "poisson", contrastReg = "gam"))
  for (method in c("twoReg", "contrastReg"))
    expect_true(all(is.finite(generic[[paste0("ate.", method)]]$ate.est.valid.high.cv)))
})

test_that("benefit transformations preserve the held-out split's row order", {
  data <- data.frame(x = seq_len(10))
  valid <- c(9, 2, 7)
  train <- setdiff(seq_len(nrow(data)), valid)
  design <- precmed:::cv_benefit_design(data, valid, ~ scale(x), NULL, NULL)
  expected <- (data$x[valid] - mean(data$x[train])) / sd(data$x[train])
  expect_equal(unname(design$valid[, 1]), expected)
})

test_that("validation preprocessing warnings omit internal argument names", {
  count <- simulate_count(n = 100, seed = 1220)
  survival <- simulate_survival(n = 100, seed = 1221)
  messages <- character()
  withCallingHandlers({
    suppressMessages(precmed:::data.preproc("crossv", y ~ age + offset(log(years)),
      trt ~ age, count, prop.cutoff = c(0, .5, 1), prop.multi = c(0, .5, 1),
      ps.method = "glm"))
    suppressMessages(precmed:::data.preproc.surv("crossv", survival::Surv(y, d) ~ age,
      trt ~ age, data = survival, prop.cutoff = c(0, .5, 1), prop.multi = c(0, .5, 1),
      ps.method = "glm", response = "survival"))
  }, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    expect_null(conditionCall(w))
    invokeRestart("muffleWarning")
  })
  expect_equal(sum(startsWith(messages, "The first element of prop.cutoff")), 2)
  expect_true(any(startsWith(messages, "No value supplied for tau0")))
})
