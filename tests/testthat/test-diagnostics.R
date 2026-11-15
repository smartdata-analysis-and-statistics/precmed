test_that("survival diagnostics produce one row per recorded message", {
  x <- structure(list(response = "survival", `errors/warnings` = list(
    randomForest = list(est.valid.low.cv = list(cv1 = list(
      warnings = list(`prop 0.2` = c("Cumulative hazard ratio is undefined", "Another warning")),
      errors = list(`prop 0.4` = "Model failed"))),
      est.train.group.cv = list(cv2 = list(errors = list(`prop 0.7` = "Fit failed")))))),
    class = "precmed")
  result <- diagnostics(x)
  expect_equal(nrow(result), 4)
  expect_equal(result$split, c("validation", "validation", "validation", "training"))
  expect_equal(result$iteration, c(1L, 1L, 1L, 2L))
  expect_equal(result$subgroup, c("prop 0.4", "prop 0.2", "prop 0.2", "prop 0.7"))
  expect_equal(result$type, c("error", "warning", "warning", "error"))
  expect_equal(result$measure[2], "cumulative hazard ratio")
  expect_true(is.na(result$measure[1]))
})

test_that("diagnostics reports stored atefit recoding information", {
  x <- structure(list(response = "survival", warning = "Recoded drug0 to 0.\n"), class = "atefit")
  result <- diagnostics(x)
  expect_equal(nrow(result), 1)
  expect_equal(result$type, "info")
  expect_equal(result$message, "Recoded drug0 to 0.")
  expect_true(all(is.na(result[c("method", "split", "iteration", "group", "subgroup", "measure")])))
  x$warning <- NULL
  expect_equal(nrow(diagnostics(x)), 0)
})

test_that("diagnostics handles whole-sample fits and empty results", {
  x <- structure(list(response = "survival", `errors/warnings` = list(
    contrastReg = list(est.high = list(warnings = list(`prop 1.0` = "Warning"))))), class = "catefit")
  result <- diagnostics(x)
  expect_equal(result$group, "high")
  expect_true(is.na(result$split) && is.na(result$iteration))
  x[["errors/warnings"]]$contrastReg$est.high$warnings <- NULL
  expect_equal(nrow(diagnostics(x)), 0)
  expect_type(diagnostics(x)$iteration, "integer")
  expect_error(diagnostics(structure(list(response = "continuous"), class = "precmed")), "count or survival")
})

test_that("count subgroup diagnostics record support and preserve warnings", {
  d <- precmed:::count_diagnostic_collector()
  callback <- d$callback("poisson", "validation", 1L, c(2, 0, 1), c(0, 0, 1))
  expect_warning(value <- precmed:::count_subgroup_log_effect(
    c(TRUE, TRUE, FALSE), c(2, 0, 1), cbind(x = 1:3), cbind(1, 1:3),
    rep(0, 3), c(0, 0, 1), "glm", .01, .99, list(), "test", callback, "high", "0.5"),
    "treatment effect unavailable")
  expect_true(is.na(value))
  expect_equal(d$support$n0, 2L)
  expect_equal(d$support$n1, 0L)
  expect_equal(d$support$events0, 2)
  x <- structure(list(response = "count", diagnostic.messages = d$messages), class = "precmed")
  expect_equal(diagnostics(x)$subgroup, "0.5")
  expect_equal(diagnostics(x)$type, "warning")
})

test_that("count CV records convergence and support without verbose output", {
  data(countExample)
  fit <- catecvcount(data = countExample[1:500, ],
    cate.model = y ~ age + previous_treatment + offset(log(years)),
    ps.model = trt ~ age + previous_treatment,
    score.method = "contrastReg", initial.predictor.method = "poisson",
    B = 1, Kfold = 2, cv.n = 1, abc = FALSE, verbose = 0, seed = 45)
  expect_length(fit$ate.contrastReg$converge.contrastReg.cv, 1)
  expect_length(fit$ate.contrastReg$converge.repetitions.cv[[1]], 1)
  expect_true(nrow(fit$subgroup.support) > 0)
  expect_true(all(c("training", "validation") %in% fit$subgroup.support$split))
  expect_s3_class(diagnostics(fit), "data.frame")
})
