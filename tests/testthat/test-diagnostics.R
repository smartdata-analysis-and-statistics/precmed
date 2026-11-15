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
  expect_error(diagnostics(structure(list(response = "count"), class = "precmed")), "survival")
})
