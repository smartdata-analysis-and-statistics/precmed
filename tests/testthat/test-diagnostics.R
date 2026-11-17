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

test_that("overlap summaries use raw scores and truncated treatment weights", {
  p <- c(.1, .4, .6, .9)
  attr(p, "raw") <- c(.01, .4, .6, .99)
  x <- precmed:::treatment_overlap(c(0, 0, 1, 1), p, .1, .9)
  expect_equal(x$ps.min, c(.01, .6))
  expect_equal(x$outside.other.range, c(2L, 2L))
  expect_equal(x$truncated.low, c(1L, 0L))
  expect_equal(x$truncated.high, c(0L, 1L))
  expect_equal(x$effective.n, rep((1/.9 + 1/.6)^2 / (1/.9^2 + 1/.6^2), 2))
})

test_that("ATE overlap is available separately from messages", {
  set.seed(123)
  d <- data.frame(y = rep(c(0, 1, 1, 0), 20), trt = rep(c(0, 0, 1, 1), 20), x = rnorm(80))
  fit <- atefitbinary(d, y ~ x, trt ~ x, n.boot = 4, seed = 7)
  overlap <- diagnostics(fit, type = "propensity", layout = "arms")
  expect_equal(overlap$n, c(40L, 40L))
  expect_true(all(overlap$effective.n > 0 & overlap$effective.n <= overlap$n))
  expect_equal(nrow(diagnostics(fit)), 0L)
  expect_error(diagnostics(structure(list(response = "binary"), class = "atefit"), "propensity"), "refit")
})

test_that("count ATE retains overlap from its treatment model", {
  data(countExample, package = "precmed")
  expect_message(fit <- atefitcount(countExample,
    y ~ age + female + offset(log(years)), trt ~ age,
    n.boot = 2, seed = 123), "recoded")
  overlap <- diagnostics(fit, "propensity", layout = "arms")
  expect_equal(sum(overlap$n), nrow(countExample))
  expect_true(all(is.finite(overlap$effective.n)))
})

test_that("metric layout preserves values and treatment labels", {
  p <- c(.1, .4, .6, .9)
  attr(p, "raw") <- c(.01, .4, .6, .99)
  raw <- precmed:::treatment_overlap(c(0, 0, 1, 1), p, .1, .9)
  fit <- structure(list(response = "binary", overlap.diagnostics = raw,
                        treatment.labels = c("usual care", "new drug")), class = "atefit")
  display <- diagnostics(fit, "propensity", layout = "metrics")
  expect_named(display, c("Diagnostic", "usual care", "new drug"))
  expect_equal(display[[2]], as.numeric(raw[1, -1]))
  expect_equal(display[[3]], as.numeric(raw[2, -1]))
  expect_equal(display$Diagnostic[1], "Number of patients")
  expect_equal(display$Diagnostic[5], "Outside other group's range (n)")
  expect_equal(diagnostics(fit, "propensity", "arms"), raw)
  fit$treatment.labels <- NULL
  expect_named(diagnostics(fit, "propensity", "metrics"),
               c("Diagnostic", "treatment 0", "treatment 1"))
  expect_error(diagnostics(fit, "propensity", "invalid"), "arg")
})


test_that("formatted overlap is compact and uses within-arm percentages", {
  p <- c(.1, .4, .6, .9)
  attr(p, "raw") <- c(.02, .4, .6, .99)
  raw <- precmed:::treatment_overlap(c(0, 0, 1, 1), p, .1, .9)
  fit <- structure(list(response = "binary", overlap.diagnostics = raw,
    treatment.labels = c("usual care", "new drug")), class = "atefit")
  formatted <- diagnostics(fit, "propensity")
  expect_identical(diagnostics(fit, "overlap"), formatted)
  expect_identical(diagnostics(fit, "overlap", "arms"), diagnostics(fit, "propensity", "arms"))
  expect_named(formatted, c("Diagnostic", "usual care", "new drug"))
  expect_equal(nrow(formatted), 6L)
  expect_equal(formatted[[2]][2], "0.21 (0.02-0.40)")
  expect_equal(formatted[[2]][3], "2 (100.0%)")
  expect_equal(formatted[[2]][4], "1 / 0")
  expect_type(formatted[[2]], "character")
  expect_type(diagnostics(fit, "propensity", "arms")$ps.min, "double")
  raw$effective.n[1] <- NA_real_
  fit$overlap.diagnostics <- raw
  expect_equal(diagnostics(fit, "propensity")[[2]][6], "NA")
})
