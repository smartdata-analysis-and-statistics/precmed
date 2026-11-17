test_that("outcome_method follows outcome_model in ATE interfaces", {
  for (name in c("atefit", "atefitbinary", "atefitcount", "atefitsurv")) {
    args <- names(formals(get(name, envir = asNamespace("precmed"))))
    expect_identical(match("outcome_method", args), match("outcome_model", args) + 1L, info = name)
  }
})

test_that("new and previous positional ATE calls agree", {
  b <- simulate_binary(n = 350, seed = 511)
  old <- atefit("binary", b, y ~ age, trt ~ age, "glm", n.boot = 3, seed = 512)
  new <- atefit("binary", b, y ~ age, "logistic", trt ~ age, "glm", n.boot = 3, seed = 512)
  expect_equal(old, new)
  old_string <- atefit("binary", b, y ~ age, "trt ~ age", "glm", n.boot = 3, seed = 512)
  expect_equal(old, old_string)
  old_named_method <- atefit("binary", b, y ~ age, trt ~ age, "glm", outcome_method = "logistic", n.boot = 3, seed = 512)
  expect_equal(old, old_named_method)
  old_direct <- atefitbinary(b, y ~ age, trt ~ age, "glm", n.boot = 3, seed = 512)
  new_direct <- atefitbinary(b, y ~ age, "logistic", trt ~ age, "glm", n.boot = 3, seed = 512)
  expect_equal(old_direct, new_direct)
  fully_positional <- atefitbinary(b, y ~ age, trt ~ age, "glm", .01, .99, 3, 512)
  expect_equal(old_direct, fully_positional)
  d <- simulate_count(n = 350, seed = 513)
  old_count <- atefitcount(d, y ~ age + offset(log(years)), trt ~ age, "glm", n.boot = 3, seed = 514)
  new_count <- atefitcount(d, y ~ age + offset(log(years)), "poisson", trt ~ age, "glm", n.boot = 3, seed = 514)
  expect_equal(old_count, new_count)
  s <- simulate_survival(n = 350, seed = 515)
  tau <- median(s$y) * 1.1
  old_surv <- suppressMessages(atefitsurv(s, survival::Surv(y, d) ~ age, trt ~ age, "glm", tau0 = tau, n.boot = 3, seed = 516))
  new_surv <- suppressMessages(atefitsurv(s, survival::Surv(y, d) ~ age, "coxph", trt ~ age, "glm", tau0 = tau, n.boot = 3, seed = 516))
  expect_equal(old_surv, new_surv)
})

test_that("an explicitly named invalid outcome method is not treated as a treatment formula", {
  d <- simulate_binary(n = 100, seed = 517)
  expect_error(atefitbinary(d, y ~ age, ps_model = trt ~ age, outcome_method = y ~ age), "outcome_method")
})
