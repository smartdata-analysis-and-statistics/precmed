make_atefit_fixture <- function(response = "count") {
  effect <- function(estimate, se) {
    data.frame(estimate = estimate, SE = se,
               CI.lower = estimate - qnorm(0.975) * se,
               CI.upper = estimate + qnorm(0.975) * se)
  }
  structure(list(response = response,
                 log.rate.ratio = effect(log(2), 0.2),
                 log.rmtl.ratio = effect(log(0.8), 0.1),
                 log.hazard.ratio = effect(log(0.7), 0.15),
                 rate0 = data.frame(estimate = 0.1)), class = "atefit")
}

test_that("count effects and intervals are extracted on the log scale", {
  fit <- make_atefit_fixture()
  expect_equal(coef(fit), c(log.rate.ratio = log(2)))
  ci <- confint(fit)
  expect_equal(dim(ci), c(1L, 2L))
  expect_equal(rownames(ci), "log.rate.ratio")
  expect_equal(colnames(ci), c("2.5 %", "97.5 %"))
  expect_equal(unname(ci[1, ]),
               unlist(fit$log.rate.ratio[c("CI.lower", "CI.upper")], use.names = FALSE))
  expect_equal(unname(exp(coef(fit))), 2)
  expect_equal(unname(exp(ci[1, ])), exp(log(2) + c(-1, 1) * qnorm(0.975) * 0.2))
  expect_equal(unname(confint(fit, level = 0.90)[1, ]),
               log(2) + c(-1, 1) * qnorm(0.95) * 0.2)
})

test_that("survival effects support selection by name and position", {
  fit <- make_atefit_fixture("survival")
  expect_equal(coef(fit), c(log.rmtl.ratio = log(0.8), log.hazard.ratio = log(0.7)))
  ci <- confint(fit)
  expect_equal(dim(ci), c(2L, 2L))
  for (term in names(coef(fit))) {
    expect_equal(unname(ci[term, ]),
                 unlist(fit[[term]][c("CI.lower", "CI.upper")], use.names = FALSE))
  }
  expect_equal(confint(fit, parm = 2), ci[2, , drop = FALSE])
  expect_equal(confint(fit, parm = "log.hazard.ratio"), ci[2, , drop = FALSE])
  expect_equal(confint(fit, parm = c(2, 1)), ci[c(2, 1), , drop = FALSE])
})

test_that("invalid confidence levels and effect selections are rejected", {
  fit <- make_atefit_fixture()
  for (level in list(0, 1, -0.1, 1.1, NA_real_, Inf, numeric(), c(0.9, 0.95), "95%")) {
    expect_error(confint(fit, level = level), "level must")
  }
  for (parm in list("rate0", 0, 2, 1.5, NA_real_, NA_character_, Inf, TRUE)) {
    expect_error(confint(fit, parm = parm), "parm must")
  }
})
