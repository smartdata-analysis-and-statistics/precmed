test_that("catefit score plots show the requested methods and effect scales", {
  fit <- structure(list(response = "binary", score.method = c("twoReg", "logistic"),
                        score.twoReg = c(-.2, .1, .3), score.logistic = c(-.1, .2, .4)), class = "catefit")
  p <- plot(fit)
  expect_s3_class(p, "ggplot")
  expect_equal(levels(p$data$method), c("twoReg", "logistic"))
  expect_equal(nrow(ggplot2::ggplot_build(p)$layout$layout), 2L)
  expect_match(p$labels$x, "risk difference")
  expect_match(p$labels$x, "treatment 1 - treatment 0", fixed = TRUE)
  expect_equal(ggplot2::ggplot_build(p)$data[[2]]$xintercept, c(0, 0))
  selected <- plot(fit, methods = "logistic", type = "density", null.value = NULL)
  expect_equal(levels(selected$data$method), "logistic")
  expect_length(selected$layers, 1L)
  expect_s3_class(ggplot2::ggplot_build(selected), "ggplot_built")
  for (response in c("count", "survival", "continuous")) {
    fit$response <- response
    expect_match(plot(fit)$labels$x, switch(response, count = "log rate ratio", survival = "log RMTL ratio", continuous = "mean difference"))
    fit$treatment.labels <- c("control", "new drug")
    expect_match(plot(fit)$labels$x, if (response == "continuous") "new drug - control" else "new drug vs control", fixed = TRUE)
  }
})

test_that("catefit retains treatment labels for score plots", {
  d <- simulate_binary(n = 250, seed = 93)
  binary <- catefit("binary", d, "logistic", y ~ age, trt ~ age, prop.cutoff = 1)
  expect_identical(binary$treatment.labels, c("drug0", "drug1"))
  expect_match(plot(binary)$labels$x, "drug1 - drug0", fixed = TRUE)
  d$trt <- d$trt == "drug1"
  logical_fit <- catefit("binary", d, "logistic", y ~ age, trt ~ age, prop.cutoff = 1)
  expect_match(plot(logical_fit)$labels$x, "TRUE - FALSE", fixed = TRUE)

  count_data <- simulate_count(n = 300, seed = 94)
  count <- suppressWarnings(catefitcount(count_data, "poisson", y ~ age + offset(log(years)), trt ~ age, prop.cutoff = 1))
  expect_identical(count$treatment.labels, c("drug0", "drug1"))
  expect_match(plot(count)$labels$x, "drug1 vs drug0", fixed = TRUE)

  surv_data <- simulate_survival(n = 350, seed = 95)
  survival <- suppressWarnings(catefitsurv(surv_data, "poisson", survival::Surv(y, d) ~ age, trt ~ age, tau0 = 300, prop.cutoff = c(.5, 1)))
  expect_identical(survival$treatment.labels, c("drug0", "drug1"))
  expect_match(plot(survival)$labels$x, "drug1 vs drug0", fixed = TRUE)
})

test_that("catefit plots handle missing values, constant scores, and bad arguments", {
  fit <- structure(list(response = "binary", score.method = "twoReg", score.twoReg = c(1, 1, NA_real_)), class = "catefit")
  expect_warning(p <- plot(fit, null.value = NA), "Non-finite")
  expect_equal(nrow(p$data), 2)
  expect_length(p$layers, 1L)
  expect_error(suppressWarnings(plot(fit, type = "density")), "two distinct")
  expect_error(plot(fit, methods = "other"), "available scoring methods")
  expect_error(plot(fit, methods = c("twoReg", "twoReg")), "distinct")
  expect_error(plot(fit, bins = 0), "bins")
  expect_error(plot(fit, null.value = Inf), "null.value")
  expect_error(plot(structure(list(response = "count"), class = "catefit")), "no scores")
  fit$score.twoReg <- c(NA_real_, Inf)
  expect_error(suppressWarnings(plot(fit)), "at least one finite")
})
