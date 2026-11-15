test_that("validation plots build without ggplot2 deprecation warnings", {
  for (response in c("count", "survival")) {
    fit <- structure(list(
      response = response, cv.n = 2, higher.y = FALSE,
      props = list(prop.onlyhigh = c(0.5, 0.75, 1)),
      formulas = list(ps.model = trt ~ 1),
      overall.ate.train = c(1, 1), overall.ate.valid = c(1, 1),
      ate.poisson = list(
        ate.est.train.high.cv = matrix(c(0.6, 0.8, 1), 3, 2),
        ate.est.valid.high.cv = matrix(c(0.7, 0.9, 1), 3, 2)
      )
    ), class = "precmed")

    for (cv.i in list(NULL, 1)) {
      expect_warning(p <- plot(fit, cv.i = cv.i, show.abc = FALSE), NA)
      expect_warning(built <- ggplot2::ggplot_build(p), NA)
      expect_equal(nrow(built$data[[1]]), 6L)
      expect_equal(nrow(built$data[[2]]), 2L)
    }
  }
})

test_that("subgroup box plots select columns without tidyselect warnings", {
  fit <- structure(list(
    ate.poisson = list(ate.est.valid.group.cv = matrix(c(0.7, 1.3), 2, 2)),
    response = "count", higher.y = FALSE,
    props = list(prop.multi = c(0, 0.5, 1)),
    formulas = list(ps.model = trt ~ 1),
    overall.ate.valid = c(1, 1)
  ), class = "precmed")

  expect_warning(p <- boxplot(fit), NA)
  expect_warning(built <- ggplot2::ggplot_build(p), NA)
  expect_equal(nrow(built$data[[1]]), 2L)
})

test_that("hazard plots default to cumulative ratios and retain explicit Cox plots", {
  m <- function(v) matrix(rep(v, 2), ncol = 2)
  fit <- structure(list(
    response = "survival", cv.n = 2, higher.y = TRUE, abc = FALSE,
    props = list(prop.onlyhigh = c(0.5, 1), prop.multi = c(0, 0.5, 1)),
    formulas = list(ps.model = trt ~ 1),
    overall.ate.train = c(1, 1), overall.ate.valid = c(1, 1),
    overall.hr.train = c(3, 3), overall.hr.valid = c(3, 3),
    overall.chr.train = c(2, 2), overall.chr.valid = c(2, 2),
    ate.poisson = list(ate.est.train.high.cv = m(c(1, 1)),
                       ate.est.valid.high.cv = m(c(1, 1)), ate.est.valid.group.cv = m(c(1, 1))),
    hr.poisson = list(hr.est.train.high.cv = m(c(3, 3)), hr.est.valid.high.cv = m(c(3, 3)),
                      hr.est.valid.group.cv = m(c(3, 3))),
    chr.poisson = list(chr.est.train.high.cv = m(c(2, 2)), chr.est.valid.high.cv = m(c(2, 2)),
                       chr.est.valid.group.cv = m(c(2, 2)))
  ), class = "precmed")
  fit <- structure(fit[c("ate.poisson", setdiff(names(fit), "ate.poisson"))], class = "precmed")
  expect_warning(cumulative <- plot(fit, plot.hr = TRUE, show.abc = FALSE), "deprecated")
  expect_warning(conditional <- plot(fit, plot.hr = TRUE, hazard.type = "conditional", show.abc = FALSE), "deprecated")
  expect_true(all(ggplot2::ggplot_build(cumulative)$data[[1]]$y == 2))
  expect_true(all(ggplot2::ggplot_build(conditional)$data[[1]]$y == 3))
  expect_match(cumulative$labels$y, "Cumulative hazard ratio")
  expect_match(conditional$labels$y, "Conditional Cox hazard ratio")
  expect_match(suppressWarnings(boxplot(fit, plot.hr = TRUE))$labels$y, "Cumulative hazard ratio")
  expect_match(suppressWarnings(boxplot(fit, plot.hr = TRUE, hazard.type = "conditional"))$labels$y, "Conditional Cox")
  expect_true(all(plot(fit, measure = "cumulative_hazard_ratio", show.abc = FALSE)$data$ate.cv == 2))
  expect_true(all(plot(fit, measure = "conditional_hazard_ratio", show.abc = FALSE)$data$ate.cv == 3))
  expect_match(boxplot(fit, measure = "cumulative_hazard_ratio")$labels$y, "Cumulative hazard ratio")
  expect_true(all(plot(fit, measure = "conditional_hr")$data$ate.cv == 3))
  expect_true(all(plot(fit, measure = "cumulative_hr")$data$ate.cv == 2))
  expect_match(boxplot(fit, measure = "conditional_hr")$labels$y, "Conditional Cox")
  both <- plot(fit, measure = c("rmtl", "cumulative_hr"))
  built <- ggplot2::ggplot_build(both)
  expect_equal(nrow(built$layout$layout), 4)
  expect_equal(levels(both$data$measure), c("RMTL ratio", "Cumulative hazard ratio"))
  expect_equal(levels(plot(fit, measure = c("cph", "rmtl"))$data$measure), c("Cumulative hazard ratio", "RMTL ratio"))
  expect_error(plot(fit, measure = c("chr", "cph")), "duplicates")
  expect_error(plot(fit, measure = character()), "measure")
  expect_true(all(both$data$ate.cv[both$data$measure == "RMTL ratio"] == 1))
  expect_true(all(both$data$ate.cv[both$data$measure == "Cumulative hazard ratio"] == 2))
  expect_equal(nrow(ggplot2::ggplot_build(plot(fit, measure = c("rmtl", "cph"), valid.only = TRUE))$layout$layout), 2)
  expect_error(plot(fit, measure = "rmtl", plot.hr = TRUE), "deprecated")
  fit$chr.poisson$chr.est.valid.high.cv[1, 1] <- NA_real_
  expect_warning(plot(fit, measure = "chr", show.abc = FALSE), "undefined")
  fit$chr.poisson$chr.est.valid.group.cv[1, 1] <- NA_real_
  expect_warning(boxplot(fit, measure = "chr"), "undefined")
  fit$overall.chr.valid <- NULL
  expect_error(plot(fit, measure = "chr", show.abc = FALSE), "refit")
  expect_s3_class(plot(fit, measure = "conditional_hazard_ratio", show.abc = FALSE), "ggplot")
})
test_that("ATE survival plot labels treatment panels and supports old fits", {
  fit <- structure(list(response = "survival", n.boot = 2,
    treatment.labels = c("control", "new drug"),
    rmst0 = data.frame(estimate = 3), rmst1 = data.frame(estimate = 4),
    log.rmtl.ratio = data.frame(estimate = 0.1),
    log.cumulative.hazard.ratio = data.frame(estimate = 0.2),
    trt.boot = data.frame(rmst1 = c(3, 4), rmst0 = c(2, 3),
                          log.rmtl.ratio = c(0.1, 0.2),
                          log.cumulative.hazard.ratio = c(0.2, 0.3))), class = "atefit")
  labels <- plot(fit)$facet$params$labeller
  keys <- data.frame(key = c("rmst1", "rmst0", "log.rmtl.ratio", "log.cumulative.hazard.ratio"))
  expect_equal(unname(unlist(labels(keys))), c("RMST (new drug)", "RMST (control)",
    "Log RMTL ratio", "Log cumulative hazard ratio"))
  fit$treatment.labels <- NULL
  expect_equal(unname(unlist(plot(fit)$facet$params$labeller(keys)))[1:2],
               c("RMST (treatment 1)", "RMST (treatment 0)"))
  expect_s3_class(ggplot2::ggplot_build(plot(fit)), "ggplot_built")
})
