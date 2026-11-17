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
  expect_match(boxplot(fit, measure = "conditional.hr")$labels$y, "Conditional Cox")
  expect_error(plot(fit, measure = "risk_difference"), "refit")
  fit$overall.rd.train <- fit$overall.rd.valid <- c(-.1, -.1)
  fit$rd.poisson <- list(rd.est.train.high.cv = m(c(-.2, -.1)),
                         rd.est.valid.high.cv = m(c(-.2, -.1)))
  rd <- plot(fit, measure = "risk.difference")
  expect_equal(rd$data$ate.cv, plot(fit, measure = "risk_difference")$data$ate.cv)
  expect_true(all(rd$data$ate.cv < 0))
  arr <- plot(fit, measure = "absolute_risk_reduction")
  expect_equal(arr$data$ate.cv, -rd$data$ate.cv)
  expect_equal(arr$layers[[2]]$data$yintercept, -rd$layers[[2]]$data$yintercept)
  combined <- plot(fit, measure = c("cumulative.hr", "risk.difference"))
  expect_equal(nrow(ggplot2::ggplot_build(combined)$layout$layout), 4L)
  both <- plot(fit, measure = c("rmtl", "cumulative_hr"))
  method <- plot(fit, measure = c("rmtl", "cumulative_hr"), layout = "method")
  expect_equal(nrow(ggplot2::ggplot_build(method)$layout$layout), 2)
  expect_match(method$labels$caption, "Dotted")
  expect_s3_class(ggplot2::ggplot_build(plot(fit, measure = c("rmtl", "cumulative_hr"), layout = "method", show.overall = FALSE)), "ggplot_built")
  fit$ate.contrastReg <- fit$ate.poisson
  fit$chr.contrastReg <- fit$chr.poisson
  by.method <- ggplot2::ggplot_build(plot(fit, measure = c("rmtl", "cumulative_hr"), layout = "method"))
  expect_equal(nrow(by.method$layout$layout), 4L)
  expect_equal(nrow(by.method$data[[2]]), 8L)
  expect_equal(by.method$layout$panel_params[[1]]$y.range, by.method$layout$panel_params[[2]]$y.range)
  expect_s3_class(ggplot2::ggplot_build(plot(fit, layout = "method", valid.only = TRUE, show.abc = FALSE)), "ggplot_built")
  fit$ate.contrastReg <- NULL
  fit$chr.contrastReg <- NULL
  expect_error(plot(fit, layout = "unknown"), "layout")
  built <- ggplot2::ggplot_build(both)
  expect_equal(nrow(built$layout$layout), 4)
  hidden <- plot(fit, measure = c("rmtl", "cumulative_hr"), show.overall = FALSE)
  expect_null(hidden$labels$caption)
  expect_false(any(vapply(hidden$layers, function(layer) inherits(layer$geom, "GeomHline"), logical(1))))
  expect_equal(nrow(ggplot2::ggplot_build(hidden)$layout$layout), 4)
  expect_error(plot(fit, show.overall = NA), "show.overall")
  expect_equal(built$layout$panel_params[[1]]$y.range, built$layout$panel_params[[2]]$y.range)
  expect_equal(built$layout$panel_params[[3]]$y.range, built$layout$panel_params[[4]]$y.range)
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

test_that("ATE bootstrap null lines use the log-ratio panels", {
  count <- structure(list(response = "count", n.boot = 2,
    log.rate.ratio = data.frame(estimate = 0.1), trt.boot = c(0.1, 0.2)), class = "atefit")
  expect_equal(tail(ggplot2::ggplot_build(plot(count))$data, 1)[[1]]$xintercept, 0)
  expect_equal(tail(ggplot2::ggplot_build(plot(count, null.value = 0.5))$data, 1)[[1]]$xintercept, 0.5)
  expect_length(plot(count, null.value = NULL)$layers, 2)
  expect_length(plot(count, null.value = NA)$layers, 2)
  expect_error(plot(count, metrics = "rate0"), "available metric")
  expect_error(plot(count, metrics = character()), "available metric")
  expect_error(plot(count, metrics = c("log.rate.ratio", "log.rate.ratio")), "available metric")
  expect_error(plot(count, null.value = Inf), "null.value")
  expect_error(plot(count, null.value = c(0, 1)), "null.value")

  surv <- structure(list(response = "survival", n.boot = 2,
    rmst0 = data.frame(estimate = 3), rmst1 = data.frame(estimate = 4),
    log.rmtl.ratio = data.frame(estimate = 0.1),
    log.cumulative.hazard.ratio = data.frame(estimate = 0.2),
    trt.boot = data.frame(rmst1 = c(3, 4), rmst0 = c(2, 3),
      log.rmtl.ratio = c(0.1, 0.2), log.cumulative.hazard.ratio = c(0.2, 0.3))), class = "atefit")
  built <- ggplot2::ggplot_build(plot(surv))
  layer <- tail(built$data, 1)[[1]]
  expect_setequal(as.character(built$layout$layout$key[as.integer(layer$PANEL)]),
                  c("log.rmtl.ratio", "log.cumulative.hazard.ratio"))
  expect_equal(layer$xintercept, c(0, 0))
  expect_length(plot(surv, null.value = NULL)$layers, 2)
  selected <- ggplot2::ggplot_build(plot(surv, metrics = c("rmst0", "rmst1")))
  expect_equal(as.character(selected$layout$layout$key), c("rmst0", "rmst1"))
  expect_length(selected$data, 2)
  surv$trt.boot <- as.matrix(surv$trt.boot)
  selected <- ggplot2::ggplot_build(plot(surv, metrics = c("log.cumulative.hazard.ratio", "rmst0")))
  expect_equal(as.character(selected$layout$layout$key), c("log.cumulative.hazard.ratio", "rmst0"))
  expect_equal(as.integer(tail(selected$data, 1)[[1]]$PANEL), 1L)
})


test_that("survival ATE facets group related quantities and honor metric order", {
  keys <- c("rmst1", "rmst0", "log.rmtl.ratio", "risk0", "risk1",
            "risk.difference", "log.cumulative.hazard.ratio")
  fit <- structure(c(list(response = "survival", n.boot = 2),
    setNames(lapply(seq_along(keys), function(i) data.frame(estimate = i / 10)), keys),
    list(trt.boot = as.data.frame(setNames(lapply(seq_along(keys), function(i) c(i / 10, i / 10 + .01)), keys)))),
    class = "atefit")
  layout <- ggplot2::ggplot_build(plot(fit))$layout$layout
  expect_equal(as.character(layout$key), keys)
  expect_equal(layout$ROW, c(1L, 1L, 1L, 2L, 2L, 2L, 3L))
  selected <- ggplot2::ggplot_build(plot(fit, metrics = c("risk_difference", "rmst0")))$layout$layout
  expect_equal(as.character(selected$key), c("risk.difference", "rmst0"))
})

test_that("binary ATE plots select treatment risks and risk differences", {
  fit <- structure(list(response = "binary", n.boot = 3,
    treatment.labels = c("control", "new drug"),
    risk0 = data.frame(estimate = .4), risk1 = data.frame(estimate = .25),
    risk.difference = data.frame(estimate = -.15), trt.boot = c(-.1, -.2, -.15),
    risk.boot = data.frame(risk0 = c(.35, .45, .4), risk1 = c(.25, .25, .25))), class = "atefit")
  p <- plot(fit)
  built <- ggplot2::ggplot_build(p)
  expect_equal(as.character(built$layout$layout$key), c("risk0", "risk1", "risk.difference"))
  expect_equal(unname(unlist(p$facet$params$labeller(data.frame(key = c("risk0", "risk1", "risk.difference"))))),
               c("Risk (control)", "Risk (new drug)", "Risk difference"))
  expect_equal(built$data[[2]]$xintercept, c(.4, .25, -.15))
  expect_equal(as.integer(built$data[[3]]$PANEL), 3L)
  selected <- ggplot2::ggplot_build(plot(fit, metrics = c("risk.difference", "risk0")))
  expect_equal(as.character(selected$layout$layout$key), c("risk.difference", "risk0"))
  expect_equal(as.integer(selected$data[[3]]$PANEL), 1L)
  expect_length(plot(fit, metrics = "risk1")$layers, 2L)
  expect_length(plot(fit, null.value = NULL)$layers, 2L)
  expect_error(plot(fit, metrics = c("risk0", "risk0")), "distinct")
  expect_equal(plot(fit, metrics = "risk_difference")$data, plot(fit, metrics = "risk.difference")$data)
  expect_error(plot(fit, metrics = "p0"), "available metric")
  fit$treatment.labels <- NULL
  expect_equal(unname(unlist(plot(fit)$facet$params$labeller(data.frame(key = c("risk0", "risk1"))))),
               c("Risk (treatment 0)", "Risk (treatment 1)"))
  fit$risk.boot <- NULL
  expect_s3_class(plot(fit), "ggplot")
  expect_error(plot(fit, metrics = "risk0"), "available metric")
})

test_that("count ATE plots show treatment-specific bootstrap rates", {
  fit <- structure(list(response = "count", n.boot = 2,
    treatment.labels = c("control", "new drug"),
    rate0 = data.frame(estimate = 2), rate1 = data.frame(estimate = 3),
    log.rate.ratio = data.frame(estimate = log(1.5)), trt.boot = c(0.3, 0.4),
    rate.boot = data.frame(rate0 = c(1, 2), rate1 = c(3, 4))), class = "atefit")
  p <- plot(fit)
  built <- ggplot2::ggplot_build(p)
  expect_equal(as.character(built$layout$layout$key), c("rate0", "rate1", "log.rate.ratio"))
  expect_equal(unname(unlist(p$facet$params$labeller(data.frame(key = c("rate0", "rate1", "log.rate.ratio"))))),
               c("Rate (control)", "Rate (new drug)", "Log rate ratio"))
  expect_equal(as.integer(tail(built$data, 1)[[1]]$PANEL), 3L)
  expect_equal(built$data[[2]]$xintercept, c(2, 3, log(1.5)))
  selected <- ggplot2::ggplot_build(plot(fit, metrics = c("log.rate.ratio", "rate0")))
  expect_equal(as.character(selected$layout$layout$key), c("log.rate.ratio", "rate0"))
  expect_equal(as.integer(tail(selected$data, 1)[[1]]$PANEL), 1L)
  expect_length(plot(fit, metrics = "rate1")$layers, 2)
  expect_error(plot(fit, metrics = "rmst0"), "available metric")
  fit$treatment.labels <- NULL
  expect_equal(unname(unlist(plot(fit)$facet$params$labeller(data.frame(key = c("rate0", "rate1"))))),
               c("Rate (treatment 0)", "Rate (treatment 1)"))
})

test_that("count ATE retains bootstrap rates and original treatment labels", {
  data(countExample, package = "precmed")
  fit <- atefit(response = "count", data = countExample,
                cate.model = y ~ age + offset(log(years)), ps.model = trt ~ age,
                n.boot = 3, seed = 123, verbose = 0)
  expect_equal(fit$treatment.labels, c("drug0", "drug1"))
  expect_equal(nrow(fit$rate.boot), 3L)
  expect_equal(log(fit$rate.boot$rate1 / fit$rate.boot$rate0), fit$trt.boot,
               tolerance = 1e-10)
})

test_that("validation and bootstrap plots use the same default theme", {
  data(catecvCountExample, package = "precmed")
  validation <- plot(catecvCountExample, show.abc = FALSE)
  subgroup <- boxplot(catecvCountExample)
  expect_equal(validation$theme$panel.grid.major.y$colour, "#E5E9ED")
  expect_equal(subgroup$theme$panel.grid.major.y$colour, "#E5E9ED")
  expect_match(validation$labels$caption, "mean overall treatment effect")
  expect_match(plot(catecvCountExample, combine = "median", show.abc = FALSE)$labels$caption,
               "median overall treatment effect")
  expect_match(plot(catecvCountExample, cv.i = 1, show.abc = FALSE)$labels$caption,
               "overall treatment effect")
  expect_null(plot(catecvCountExample, show.overall = FALSE)$labels$caption)
  expect_match(subgroup$labels$caption, "validation set")
  expect_equal(subgroup$layers[[1]]$aes_params$fill, "#7293AD")
  subgroup.hidden <- boxplot(catecvCountExample, show.overall = FALSE)
  expect_null(subgroup.hidden$labels$caption)
  expect_false(any(vapply(subgroup.hidden$layers, function(layer) inherits(layer$geom, "GeomHline"), logical(1))))
  expect_error(boxplot(catecvCountExample, show.overall = NA), "show.overall")
  expect_error(boxplot(catecvCountExample, ncol = 0), "ncol")
  expect_error(boxplot(catecvCountExample, ncol = 1.5), "ncol")
  expect_equal(subgroup$facet$params$ncol, 2)
  expect_equal(boxplot(catecvCountExample, ncol = 1)$facet$params$ncol, 1)
  custom <- plot(catecvCountExample, show.abc = FALSE, theme = ggplot2::theme_classic())
  expect_s3_class(custom$theme$panel.grid, "element_blank")
})
