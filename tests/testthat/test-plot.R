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
