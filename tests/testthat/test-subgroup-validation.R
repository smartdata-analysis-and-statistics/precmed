test_that("valid subgroup cutoffs do not trigger vector-recycling warnings", {
  check <- function(cutoffs) precmed:::arg.checks.common(
    fun = "crossv", ps.method = "glm", minPS = 0.01, maxPS = 0.99,
    higher.y = FALSE, abc = TRUE, prop.cutoff = c(0.5, 1),
    prop.multi = cutoffs, B = 1, Kfold = 2, plot.gbmperf = FALSE,
    tree.depth = 2, n.trees.boosting = 20, error.maxNR = 1e-3,
    max.iterNR = 150, tune = c(0.5, 2), train.prop = 0.75,
    cv.n = 2, error.max = 0.1, max.iter = 100)
  expect_warning(check(c(0, 0.5, 1)), NA)
  expect_warning(check(c(0, 0.25, 0.5, 0.75, 1)), NA)
  expect_error(check(c(0, 1)), "at least 2 subgroups")
})
