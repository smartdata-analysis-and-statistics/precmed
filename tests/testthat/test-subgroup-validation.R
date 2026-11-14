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

test_that("count subgroup estimation reports empty and single-treatment groups", {
  n <- 80L
  x <- matrix(seq_len(n) / n, ncol = 1)
  trt <- rep(0:1, each = n / 2)
  y <- rep(1:4, length.out = n)
  base <- list(y = y, x.cate = x, x.ps = cbind(1, x), time = rep(0, n),
               trt = trt, ps.method = "glm", minPS = 0.01, maxPS = 0.99,
               ps.control = list(), label = "test subgroup")
  call <- function(selected) do.call(precmed:::count_subgroup_log_effect,
                                      c(list(selected = selected), base))
  expect_warning(empty <- call(rep(FALSE, n)), "treatment 0: 0, treatment 1: 0")
  expect_true(is.na(empty))
  expect_warning(one_arm <- call(trt == 0), "treatment 1: 0")
  expect_true(is.na(one_arm))
  expect_true(is.finite(call(rep(TRUE, n))))
  expect_warning(tied <- precmed:::estcount.bilevel.subgroups(
    y = y, x.cate = x, x.ps = cbind(1, x), time = rep(0, n), trt = trt,
    score = rep(0, n), higher.y = TRUE, prop = 0.5, onlyhigh = FALSE,
    context = "test iteration 1"), "test iteration 1, low subgroup")
  expect_true(is.finite(tied$ate.est.high))
  expect_true(is.na(tied$ate.est.low))
})
