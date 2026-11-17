test_that("boosting plots are opt-in across fitting and validation", {
  functions <- c("catefit", "catefitbinary", "catefitcount", "catefitmean",
                 "catefitsurv", "intxcount", "intxmean", "intxsurv",
                 "catecv", "catecvbinary", "catecvcount", "catecvmean",
                 "catecvsurv")
  for (name in functions) {
    expect_identical(formals(get(name, envir = asNamespace("precmed")))$plot.gbmperf,
                     FALSE, info = name)
  }

  seen <- logical()
  original <- precmed:::best_boosting_iteration
  local_mocked_bindings(
    best_boosting_iteration = function(object, method, plot.it) {
      seen <<- c(seen, plot.it)
      original(object, method, plot.it = FALSE)
    }, .package = "precmed")
  d <- simulate_binary(n = 300, seed = 42)
  args <- list(response = "binary", data = d, benefit_method = "boosting",
               outcome_model = y ~ age, ps.model = trt ~ 1, prop.cutoff = 1,
               n.trees.boosting = 10, n.minobsinnode = 5, n.cores = 1, seed = 42)
  do.call(catefit, args)
  expect_identical(seen, c(FALSE, FALSE))
  seen <- logical()
  do.call(catefit, c(args, list(plot.gbmperf = TRUE)))
  expect_identical(seen, c(TRUE, TRUE))
})
