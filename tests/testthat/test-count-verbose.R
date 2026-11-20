test_that("verbose count fitting does not report runtime", {
  output <- capture.output(
    fit <- catefitcount(
      data = countExample,
      outcome_model = y ~ age + offset(log(years)),
      ps_model = trt ~ age,
      seed = 999, verbose = 1
    )
  )
  expect_s3_class(fit, "catefit")
  expect_false(any(grepl("Total runtime", output, fixed = TRUE)))
})
