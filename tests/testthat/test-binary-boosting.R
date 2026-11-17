test_that("binary boosting predicts treatment-specific response differences", {
  d <- simulate_binary(n = 500, seed = 91)
  args <- list(response = "binary", data = d, benefit_method = c("logistic", "boosting"),
               outcome_model = y ~ age + female + previous_treatment, ps.model = trt ~ age,
               prop.cutoff = 1, n.trees.boosting = 20, plot.gbmperf = FALSE,
               n.minobsinnode = 5, n.cores = 1, seed = 50)
  fit <- do.call(catefit, args)
  expect_equal(length(fit$score.boosting), nrow(d))
  expect_true(all(is.finite(fit$score.boosting) & abs(fit$score.boosting) <= 1))
  expect_equal(unname(fit$ate.boosting), unname(fit$ate.logistic))
  expect_s3_class(fit$fit$boosting$fit0.boosting, "GBMFit")
  expect_s3_class(fit$fit$boosting$fit1.boosting, "GBMFit")
  expect_true(all(c(fit$fit$boosting$best0.iter, fit$fit$boosting$best1.iter) %in% 1:20))
  x <- model.matrix(y ~ age + female + previous_treatment, d)[, -1, drop = FALSE]
  newdata <- as.data.frame(x)
  names(newdata) <- paste0("x", seq_len(ncol(x)))
  p0 <- predict(fit$fit$boosting$fit0.boosting, newdata, n.trees = fit$fit$boosting$best0.iter, type = "response")
  p1 <- predict(fit$fit$boosting$fit1.boosting, newdata, n.trees = fit$fit$boosting$best1.iter, type = "response")
  expect_equal(fit$score.boosting, as.numeric(p1 - p0))
  expect_equal(do.call(catefit, args)$score.boosting, fit$score.boosting)
  expect_output(print(fit), "boosting")
})

test_that("binary boosting validates unsupported and small-sample inputs", {
  d <- simulate_binary(n = 200, seed = 12)
  common <- list(response = "binary", data = d, benefit_method = "boosting", ps.model = trt ~ 1,
                 n.trees.boosting = 10, plot.gbmperf = FALSE, n.cores = 1)
  expect_error(do.call(catefit, c(common, list(outcome_model = y ~ 1))), "at least one outcome-model covariate")
  d$y[d$trt == "drug0"] <- 0
  common$data <- d
  expect_error(do.call(catefit, c(common, list(outcome_model = y ~ age))), "both outcome values")
  expect_warning(cv <- catecv("binary", d, "boosting", y ~ age, trt ~ 1,
                              cv.n = 1, abc = FALSE, seed = 1), "score fitting failed")
  expect_true(all(is.na(cv$ate.boosting$ate.est.valid.high.cv)))
})
