test_that("ATE and CATE entry points reject invalid response values", {
  for (response in list(NULL, character(), NA_character_, 1, c("count", "survival"), "other")) {
    expect_error(atefit(response), "response must be one of")
    expect_error(catefit(response), "response must be one of")
    expect_error(catecv(response), "response must be one of")
  }
})

test_that("score methods match the outcome and are nonempty character values", {
  valid <- list(count = c("boosting", "poisson", "twoReg", "contrastReg", "negBin"),
                survival = c("boosting", "poisson", "twoReg", "contrastReg", "randomForest"),
                continuous = c("boosting", "gaussian", "twoReg", "contrastReg", "randomForest", "gam"))
  for (response in names(valid)) {
    expect_silent(precmed:::check_score_method(response, valid[[response]]))
    for (method in list(NULL, character(), NA_character_, 1, c(valid[[response]][1], NA_character_), "invalid")) {
      expect_error(precmed:::check_score_method(response, method), "score.method")
    }
  }
  for (response in names(valid)) {
    invalid <- switch(response, count = "randomForest", survival = "negBin", continuous = "poisson")
    expect_error(catefit(response, score.method = invalid), paste0("score.method for ", response))
    expect_error(catecv(response, score.method = invalid), paste0("score.method for ", response))
  }
  expect_error(catefitmean(data = data.frame(), score.method = "poisson"), "score.method for continuous")
  expect_error(catecvmean(data = data.frame(), score.method = character()), "score.method for continuous")
  expect_error(catefitcount(data = data.frame(), score.method = "randomForest"), "score.method for count")
  expect_error(catecvsurv(data = data.frame(), score.method = "negBin"), "score.method for survival")
})
