test_that("binary simulator thresholds the count DGM at a common horizon", {
  d <- simulate_binary(n = 200, seed = 17, followup = 2, threshold = 1, dispersion = 2)
  counts <- simulate_count(n = 200, seed = 17, followup_range = c(2, 2), dispersion = 2)
  expect_identical(d$relapse_count, counts$y)
  expect_identical(d$y, as.integer(counts$y > 1))
  expect_true(all(d$years == 2))
  expect_identical(d, simulate_binary(n = 200, seed = 17, followup = 2, threshold = 1, dispersion = 2))
  expect_identical(simulate_binary(n = 200, seed = 17)$y,
                   as.integer(simulate_count(n = 200, seed = 17, followup_range = c(1, 1))$y > 0))
  expect_error(simulate_binary(followup = 0), "followup")
  expect_error(simulate_binary(threshold = -1), "threshold")
  expect_error(simulate_binary(threshold = 1.5), "threshold")
})
