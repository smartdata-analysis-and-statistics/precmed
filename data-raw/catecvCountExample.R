# Regenerate the pre-fitted count validation example from the repository root.
pkgload::load_all(".")
catecvCountExample <- catecv(
  response = "count", data = countExample, score.method = "poisson",
  cate.model = y ~ age + female + previous_treatment + previous_cost +
    previous_number_relapses + offset(log(years)),
  ps.model = trt ~ age + previous_treatment,
  higher.y = FALSE, cv.n = 5, seed = 999, verbose = 0
)
save(catecvCountExample, file = "data/catecvCountExample.rda",
     compress = "xz", version = 2)
