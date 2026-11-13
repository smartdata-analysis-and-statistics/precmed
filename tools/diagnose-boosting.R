# Run with the GBM3 branch installed and both gbm and gbm3 available:
# Rscript tools/diagnose-boosting.R diagnosis.csv
# Compare fixed tree counts to separate fitting differences from CV selection.
# The positive-count variant is a diagnostic perturbation, not an analysis of
# the original outcome. bag.fraction=1 disables random row subsampling.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript tools/diagnose-boosting.R OUTPUT.csv")
fit.boosting <- getFromNamespace("fit_boosting", "precmed")
if (!requireNamespace("gbm", quietly = TRUE) || !requireNamespace("gbm3", quietly = TRUE)) {
  stop("Install both gbm and gbm3 before running this diagnostic.")
}
data.env <- new.env()
utils::data("countExample", package = "precmed", envir = data.env)
set.seed(2026)
dat <- data.env$countExample[sort(sample.int(nrow(data.env$countExample), 1000)), ]
dat$trt <- as.integer(dat$trt == "drug1")
rows <- list()
for (arm in 0:1) for (outcome in c("original", "positive")) {
  d <- dat[dat$trt == arm, ]
  if (outcome == "positive") d$y <- d$y + 1
  for (bag in c(1, 0.5)) for (cv in c(0, 5)) {
    common <- list(formula = y ~ age + female + previous_cost + offset(log(years)),
      data = d, distribution = "poisson", n.trees = 100,
      interaction.depth = 2, shrinkage = 0.1,
      bag.fraction = bag, cv.folds = cv, n.cores = 1)
    set.seed(11)
    invisible(capture.output(old <- do.call(gbm::gbm, common)))
    set.seed(11)
    new <- do.call(fit.boosting, common)
    stopifnot(old$shrinkage == new$params$shrinkage,
              old$nTrain == new$params$num_train)
    for (trees in c(1, 10, 100)) {
      a <- suppressWarnings(predict(old, d, n.trees = trees))
      b <- suppressWarnings(predict(new, d, n.trees = trees))
      rows[[length(rows) + 1L]] <- data.frame(
        arm = arm, outcome = outcome, bag.fraction = bag, cv.folds = cv,
        trees = trees, correlation = cor(a, b),
        rmse = sqrt(mean((a - b)^2)), max.absolute.difference = max(abs(a - b)),
        gbm.version = as.character(utils::packageVersion("gbm")),
        gbm3.version = as.character(utils::packageVersion("gbm3")))
    }
  }
}
results <- do.call(rbind, rows)
utils::write.csv(results, args[1], row.names = FALSE)
print(results, row.names = FALSE)
