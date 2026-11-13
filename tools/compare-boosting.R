# Install main and codex/gbm3-transition into separate R libraries first, then:
# Rscript tools/compare-boosting.R /path/to/main-library /path/to/gbm3-library comparison.rds
# Dependency libraries can be supplied through the usual R_LIBS_USER setting.
# Each version runs in a fresh process with identical balanced outer splits.
# Inner model-fitting randomness remains specific to each backend.
# Neither checkout is modified.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript tools/compare-boosting.R MAIN_LIBRARY GBM3_LIBRARY OUTPUT.rds")
}
libraries <- vapply(args[1:2], normalizePath, character(1), mustWork = TRUE)
worker <- tempfile(fileext = ".R")
writeLines(c(
  'args <- commandArgs(trailingOnly = TRUE)',
  '.libPaths(c(args[1], .libPaths()))',
  'library(precmed)',
  'data(countExample, package = "precmed")',
  'set.seed(2026)',
  'dat <- countExample[sort(sample.int(nrow(countExample), 1000)), ]',
  'dat$trt <- as.integer(dat$trt == "drug1")',
  'results <- lapply(c(11L, 29L, 51L), function(seed) {',
  '  common <- list(data = dat, score.method = "boosting",',
  '    cate.model = y ~ age + female + previous_cost + offset(log(years)),',
  '    ps.model = trt ~ age + previous_cost, higher.y = FALSE,',
  '    prop.cutoff = c(0.5, 0.75, 1), seed = seed,',
  '    n.trees.boosting = 100, plot.gbmperf = FALSE, n.cores = 1)',
  '  elapsed <- system.time(fit <- do.call(catefitcount, common))[["elapsed"]]',
  '  # Reset the splitter seed independently of RNG consumed by model fitting.',
  '  splitter <- getFromNamespace("balance.split", "precmed")',
  '  split.count <- 0L',
  '  validation.rows <- list()',
  '  fixed.cv <- catecvcount',
  '  environment(fixed.cv) <- list2env(list(balance.split = function(...) {',
  '    split.count <<- split.count + 1L',
  '    set.seed(100000L + seed + split.count)',
  '    split <- splitter(...)',
  '    validation.rows[[split.count]] <<- rownames(split$x.cate.valid)',
  '    split',
  '  }), parent = environment(fixed.cv))',
  '  cv <- do.call(fixed.cv, c(common, list(cv.n = 3, max.iter = 100, error.max = 0.2)))',
  '  stopifnot(length(validation.rows) == 3L, all(lengths(validation.rows) == 250L))',
  '  list(seed = seed, elapsed = elapsed, scores = fit$score.boosting,',
  '    ate = fit$ate.boosting, abc = abc(cv), fit = fit, cv = cv,',
  '    validation.rows = validation.rows)',
  '})',
  'saveRDS(list(version = as.character(packageVersion("precmed")),',
  '  imports = packageDescription("precmed")$Imports,',
  '  session = sessionInfo(), results = results), args[2])'
), worker)
outputs <- vapply(seq_along(libraries), function(i) tempfile(fileext = ".rds"), character(1))
for (i in seq_along(libraries)) {
  status <- system2(file.path(R.home("bin"), "Rscript"),
    shQuote(c(worker, libraries[i], outputs[i])))
  if (status != 0L) stop("Comparison worker failed for ", libraries[i])
}
baseline <- readRDS(outputs[1])
candidate <- readRDS(outputs[2])
stopifnot(identical(lapply(baseline$results, `[[`, "validation.rows"),
                    lapply(candidate$results, `[[`, "validation.rows")))
summary <- do.call(rbind, lapply(seq_along(baseline$results), function(i) {
  a <- baseline$results[[i]]
  b <- candidate$results[[i]]
  data.frame(seed = a$seed,
    score.correlation = cor(a$scores, b$scores),
    score.rmse = sqrt(mean((a$scores - b$scores)^2)),
    abc.main = mean(a$abc), abc.gbm3 = mean(b$abc),
    elapsed.main = a$elapsed, elapsed.gbm3 = b$elapsed)
}))
saveRDS(list(summary = summary, main = baseline, gbm3 = candidate), args[3])
print(summary, row.names = FALSE)
cat("Saved fitted results, CV results, and session information to", args[3], "\n")
unlink(c(worker, outputs))
