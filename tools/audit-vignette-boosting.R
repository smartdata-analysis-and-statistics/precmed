# Usage:
# Rscript tools/audit-vignette-boosting.R MAIN_LIBRARY GBM3_LIBRARY VIGNETTE_SOURCE_DIRECTORY OUTPUT_DIRECTORY [original|deterministic|both]
# Use fresh installations from the two branches and a shared dependency library.
# No cached vignette results are used. Unevaluated boosting examples are included.
# Original = historical fitting settings (with one core for reproducibility).
# Deterministic = both engines use bag.fraction=1, n.minobsinnode=50, and
# per-fit seeds with RNG restoration, so model fitting cannot change outer splits.
# Deterministic agreement does NOT reproduce the historical default outputs.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(4L, 5L)) stop("See usage at the top of this script.")
locations <- vapply(args[1:3], normalizePath, character(1), mustWork = TRUE)
dir.create(args[4], recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(args[4], mustWork = TRUE)
profile <- if (length(args) == 5L) match.arg(args[5], c("original", "deterministic", "both")) else "both"
self <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
helpers <- file.path(dirname(normalizePath(self)), "vignette-boosting")
run <- function(script, arguments) {
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    shQuote(c(file.path(helpers, script), arguments)))
  if (status != 0L) stop("Vignette audit step failed: ", script)
}
run("inventory.R", c(locations[3], output))
profiles <- if (profile == "both") c("original", "deterministic") else profile
for (setting in profiles) {
  paths <- file.path(output, paste0(c("gbm", "gbm3"), "-", setting, ".rds"))
  for (i in 1:2) run("worker.R", c(locations[i], c("gbm", "gbm3")[i], paths[i], setting))
  run("summarize.R", c(paths, file.path(output, paste0("outputs-", setting, ".csv"))))
  if (setting == "original") {
    run("replay.R", c(paths[1], locations[2], file.path(output, "configuration-replay.csv")))
  }
}
cat("Audit results saved to", output, "\n")
