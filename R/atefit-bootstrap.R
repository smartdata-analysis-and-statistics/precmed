# Record each undefined bootstrap metric with its replicate and reason.
# Author: Thomas Debray <tdebray@fromdatatowisdom.com>
ate_bootstrap_diagnostics <- function(estimates, errors, reasons = character()) {
  issues <- data.frame(iteration = integer(), measure = character(),
                       reason = character(), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(estimates))) {
    bad <- which(!is.finite(estimates[i, ]))
    if (!length(bad)) next
    for (j in bad) {
      measure <- colnames(estimates)[j]
      reason <- if (!is.na(errors[i])) errors[i] else if (measure %in% names(reasons))
        unname(reasons[measure]) else paste0("Non-finite ", measure, " estimate")
      issues[nrow(issues) + 1L, ] <- list(i, measure, reason)
    }
  }
  issues
}

# Report distinct reasons and affected estimates without suppressing failed samples.
# Author: Thomas Debray <tdebray@fromdatatowisdom.com>
ate_bootstrap_warn <- function(issues) {
  if (!nrow(issues)) return(invisible(NULL))
  summary <- unique(issues[c("measure", "reason")])
  detail <- vapply(seq_len(nrow(summary)), function(i) {
    same <- issues$measure == summary$measure[i] & issues$reason == summary$reason[i]
    paste0(summary$measure[i], ": ", summary$reason[i], " (", sum(same), ")")
  }, character(1))
  warning(length(unique(issues$iteration)), " bootstrap replicate(s) failed or had undefined estimates; ",
          "inference for affected measures is unavailable. Reasons by measure: ",
          paste(detail, collapse = "; "), ". See bootstrap.diagnostics for replicate details.", call. = FALSE)
  invisible(NULL)
}
