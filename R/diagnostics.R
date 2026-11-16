#' Summarize recorded fitting diagnostics
#'
#' Collect treatment-recoding information from \code{atefit()}, or recorded
#' subgroup errors and warnings from count or survival \code{catefit()} and binary, count or survival
#' \code{catecv()}, into a data frame. Each row represents one recorded message.
#' Survival subgroup messages also remain in \code{object[["errors/warnings"]]}.
#'
#' @param type Choose \code{"messages"} (the default) for recorded issues or
#'   \code{"propensity"} for full-sample treatment-model diagnostics from a
#'   binary, count or survival \code{atefit} result. \code{"overlap"} is an
#'   alias for \code{"propensity"}.
#' @param layout For overlap diagnostics, \code{"formatted"} (the default)
#'   returns a compact, rounded character table; \code{"arms"} returns one
#'   numeric row per treatment; \code{"metrics"} returns one numeric column per
#'   treatment and one row per diagnostic. Not used for messages.
#' @param object A binary, count or survival \code{atefit}, a count or survival \code{catefit}, or a binary, count or survival \code{catecv} result.
#' @return A data frame with columns \code{method}, \code{split},
#'   \code{iteration}, \code{group}, \code{subgroup}, \code{type},
#'   \code{measure}, and \code{message}. For count mutually exclusive groups,
#'   \code{subgroup} gives the lower and upper boundaries; survival uses the
#'   upper percentile label recorded by the estimator.
#'   For \code{atefit}, bootstrap failures have split \code{bootstrap},
#'   the replicate number in iteration, and the affected estimate in measure.
#'   Split and iteration are \code{NA} for \code{catefit}.
#'   Otherwise, measure is \code{NA}
#'   unless the recorded message identifies it. Returns a zero-row data frame
#'   with the same columns when there are no recorded messages.
#' For \code{type = "propensity"}, the default formatted data frame has one row
#' per diagnostic and one character column per treatment. Its percentages use
#' the number of patients in the corresponding treatment group. Use
#' \code{layout = "arms"} or \code{layout = "metrics"} for unrounded numeric
#' results. The underlying diagnostics report sample size, untruncated propensity score range and median, number
#' outside the other arm's observed score range, number truncated below and
#' above the specified limits, and the median, maximum and effective sample
#' size of the inverse treatment probability weights after truncation. The
#' effective sample size is \eqn{(\sum w)^2 / \sum w^2}. These are treatment
#' weights only, not censoring weights or the complete weights of a particular
#' estimator. With \code{layout = "formatted"} or \code{layout = "metrics"},
#' treatment columns use the
#' original labels when available, or treatment 0 and treatment 1 otherwise.
#' Older fitted objects do not contain this information.
#' @details Overlap and weights describe the full analysis sample and the fitted
#' treatment model; they do not establish covariate balance, positivity, or
#' absence of unmeasured confounding. Extreme weights or little common support
#' merit investigation, but no universal cutoff guarantees validity.
#' The message view reports only diagnostics stored in the result, not every
#'   warning emitted during fitting. In particular, \code{atefit()} stores the
#'   treatment-recoding information and failed or undefined bootstrap estimates,
#'   but not all bootstrap warnings. Count
#'   fits record subgroup messages and support in \code{object$subgroup.support}.
#'   Continuous results are not supported by this helper.
#' @examples
#' \dontrun{
#' diagnostics(output_atefit)
#' diagnostics(output_atefit, type = "propensity")
#' diagnostics(output_atefit, type = "propensity", layout = "arms")
#' diagnostics(output_catecv)
#' subset(diagnostics(output_catecv), split == "validation")
#' }
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#' @export
diagnostics <- function(object, type = c("messages", "propensity", "overlap"),
                        layout = c("formatted", "arms", "metrics")) {
  type <- match.arg(type)
  layout <- match.arg(layout)
  if (type %in% c("propensity", "overlap")) {
    if (!inherits(object, "atefit") || is.null(object$overlap.diagnostics))
      stop("Overlap diagnostics are available for newly fitted binary, count or survival atefit objects only; refit older objects.")
    if (layout == "formatted") return(overlap_formatted_table(object))
    if (layout == "metrics") return(overlap_metrics_table(object))
    return(object$overlap.diagnostics)
  }
  if (!inherits(object, c("atefit", "precmed", "catefit")) ||
      !(object$response %in% c("survival", "count") ||
        (inherits(object, c("atefit", "precmed")) && identical(object$response, "binary"))))
    stop("object must be a binary, count or survival atefit, a count or survival catefit, or a binary, count or survival catecv result.")

  out <- data.frame(method = character(), split = character(), iteration = integer(),
                    group = character(), subgroup = character(), type = character(),
                    measure = character(), message = character(), stringsAsFactors = FALSE)
  if (inherits(object, "atefit")) {
    for (message in object$warning) {
      out[nrow(out) + 1L, ] <- list(NA_character_, NA_character_, NA_integer_,
                                     NA_character_, NA_character_, "info",
                                     NA_character_, trimws(message))
    }
    issues <- object$bootstrap.diagnostics
    if (!is.null(issues) && nrow(issues)) {
      for (i in seq_len(nrow(issues))) {
        out[nrow(out) + 1L, ] <- list(NA_character_, "bootstrap", issues$iteration[i],
          NA_character_, NA_character_, "error", issues$measure[i], issues$reason[i])
      }
    }
    return(out)
  }
  if (identical(object$response, "count") || identical(object$response, "binary")) {
    if (!is.null(object$diagnostic.messages)) return(object$diagnostic.messages)
    return(out)
  }
  if (is.null(object[["errors/warnings"]]))
    stop("object must be a count or survival atefit, catefit, or catecv result with recorded diagnostics.")
  add_messages <- function(entry, method, split, iteration, group) {
    for (type in c("errors", "warnings")) {
      entries <- entry[[type]]
      if (is.null(entries)) next
      for (i in seq_along(entries)) {
        label <- names(entries)[i]
        if (is.null(label) || is.na(label) || !nzchar(label)) label <- NA_character_
        messages <- as.character(entries[[i]])
        for (message in messages) {
          measure <- if (grepl("Cumulative hazard ratio", message, fixed = TRUE))
            "cumulative hazard ratio" else NA_character_
          out[nrow(out) + 1L, ] <<- list(method, split, iteration, group, label,
                                           if (type == "errors") "error" else "warning",
                                           measure, message)
        }
      }
    }
  }
  for (method in names(object[["errors/warnings"]])) {
    for (section in names(object[["errors/warnings"]][[method]])) {
      section_data <- object[["errors/warnings"]][[method]][[section]]
      if (is.null(section_data)) next
      group <- if (grepl("high", section, fixed = TRUE)) "high" else if (grepl("low", section, fixed = TRUE)) "low" else "group"
      if (identical(section, "est.high")) {
        add_messages(section_data, method, NA_character_, NA_integer_, group)
      } else {
        split <- if (grepl(".train.", section, fixed = TRUE)) "training" else "validation"
        for (iteration in names(section_data)) {
          add_messages(section_data[[iteration]], method, split,
                       as.integer(sub("^cv", "", iteration)), group)
        }
      }
    }
  }
  out
}

count_diagnostic_collector <- function() {
  e <- new.env(parent = emptyenv())
  e$messages <- data.frame(method = character(), split = character(), iteration = integer(),
    group = character(), subgroup = character(), type = character(), measure = character(),
    message = character(), stringsAsFactors = FALSE)
  e$support <- data.frame(method = character(), split = character(), iteration = integer(),
    group = character(), subgroup = character(), n0 = integer(), n1 = integer(),
    events0 = numeric(), events1 = numeric(), stringsAsFactors = FALSE)
  e$callback <- function(method, split, iteration, y, trt, include.support = TRUE) {
    force(method); force(split); force(iteration); force(y); force(trt); force(include.support)
    function(type, group, subgroup, selected, message) {
      if (type == "support") {
        if (include.support) e$support[nrow(e$support) + 1L, ] <- list(method, split, iteration, group, subgroup,
          sum(selected & trt == 0), sum(selected & trt == 1),
          sum(y[selected & trt == 0]), sum(y[selected & trt == 1]))
      } else e$messages[nrow(e$messages) + 1L, ] <- list(method, split, iteration,
        group, subgroup, type, "rate ratio", message)
    }
  }
  e
}

# Summarize the actual full-sample treatment probabilities and ATE weights.
# Author: Thomas Debray <tdebray@fromdatatowisdom.com>
treatment_overlap <- function(trt, propensity, minPS, maxPS) {
  raw <- attr(propensity, "raw")
  p <- as.numeric(propensity)
  if (is.null(raw)) stop("Untruncated propensity scores are unavailable.")
  do.call(rbind, lapply(0:1, function(arm) {
    idx <- trt == arm
    other <- raw[!idx]
    w <- if (arm == 1) 1 / p[idx] else 1 / (1 - p[idx])
    data.frame(treatment = arm, n = sum(idx),
      ps.min = min(raw[idx]), ps.median = stats::median(raw[idx]),
      ps.max = max(raw[idx]),
      outside.other.range = sum(raw[idx] < min(other) | raw[idx] > max(other)),
      truncated.low = sum(raw[idx] < minPS), truncated.high = sum(raw[idx] > maxPS),
      weight.median = stats::median(w), weight.max = max(w),
      effective.n = if (all(is.finite(w))) sum(w)^2 / sum(w^2) else NA_real_)
  }))
}

# Arrange overlap diagnostics for compact display without rounding values.
# Author: Thomas Debray <tdebray@fromdatatowisdom.com>
overlap_metrics_table <- function(object) {
  overlap <- object$overlap.diagnostics
  labels <- object$treatment.labels
  if (is.null(labels) || length(labels) != 2L || anyNA(labels) ||
      any(!nzchar(labels)) || anyDuplicated(labels) || any(labels == "Diagnostic"))
    labels <- c("treatment 0", "treatment 1")
  metrics <- c("Number of patients", "Propensity score: minimum",
    "Propensity score: median", "Propensity score: maximum",
    "Outside other group's range (n)", "Scores truncated below limit",
    "Scores truncated above limit", "Median treatment weight",
    "Maximum treatment weight", "Effective sample size")
  values <- as.matrix(overlap[, c("n", "ps.min", "ps.median", "ps.max",
    "outside.other.range", "truncated.low", "truncated.high",
    "weight.median", "weight.max", "effective.n")])
  result <- data.frame(Diagnostic = metrics, values[1, ], values[2, ],
                       check.names = FALSE, row.names = NULL)
  names(result) <- c("Diagnostic", as.character(labels))
  result
}

# Format full-sample overlap diagnostics for presentation.
# Author: Thomas Debray <tdebray@fromdatatowisdom.com>
overlap_formatted_table <- function(object) {
  raw <- object$overlap.diagnostics
  labels <- names(overlap_metrics_table(object))[-1L]
  format_arm <- function(i) {
    x <- raw[i, ]
    c(as.character(x$n),
      sprintf("%.2f (%.2f-%.2f)", x$ps.median, x$ps.min, x$ps.max),
      sprintf("%d (%.1f%%)", x$outside.other.range,
              100 * x$outside.other.range / x$n),
      sprintf("%d / %d", x$truncated.low, x$truncated.high),
      sprintf("%.2f (%.2f)", x$weight.median, x$weight.max),
      if (is.finite(x$effective.n)) sprintf("%.0f", x$effective.n) else "NA")
  }
  result <- data.frame(Diagnostic = c("Patients, n",
    "Propensity score, median (min-max)", "Outside other group's range, n (%)",
    "Truncated below / above, n", "Treatment weight, median (max)",
    "Effective sample size"), format_arm(1), format_arm(2),
    check.names = FALSE, stringsAsFactors = FALSE)
  names(result) <- c("Diagnostic", labels)
  result
}
