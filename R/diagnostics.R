#' Summarize recorded fitting diagnostics
#'
#' Collect treatment-recoding information from \code{atefit()}, or recorded
#' subgroup errors and warnings from count or survival \code{catefit()} and binary, count or survival
#' \code{catecv()}, into a data frame. Each row represents one recorded message.
#' Survival subgroup messages also remain in \code{object[["errors/warnings"]]}.
#'
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
#' @details This reports only diagnostics stored in the result, not every
#'   warning emitted during fitting. In particular, \code{atefit()} stores the
#'   treatment-recoding information and failed or undefined bootstrap estimates,
#'   but not all bootstrap warnings. Count
#'   fits record subgroup messages and support in \code{object$subgroup.support}.
#'   Continuous results are not supported by this helper.
#' @examples
#' \dontrun{
#' diagnostics(output_atefit)
#' diagnostics(output_catecv)
#' subset(diagnostics(output_catecv), split == "validation")
#' }
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#' @export
diagnostics <- function(object) {
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
