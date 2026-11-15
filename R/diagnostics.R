#' Summarize recorded survival diagnostics
#'
#' Collect the treatment-recoding information stored by survival \code{atefit()},
#' or the subgroup errors and warnings stored by survival \code{catefit()} or
#' \code{catecv()}, into a data frame. Each row represents one recorded message.
#' The original subgroup list remains available as \code{object[["errors/warnings"]]}.
#'
#' @param object A survival \code{atefit}, \code{catefit}, or \code{catecv} result.
#' @return A data frame with columns \code{method}, \code{split},
#'   \code{iteration}, \code{group}, \code{subgroup}, \code{type},
#'   \code{measure}, and \code{message}. For mutually exclusive groups,
#'   \code{subgroup} is the upper percentile label recorded by the estimator.
#'   Method, split, iteration, group, subgroup, and measure are \code{NA} for
#'   \code{atefit}; split and iteration are \code{NA} for \code{catefit}.
#'   Otherwise, measure is \code{NA}
#'   unless the recorded message identifies it. Returns a zero-row data frame
#'   with the same columns when there are no recorded messages.
#' @details This reports only diagnostics stored in the result, not every
#'   warning emitted during fitting. In particular, \code{atefit()} stores the
#'   treatment-recoding information, but not all bootstrap warnings. Count and
#'   continuous results are not supported by this helper.
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
      !identical(object$response, "survival") ||
      (!inherits(object, "atefit") && is.null(object[["errors/warnings"]])))
    stop("object must be a survival atefit, catefit, or catecv result.")

  out <- data.frame(method = character(), split = character(), iteration = integer(),
                    group = character(), subgroup = character(), type = character(),
                    measure = character(), message = character(), stringsAsFactors = FALSE)
  if (inherits(object, "atefit")) {
    for (message in object$warning) {
      out[nrow(out) + 1L, ] <- list(NA_character_, NA_character_, NA_integer_,
                                     NA_character_, NA_character_, "info",
                                     NA_character_, trimws(message))
    }
    return(out)
  }
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
