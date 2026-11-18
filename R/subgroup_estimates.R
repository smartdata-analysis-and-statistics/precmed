#' Extract estimates for score-defined subgroups
#'
#' Extract treatment effects or treatment-specific outcomes for the nested
#' subgroups in a binary, count or survival \code{catefit} object.
#'
#' @param object A binary, count or survival \code{catefit} object.
#' @param method One scoring method used to fit the object, such as
#' \code{"contrast_regression"} or a T-learner outcome method such as
#' \code{"poisson"}. The same method selectors as \code{predict()} are
#' accepted; see \code{\link{catefit-extract}}. If omitted, it is inferred
#' when exactly one method is available.
#' @param measure For binary outcomes, \code{"risk_difference"} (risk under
#' treatment 1 minus treatment 0). For count outcomes, \code{"rate_ratio"}. For survival outcomes,
#' one of \code{"rmtl"} (RMTL ratio), \code{"cumulative_hr"}
#' (marginal cumulative hazard ratio), \code{"conditional_hr"} (conditional
#' Cox hazard ratio), \code{"risk0"}, \code{"risk1"}, \code{"risk_difference"}
#' (risk under treatment 1 minus treatment 0), \code{"absolute_risk_reduction"}
#' (the opposite contrast), \code{"rmst0"}, or \code{"rmst1"} at \code{tau0}.
#' The default is \code{"risk_difference"} for binary, \code{"rate_ratio"} for count,
#' and \code{"rmtl"} for survival.
#' @return A named numeric vector, with one estimate per nested subgroup
#' in increasing order of \code{prop.cutoff}. Treatment-specific risks and RMSTs are
#' outcomes, not treatment effects. Undefined estimates are \code{NA}.
#' @details The subgroups are ranked by the selected score and estimated using
#' the same data used to fit it; these results do not validate the score. The
#' cumulative hazard ratio and risk estimates use outcome, treatment, and
#' censoring models; a causal interpretation requires appropriate treatment
#' and censoring assumptions. Older fits without the requested estimates must
#' be refitted. No confidence intervals are returned. The
#' measure names are also accepted in dotted form.
#' @author Thomas Debray
#' @export
subgroup_estimates <- function(object, method = NULL, measure = NULL) {
  if (!inherits(object, "catefit") || !object$response %in% c("binary", "count", "survival"))
    stop("object must be a binary, count or survival catefit object.")
  if (is.null(measure)) measure <- switch(object$response, binary = "risk_difference",
                                           count = "rate_ratio", survival = "rmtl")
  methods <- catefit_score_methods(object)
  if (is.null(method) && length(methods) == 1L) method <- methods
  method <- resolve_cate_score_method(method, methods)
  if (!is.character(method) || length(method) != 1L || is.na(method) || !method %in% methods)
    stop("Specify one available scoring method: ", paste(methods, collapse = ", "), ".")
  aliases <- c(rate.ratio = "rate_ratio", cumulative.hr = "cumulative_hr", conditional.hr = "conditional_hr",
               risk.difference = "risk_difference", absolute.risk.reduction = "absolute_risk_reduction")
  if (is.character(measure) && length(measure) == 1L && !is.na(measure) && measure %in% names(aliases))
    measure <- unname(aliases[[measure]])
  prefixes <- if (object$response == "binary") c(risk_difference = "ate") else
    if (object$response == "count") c(rate_ratio = "ate") else
    c(rmtl = "ate", cumulative_hr = "chr", conditional_hr = "hr",
                risk0 = "risk0", risk1 = "risk1", risk_difference = "rd",
                absolute_risk_reduction = "rd", rmst0 = "rmst0", rmst1 = "rmst1")
  if (!is.character(measure) || length(measure) != 1L || is.na(measure) || !measure %in% names(prefixes))
    stop("measure must be one of: ", paste(names(prefixes), collapse = ", "), ".")
  value <- object[[paste0(prefixes[[measure]], ".", method)]]
  if (is.null(value)) stop("The requested estimates are unavailable; refit the object if it predates this measure.")
  if (measure == "absolute_risk_reduction") value <- -value
  value
}

