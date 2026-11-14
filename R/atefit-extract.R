#' Extract average treatment effect estimates and confidence intervals
#'
#' Extract treatment effects from an \code{atefit} object on their native
#' estimation scale: the log rate ratio for count outcomes, or the log RMTL
#' ratio, conditional log hazard ratio, and log cumulative hazard ratio for survival outcomes. Exponentiate the returned
#' estimates and interval limits to obtain ratios.
#'
#' @param object An object of class \code{atefit}.
#' @param ... Additional arguments; currently unused.
#' @return \code{coef()} returns a named numeric vector of treatment effect
#' estimates. \code{confint()} returns a matrix with one row per selected
#' treatment effect and columns for the lower and upper confidence limits.
#' @details Confidence intervals are normal (Wald) intervals calculated using
#' the bootstrap standard errors stored in the fitted object. Changing
#' \code{level} does not repeat the bootstrap. Treatment-specific event rates
#' and restricted mean survival times are not included.
#' @examples
#' \dontrun{
#' fit <- atefitcount(data = countExample,
#'   cate.model = y ~ age + offset(log(years)), ps.model = trt ~ age)
#' exp(coef(fit))
#' exp(confint(fit))
#' exp(confint(fit, level = 0.90))
#' }
#' @name atefit-extract
NULL

atefit_effect_terms <- function(object) {
  switch(object$response,
         count = "log.rate.ratio",
         survival = c("log.rmtl.ratio", "log.hazard.ratio",
                      if (!is.null(object$log.cumulative.hazard.ratio)) "log.cumulative.hazard.ratio"),
         stop("Unsupported response in atefit object."))
}

#' @rdname atefit-extract
#' @export
coef.atefit <- function(object, ...) {
  terms <- atefit_effect_terms(object)
  vapply(terms, function(term) object[[term]]$estimate, numeric(1))
}

#' @rdname atefit-extract
#' @param parm A character vector of effect names or a numeric vector of
#' effect positions. If omitted, all treatment effects are included.
#' @param level Confidence level, a single number strictly between 0 and 1.
#' Default is \code{0.95}.
#' @export
confint.atefit <- function(object, parm, level = 0.95, ...) {
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("level must be a single number strictly between 0 and 1.")
  }
  estimates <- coef(object)
  terms <- names(estimates)
  if (missing(parm)) parm <- terms
  if (is.numeric(parm)) {
    if (anyNA(parm) || any(!is.finite(parm) | parm != trunc(parm) |
                           parm < 1 | parm > length(terms))) {
      stop("parm must contain valid effect names or positions.")
    }
    parm <- terms[parm]
  }
  if (!is.character(parm) || anyNA(parm) || any(!parm %in% terms)) {
    stop("parm must contain valid effect names or positions.")
  }
  se <- vapply(parm, function(term) object[[term]]$SE, numeric(1))
  probabilities <- c((1 - level) / 2, (1 + level) / 2)
  limits <- outer(se, qnorm(probabilities), "*") + estimates[parm]
  dimnames(limits) <- list(parm, paste0(format(100 * probabilities, trim = TRUE), " %"))
  limits
}
