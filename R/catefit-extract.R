#' Extract CATE scores and score coefficients
#'
#' Extract scores for the original analysis data, coefficients, and, where
#' available, confidence intervals for score coefficients from a \code{catefit}.
#'
#' @param object An object of class \code{catefit}.
#' @param method A single score to extract from the fitted object. Use
#' \code{"two_regression"} or
#' \code{"contrast_regression"} for those procedures, or the outcome method
#' (such as \code{"logistic"}, \code{"boosting"}, \code{"poisson"},
#' \code{"coxph"}, \code{"weibull"}, \code{"log_normal"},
#' \code{"negative_binomial"}, or \code{"random_forest"}) for a T-learner.
#' The established stored names \code{"twoReg"}, \code{"contrastReg"},
#' \code{"negBin"}, and \code{"randomForest"} are also accepted; contrast
#' regression can also be selected with \code{"contrastreg"}. If exactly
#' one T-learner score is present, \code{"t_learner"} selects it. If exactly
#' one score of any kind is present, \code{method} may be omitted from
#' \code{predict()}. For \code{coef()}, omission returns the coefficient table.
#' For \code{confint()}, the default is \code{"contrastReg"}. Coefficients
#' and confidence intervals are available only for certain score methods.
#' @param newdata New data are not currently supported; omit this argument to
#' retrieve scores for the original analysis sample in its original order.
#' @param parm Coefficient names or positions; omitted selects all coefficients.
#' @param level Confidence level between 0 and 1; defaults to 0.95.
#' @param ... Additional arguments; currently unused.
#' @return \code{predict()} returns a numeric vector of stored scores.
#' \code{coef()} returns the stored coefficient table or a named numeric vector
#' for a selected method. When score models use different covariates, the
#' table has the union of coefficient names and missing entries for terms not
#' used by a method; selecting a method returns only its terms.
#' \code{confint()} returns a matrix of Wald confidence
#' limits for the selected coefficients.
#' @details Scores for count and survival outcomes are on the log ratio scale;
#' binary scores are on the risk-difference scale. They describe fitted scores,
#' not validated individual treatment effects. Pass the fitted \code{catefit}
#' object to \code{predict()}, not an individual outcome model in its \code{fit}
#' component. When multiple T-learners are fitted, select one by its outcome
#' method. Confidence intervals require
#' stored coefficient standard errors, currently available for contrast
#' regression in count and continuous fits only. They are not available for
#' survival or binary fits. Methods that do not fit linear scores, such as
#' random forests, boosting, Cox, Weibull and log-normal T-learners, do not
#' provide linear score coefficients.
#' @author Thomas Debray
#' @name catefit-extract
NULL

# Normalize score selectors consistently across CATE extraction methods.
# Author: Thomas Debray
resolve_cate_score_method <- function(method, methods) {
  if (is.character(method) && length(method) == 1L && !is.na(method)) {
    aliases <- c(two_regression = "twoReg", contrast_regression = "contrastReg",
                 contrastreg = "contrastReg", negative_binomial = "negBin",
                 random_forest = "randomForest")
    if (method %in% names(aliases)) method <- unname(aliases[[method]])
    if (method == "t_learner") {
      candidates <- intersect(methods, c("logistic", "poisson", "negBin", "boosting", "randomForest", "gam", "gaussian", "coxph", "weibull", "log_normal"))
      if (length(candidates) != 1L)
        stop("Select one T-learner by its outcome method. Available T-learner scores: ",
             if (length(candidates)) paste(candidates, collapse = ", ") else "none", ".", call. = FALSE)
      method <- candidates
    }
  }
  method
}

# Author: Thomas Debray
catefit_score_methods <- function(object) {
  keys <- setdiff(grep("^score\\.", names(object), value = TRUE), "score.method")
  keys <- keys[!vapply(object[keys], is.null, logical(1))]
  sub("^score\\.", "", keys)
}

#' @rdname catefit-extract
#' @export
predict.catefit <- function(object, newdata = NULL, method = NULL, ...) {
  if (!is.null(newdata)) stop("Prediction for new data is not supported for catefit objects.")
  methods <- catefit_score_methods(object)
  if (is.null(method) && length(methods) == 1L) method <- methods
  method <- resolve_cate_score_method(method, methods)
  if (!is.character(method) || length(method) != 1L || is.na(method) || !method %in% methods)
    stop("Specify one available scoring method: ", paste(methods, collapse = ", "), ".")
  object[[paste0("score.", method)]]
}

#' @rdname catefit-extract
#' @export
coef.catefit <- function(object, method = NULL, ...) {
  cf <- object$coefficients
  if (is.null(method)) return(cf)
  methods <- catefit_score_methods(object)
  if (!length(methods)) methods <- setdiff(names(cf), "SE_contrastReg")
  method <- resolve_cate_score_method(method, methods)
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      is.null(cf) || !method %in% setdiff(names(cf), "SE_contrastReg"))
    stop("No score coefficients are available for the requested method.")
  result <- stats::setNames(cf[[method]], rownames(cf))
  if (!is.null(object$coefficient.terms[[method]])) result <- result[object$coefficient.terms[[method]]]
  result
}

#' @rdname catefit-extract
#' @export
confint.catefit <- function(object, parm, level = 0.95, method = "contrastReg", ...) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("level must be a single number strictly between 0 and 1.")
  methods <- catefit_score_methods(object)
  if (!length(methods)) methods <- setdiff(names(object$coefficients), "SE_contrastReg")
  method <- resolve_cate_score_method(method, methods)
  estimates <- coef(object, method = method)
  if (!identical(method, "contrastReg") || !object$response %in% c("count", "continuous") ||
      is.null(object$coefficients$SE_contrastReg))
    stop("Coefficient standard errors are not available for this method and outcome; confidence intervals cannot be calculated.")
  se <- stats::setNames(object$coefficients$SE_contrastReg, rownames(object$coefficients))[names(estimates)]
  if (missing(parm)) parm <- seq_along(estimates)
  if (!(is.character(parm) || is.numeric(parm)) || anyNA(parm) ||
      (is.character(parm) && any(!parm %in% names(estimates))) ||
      (is.numeric(parm) && (any(!is.finite(parm)) || any(parm != floor(parm)) || any(parm < 1 | parm > length(estimates)))))
    stop("parm must identify available coefficients by name or position.")
  z <- stats::qnorm(1 - (1 - level) / 2)
  result <- cbind(estimates[parm] - z * se[parm], estimates[parm] + z * se[parm])
  colnames(result) <- paste0(format(c((1 - level) / 2, 1 - (1 - level) / 2) * 100, trim = TRUE), "%")
  result
}
