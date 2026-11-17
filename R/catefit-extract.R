#' Extract CATE scores and score coefficients
#'
#' Extract scores for the original analysis data, coefficients, and, where
#' available, confidence intervals for score coefficients from a \code{catefit}.
#'
#' @param object An object of class \code{catefit}.
#' @param method A single scoring method used to fit the object, such as
#' \code{"contrastReg"}. For \code{predict()}, it may be omitted if exactly
#' one method is available. For \code{coef()}, omission returns the stored
#' coefficient table. For \code{confint()}, the default is \code{"contrastReg"}.
#' @param newdata New data are not currently supported; omit this argument to
#' retrieve scores for the original analysis sample in its original order.
#' @param parm Coefficient names or positions; omitted selects all coefficients.
#' @param level Confidence level between 0 and 1; defaults to 0.95.
#' @param ... Additional arguments; currently unused.
#' @return \code{predict()} returns a numeric vector of stored scores.
#' \code{coef()} returns the stored coefficient table or a named numeric vector
#' for a selected method. \code{confint()} returns a matrix of Wald confidence
#' limits for the selected coefficients.
#' @details Scores for count and survival outcomes are on the log ratio scale;
#' binary scores are on the risk-difference scale. They describe fitted scores,
#' not validated individual treatment effects. Confidence intervals require
#' stored coefficient standard errors, currently available for contrast
#' regression in count and continuous fits only. They are not available for
#' survival or binary fits. Methods that do not fit linear scores, such as
#' random forests and boosting, do not provide score coefficients.
#' @author Thomas Debray
#' @name catefit-extract
NULL

#' @rdname catefit-extract
#' @export
predict.catefit <- function(object, newdata = NULL, method = NULL, ...) {
  if (!is.null(newdata)) stop("Prediction for new data is not supported for catefit objects.")
  keys <- setdiff(grep("^score\\.", names(object), value = TRUE), "score.method")
  keys <- keys[!vapply(object[keys], is.null, logical(1))]
  methods <- sub("^score\\.", "", keys)
  if (is.null(method) && length(methods) == 1L) method <- methods
  if (!is.character(method) || length(method) != 1L || is.na(method) || !method %in% methods)
    stop("Specify one available scoring method: ", paste(methods, collapse = ", "), ".")
  object[[paste0("score.", method)]]
}

#' @rdname catefit-extract
#' @export
coef.catefit <- function(object, method = NULL, ...) {
  cf <- object$coefficients
  if (is.null(method)) return(cf)
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      is.null(cf) || !method %in% setdiff(names(cf), "SE_contrastReg"))
    stop("No score coefficients are available for the requested method.")
  stats::setNames(cf[[method]], rownames(cf))
}

#' @rdname catefit-extract
#' @export
confint.catefit <- function(object, parm, level = 0.95, method = "contrastReg", ...) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("level must be a single number strictly between 0 and 1.")
  estimates <- coef(object, method = method)
  if (!identical(method, "contrastReg") || !object$response %in% c("count", "continuous") ||
      is.null(object$coefficients$SE_contrastReg))
    stop("Coefficient standard errors are not available for this method and outcome; confidence intervals cannot be calculated.")
  se <- stats::setNames(object$coefficients$SE_contrastReg, names(estimates))
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
