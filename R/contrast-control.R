#' @useDynLib precmed, .registration = TRUE
NULL

# Validate shared numerical controls before either backend is called.
validate_contrast_controls <- function(error.maxNR, max.iterNR, tune, optim.control) {
  if (!is.numeric(error.maxNR) || length(error.maxNR) != 1L ||
      !is.finite(error.maxNR) || error.maxNR <= 0) {
    stop("error.maxNR must be a finite positive number.")
  }
  if (!is.numeric(max.iterNR) || length(max.iterNR) != 1L ||
      !is.finite(max.iterNR) || max.iterNR < 0 ||
      max.iterNR != floor(max.iterNR) || max.iterNR > .Machine$integer.max) {
    stop("max.iterNR must be a non-negative integer.")
  }
  if (!is.numeric(tune) || length(tune) != 2L ||
      any(!is.finite(tune)) || any(tune <= 0)) {
    stop("tune must contain two finite positive numbers.")
  }
  allowed <- c("trace", "fnscale", "parscale", "ndeps", "maxit", "abstol",
               "reltol", "alpha", "beta", "gamma", "REPORT", "warn.1d.NelderMead")
  if (!is.list(optim.control) || (length(optim.control) &&
      (is.null(names(optim.control)) || anyNA(names(optim.control)) ||
       any(!nzchar(names(optim.control))) || anyDuplicated(names(optim.control)) ||
       any(!names(optim.control) %in% allowed)))) {
    stop("optim.control must be a named list of supported optim controls.")
  }
}

# Contrast-regression controls are separated from GBM's ... arguments.
resolve_contrast_control <- function(contrast.control, error.maxNR, max.iterNR, tune) {
  defaults <- list(solver = "C", error.maxNR = error.maxNR,
                   max.iterNR = max.iterNR, tune = tune, optim.control = list())
  if (!is.list(contrast.control) || (length(contrast.control) &&
      (is.null(names(contrast.control)) || anyNA(names(contrast.control)) ||
       any(!nzchar(names(contrast.control))) || anyDuplicated(names(contrast.control)) ||
       any(!names(contrast.control) %in% names(defaults))))) {
    stop("contrast.control must be a named list containing solver, error.maxNR, max.iterNR, tune or optim.control.")
  }
  defaults[names(contrast.control)] <- contrast.control
  defaults$solver <- match.arg(defaults$solver, c("C", "R"))
  validate_contrast_controls(defaults$error.maxNR, defaults$max.iterNR,
                             defaults$tune, defaults$optim.control)
  defaults
}
