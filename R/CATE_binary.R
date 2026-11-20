# Fit Bernoulli boosting separately in each treatment arm on a common design.
# Author: Thomas Debray <tdebray@fromdatatowisdom.com>
binary_boosting <- function(d, tree.depth, n.trees.boosting, plot.gbmperf, ...) {
  if (ncol(d$x) < 2L) stop("Binary boosting requires at least one outcome-model covariate.")
  if (!is.numeric(tree.depth) || length(tree.depth) != 1L || !is.finite(tree.depth) ||
      tree.depth < 1 || tree.depth != floor(tree.depth)) stop("tree.depth must be a positive integer.")
  if (!is.numeric(n.trees.boosting) || length(n.trees.boosting) != 1L || !is.finite(n.trees.boosting) ||
      n.trees.boosting < 1 || n.trees.boosting != floor(n.trees.boosting)) stop("n.trees.boosting must be a positive integer.")
  if (!is.logical(plot.gbmperf) || length(plot.gbmperf) != 1L || is.na(plot.gbmperf)) stop("plot.gbmperf must be TRUE or FALSE.")
  predictors <- as.data.frame(d$x[, -1, drop = FALSE])
  names(predictors) <- paste0("x", seq_len(ncol(predictors)))
  train <- data.frame(y = d$y, predictors, check.names = FALSE)
  fits <- vector("list", 2L); best <- integer(2L); predictions <- vector("list", 2L)
  for (arm in 0:1) {
    selected <- d$trt == arm
    if (sum(selected) < 5L || length(unique(d$y[selected])) != 2L)
      stop("Binary boosting requires at least five patients and both outcome values in each treatment group.")
    fits[[arm + 1L]] <- tryCatch(
      fit_boosting(y ~ ., data = train[selected, , drop = FALSE], distribution = "bernoulli",
                   interaction.depth = tree.depth, n.trees = n.trees.boosting, cv.folds = 5, ...),
      error = function(e) stop("Binary boosting failed for treatment ", arm, ": ", conditionMessage(e), call. = FALSE))
    best[[arm + 1L]] <- min(n.trees.boosting, max(10, best_boosting_iteration(fits[[arm + 1L]], "cv", plot.it = plot.gbmperf)))
    predictions[[arm + 1L]] <- as.numeric(predict(fits[[arm + 1L]], newdata = predictors,
                                                  n.trees = best[[arm + 1L]], type = "response"))
  }
  if (any(!is.finite(unlist(predictions))) || any(unlist(predictions) < 0 | unlist(predictions) > 1))
    stop("Binary boosting produced invalid predicted probabilities.")
  list(mu0 = predictions[[1]], mu1 = predictions[[2]],
       fit = list(fit0.boosting = fits[[1]], fit1.boosting = fits[[2]],
                  best0.iter = best[[1]], best1.iter = best[[2]],
                  predictor.names = colnames(d$x)[-1]))
}

#' Estimate conditional treatment-effect scores for binary outcomes
#'
#' Fit binary-outcome CATE scores and estimate risk differences in nested
#' subgroups. This is the binary-outcome interface used by \code{\link{catefit}}.
#'
#' @param data A data frame containing the variables in the outcome and treatment models; a data frame with \code{n} rows (1 row per observation).
#' @param outcome_model A formula specifying the outcome and baseline covariates.
#' @param outcome_method Outcome-model method paired with \code{benefit_method}.
#' For a T-learner, \code{"logistic"} (default) or \code{"boosting"};
#' for two-regression, \code{"logistic"}.
#' @param benefit_model Optional one-sided formula specifying the covariates and
#' terms for the two-regression or contrast-regression score. The default
#' \code{NULL} uses the covariates and terms of \code{outcome_model}.
#' An intercept is required and offsets are not supported. T-learner scores use
#' the outcome model directly and do not use this formula. In a mixed fit it
#' applies only to two-regression and contrast regression. Continuous-outcome
#' fitting is not implemented.
#' @param benefit_method Score-construction methods; defaults to \code{"t_learner"}.
#' See \code{\link{catefit}}.
#' @param ps_model A formula describing the treatment model (propensity-score model) to be fitted. The treatment must
#' appear on the left-hand side. For binary outcomes, logical treatment values
#' \code{FALSE} and \code{TRUE} are mapped to 0 and 1, respectively.
#' Numeric 0/1 values or two distinct treatment labels are also accepted.
#' If data are from a randomized controlled trial, specify \code{ps_model = trt ~ 1} as an intercept-only model.
#' @param ps_method A character value for the method to estimate the propensity score.
#' Allowed values include one of:
#' \code{'glm'} for logistic regression with main effects only (default), or
#' \code{'lasso'} for a logistic regression with main effects and LASSO penalization on
#' two-way interactions (added to the model if interactions are not specified in \code{ps_model}).
#' Relevant only when \code{ps_model} has more than one variable.
#' @param minPS A numerical value (in `[0, 1]`) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in `(0, 1]`) above which estimated propensity scores should be
#' truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.
#' @param higher.y A logical value indicating whether higher (\code{TRUE}) or
#' lower (\code{FALSE}) values of the outcome are more desirable. Default is \code{TRUE}.
#' @param prop.cutoff A vector of numerical values (in `(0, 1]`) specifying percentiles of the
#' estimated risk-difference CATE scores to define nested subgroups. Each element represents the cutoff to
#' separate observations in nested subgroups (below vs above cutoff).
#' The length of \code{prop.cutoff} is the number of nested subgroups.
#' An equally-spaced sequence of proportions ending with 1 is recommended.
#' Default is \code{seq(0.5, 1, length = 6)}.
#' @param seed An optional integer specifying an initial randomization seed for reproducibility.
#' Default is \code{NULL}, corresponding to no seed.
#' @param ps.control Named options for LASSO propensity-score interaction screening.
#' Supported entry: \code{solver = "C"} (default) or \code{solver = "R"}.
#' Only affects screening when \code{ps_method = "lasso"}; the LASSO fit and
#' interaction selection rule are unchanged. Kept separate from boosting's
#' \code{...} and other numerical controls.
#' @param B Number of repetitions of cross-fitting for \code{benefit_method =
#' 'two_regression'} or \code{'contrast_regression'}. Default is \code{3}.
#' @param Kfold A positive integer specifying the number of folds used in cross-fitting
#' to partition the data in \code{benefit_method = 'two_regression'} and \code{'contrast_regression'}.
#' Default is \code{5}.
#' @param tree.depth A positive integer specifying the maximum depth of boosting trees.
#' Used when \code{outcome_method = 'boosting'}. Default is \code{2}.
#' @param n.trees.boosting The maximum number of boosting trees, used when
#' \code{outcome_method = 'boosting'}. Default is \code{200}.
#' @param plot.gbmperf Whether to plot boosting performance measures when
#' \code{outcome_method = 'boosting'}. Default is \code{FALSE}.
#' @param ... Additional boosting options. See \code{\link{boosting-control}}.
#' Deprecated argument names from precmed 1.1.0 remain supported for
#' backward compatibility and emit warnings.
#' @return A \code{catefit} object containing treatment labels, fitted models,
#'   predicted risk-difference scores and nested subgroup risk differences for
#'   each requested method. Linear \code{twoReg} score coefficients are stored
#'   in \code{coefficients}.
#' @details
#' See the Binary outcomes section of \code{\link{catefit}} for the score
#' methods, their assumptions, and interpretation. Subgroup effects are estimated
#' on the fitting data and do not provide independent validation.
#' @examples
#' \donttest{
#' dat <- simulate_binary(n = 400, seed = 2026)
#' fit <- catefitbinary(data = dat,
#'                      outcome_model = y ~ age + previous_treatment,
#'                      outcome_method = "logistic",
#'                      benefit_method = "t_learner",
#'                      ps_model = trt ~ age + previous_treatment,
#'                      seed = 2026)
#' }
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#' @export
catefitbinary <- function(data, outcome_model = NULL, outcome_method = "logistic", benefit_model = NULL, benefit_method = "t_learner", ps_model, ps_method = "glm", minPS = 0.01, maxPS = 0.99, higher.y = TRUE, prop.cutoff = seq(0.5, 1, length = 6), seed = NULL, ps.control = list(), B = 3, Kfold = 5,
                          tree.depth = 2, n.trees.boosting = 200, plot.gbmperf = FALSE, ...) {
  previous_call <- rematch_pre_benefit_call(sys.function(), sys.call(), parent.frame())
  if (!is.null(previous_call)) return(do.call(sys.function(), previous_call, envir = parent.frame()))
  old_call <- rematch_catefit_legacy_call(sys.function(), sys.call(), parent.frame())
  if (!is.null(old_call)) return(do.call(sys.function(), old_call, envir = parent.frame()))
  .args <- as.list(match.call())[-1L]
  aliases <- intersect(names(.args), c("ps.model", "ps.method"))
  .args <- resolve_cate_aliases(.args)
  if (length(aliases)) return(do.call(sys.function(), .args, envir = parent.frame()))
  if (!missing(ps_model)) ps.model <- ps_model
  ps.method <- ps_method
  .dots <- list(...)
  if (missing(outcome_method) && !identical(benefit_method, "t_learner")) outcome_method <- NULL
  check_cate_initial_aliases(.args)
  initial_missing <- !"initial.predictor.method" %in% names(.dots)
  initial <- NULL
  if (!initial_missing) {
    initial <- .dots$initial.predictor.method
    if (is.null(initial)) initial <- "logistic"
    if (!is.null(initial) && !identical(initial, "logistic"))
      stop("Binary initial outcome models support only logistic regression.", call. = FALSE)
    .dots$initial.predictor.method <- NULL
  }
  init.model <- .dots$init.model
  .dots$init.model <- NULL
  if (!is.null(init.model) && !is.null(benefit_model)) stop("Specify only benefit_model or init.model, not both.", call. = FALSE)
  config <- resolve_catefit_arguments("binary", outcome_model, outcome_method, benefit_method, NULL, NULL, initial, initial_missing)
  cate.model <- config$outcome_model
  score.method <- config$scores
  check_score_method("binary", score.method)
  if (!ps.method %in% c("glm", "lasso") || minPS < 0 || maxPS > 1 || minPS >= maxPS) stop("Invalid propensity score settings.")
  resolve_ps_control(ps.control)
  if (!is.logical(higher.y) || length(higher.y) != 1L || is.na(higher.y)) stop("higher.y must be TRUE or FALSE.")
  if (!is.numeric(prop.cutoff) || !length(prop.cutoff) || anyNA(prop.cutoff) || any(prop.cutoff <= 0 | prop.cutoff > 1)) stop("prop.cutoff must contain values in (0, 1].")
  d <- binary_data(data, cate.model, ps.model, if ("twoReg" %in% score.method) init.model else NULL)
  x.benefit <- benefit_design(data, benefit_model, d$x, intercept = TRUE)
  if (!is.null(seed)) set.seed(seed)
  prop <- sort(prop.cutoff); out <- list(response = "binary", score.method = score.method, fit = list(), coefficients = NULL, treatment.labels = d$labels)
  for (method in score.method) {
    pred <- if (method == "twoReg") binary_tworeg_dr(within(d, x <- x.benefit), ps.method, minPS, maxPS, ps.control, B, Kfold) else if (method == "boosting") do.call(binary_boosting, c(list(d = d, tree.depth = tree.depth, n.trees.boosting = n.trees.boosting, plot.gbmperf = plot.gbmperf), .dots)) else
      binary_predictions(d$y, d$trt, d$x, "logistic_arms")
    score <- if (method == "twoReg") pred$score else pred$mu1 - pred$mu0
    out[[paste0("score.", method)]] <- score
    out$fit[[method]] <- pred$fit
    if (method == "twoReg") {
      beta <- pred$fit$coefficients
      out$coefficients <- data.frame(twoReg = beta[, "treated"] - beta[, "control"],
                                     row.names = rownames(beta))
    }
    ranking <- if (higher.y) score else -score
    effects <- vapply(prop, function(p) {
      selected <- ranking >= unname(stats::quantile(ranking, 1 - p))
      tryCatch(unname(binary_effect(binary_subset(d, selected), ps.method, minPS, maxPS, ps.control)[3]), error = function(e) { warning("Subgroup effect unavailable: ", conditionMessage(e)); NA_real_ })
    }, numeric(1))
    names(effects) <- paste0("prop", round(prop, 2))
    out[[paste0("ate.", method)]] <- effects
  }
  structure(out, class = "catefit")
}
