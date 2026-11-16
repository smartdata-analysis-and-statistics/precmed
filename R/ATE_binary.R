#' Estimate an average treatment effect for a binary outcome
#'
#' Estimate the population-average risk difference between two treatments
#' using augmented inverse probability weighting (AIPW). Bootstrap resampling
#' provides a standard error, confidence interval, and p-value.
#'
#' @param data A data frame containing the outcome, treatment, and covariates.
#' @param cate.model A formula with the binary outcome on the left and outcome
#'   model covariates on the right. The outcome must be numeric or logical and
#'   coded 0/1; missing values and offsets are not supported.
#' @param ps.model A formula with treatment on the left and treatment-model
#'   covariates on the right. The model estimates the propensity score. Both treatment groups must be present. Numeric
#'   0/1 or logical FALSE/TRUE map to treatment 0/1; other two-valued labels are
#'   sorted and mapped to 0/1 in that order. Use \code{trt ~ 1} for an
#'   intercept-only model.
#' @param ps.method Propensity score method: \code{"glm"} (default) or
#'   \code{"lasso"}, using the package's propensity-score LASSO procedure.
#' @param minPS,maxPS Lower and upper truncation limits for propensity scores.
#'   Defaults are 0.01 and 0.99; \code{minPS} must be smaller than \code{maxPS}.
#' @param interactions If \code{TRUE} (default), include interactions between
#'   treatment and all outcome-model covariates. If \code{FALSE}, fit a pooled
#'   logistic outcome model with a treatment main effect and no interactions.
#' @param n.boot Number of bootstrap samples, an integer of at least 2.
#'   Default is 500.
#' @param seed An optional integer seed; \code{NULL} uses the current random
#'   number state. Supplying a seed sets R's random-number state.
#' @param verbose Accepted for compatibility with the other ATE functions;
#'   the binary implementation does not currently display progress messages.
#' @param ps.control Named options for LASSO propensity-score interaction
#'   screening. Supported entry: \code{solver = "C"} (default) or \code{"R"}.
#'
#' @return An object of class \code{atefit} containing:
#' \itemize{
#'   \item \code{risk.difference}: A data frame with the estimate for treatment
#'   1 minus treatment 0, bootstrap standard error, 95\% Wald confidence limits,
#'   and two-sided Wald p-value.
#'   \item \code{risk0}, \code{risk1}: Data frames with estimated
#'   treatment-specific response probabilities.
#'   \item \code{trt.boot}: Bootstrap risk differences; \code{risk.boot}:
#'   bootstrap treatment-specific risks; \code{n.boot}: number of replicates.
#'   \item \code{treatment.labels}: Original treatment labels in 0/1 order.
#'   \item \code{bootstrap.diagnostics}: A data frame containing the replicate
#'   number, affected measure and reason for each failed or non-finite estimate.
#' }
#'
#' @details The estimator combines predicted risks from a logistic outcome model
#' with inverse-propensity-weighted residual corrections. Under appropriate
#' identification and regularity assumptions, it is consistent if either the
#' outcome model or treatment model is correctly specified. A causal
#' interpretation also requires consistency, no unmeasured confounding and
#' adequate treatment overlap. Propensity score truncation can affect the
#' guarantee based on the treatment model. AIPW risk estimates need not lie
#' in [0, 1]. Outcomes should be fully observed for the intended follow-up
#' horizon; this function does not adjust for censoring.
#'
#' Each bootstrap replicate refits both models. If any replicate for the risk
#' difference fails or is non-finite, inference is unavailable; failures are
#' summarized in a warning and stored in \code{bootstrap.diagnostics}.
#' \code{diagnostics()} can be used to inspect the records. The point estimate
#' and standard error are returned by \code{coef()} and through the result
#' object; \code{confint()} calculates Wald intervals at a requested level.
#'
#' @references Hoogland, J., IntHout, J., Belias, M., et al. (2021). A tutorial
#' on individualized treatment effect prediction from randomized trials with a
#' binary endpoint. \emph{Statistics in Medicine}, 40(26), 5961--5981.
#' \doi{10.1002/sim.9154}
#'
#' @seealso \code{\link{atefit}}, \code{\link{simulate_binary}},
#'   \code{\link{diagnostics}}
#' @examples
#' d <- simulate_binary(n = 500, seed = 42)
#' fit <- atefitbinary(d, cate.model = y ~ age + female,
#'                     ps.model = trt ~ age, n.boot = 10, seed = 42)
#' coef(fit)
#' confint(fit)
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#' @export
# Author: Thomas Debray <tdebray@fromdatatowisdom.com>
atefitbinary <- function(data, cate.model, ps.model, ps.method = "glm", minPS = 0.01, maxPS = 0.99, interactions = TRUE, n.boot = 500, seed = NULL, verbose = 0, ps.control = list()) {
  resolve_ps_control(ps.control)
  if (!ps.method %in% c("glm", "lasso") || minPS < 0 || maxPS > 1 || minPS >= maxPS) stop("Invalid propensity score settings.")
  if (length(n.boot) != 1L || is.na(n.boot) || n.boot < 2 || n.boot != floor(n.boot)) stop("n.boot must be an integer of at least 2.")
  if (!is.logical(interactions) || length(interactions) != 1L || is.na(interactions)) stop("interactions must be TRUE or FALSE.")
  d <- binary_data(data, cate.model, ps.model)
  est <- binary_effect(d, ps.method, minPS, maxPS, ps.control, interactions)
  if (!is.null(seed)) set.seed(seed)
  boot <- matrix(NA_real_, nrow = n.boot, ncol = 3L,
                 dimnames = list(NULL, c("risk0", "risk1", "risk.difference")))
  errors <- rep(NA_character_, n.boot)
  for (i in seq_len(n.boot)) {
    sampled <- binary_subset(d, sample.int(length(d$y), replace = TRUE))
    value <- tryCatch(binary_effect(sampled, ps.method, minPS, maxPS, ps.control, interactions),
                      error = function(e) { errors[i] <<- conditionMessage(e); NULL })
    if (!is.null(value)) boot[i, ] <- value[colnames(boot)]
  }
  issues <- ate_bootstrap_diagnostics(boot, errors)
  ate_bootstrap_warn(issues)
  if (!is.finite(est[3])) warning("The original-sample risk difference is non-finite; its inference is unavailable.", call. = FALSE)
  se <- if (!is.finite(est[3]) || any(!is.finite(boot[, 3]))) NA_real_ else stats::sd(boot[, 3])
  effect <- data.frame(estimate = unname(est[3]), SE = se, CI.lower = est[3] - qnorm(.975) * se, CI.upper = est[3] + qnorm(.975) * se, pvalue = if (is.na(se) || se == 0) NA_real_ else 2 * pnorm(-abs(est[3] / se)), row.names = "risk.difference")
  structure(list(response = "binary", risk.difference = effect, risk0 = data.frame(estimate = unname(est[1])), risk1 = data.frame(estimate = unname(est[2])), n.boot = n.boot, trt.boot = boot[, 3], risk.boot = data.frame(risk0 = boot[, 1], risk1 = boot[, 2]), bootstrap.diagnostics = issues, treatment.labels = d$labels), class = "atefit")
}

