#' Validate predicted treatment effects
#'
#' This experimental internal function is under development and issues a warning
#' when called. It is not called by
#' \code{catecv()}. Summarize calibration of predicted treatment effects and the population
#' effect of a treatment rule using independent validation data.
#'
#' @param y Numeric observed binary or continuous outcome.
#' @param trt Numeric treatment assignment, coded 0 for control and 1 for treatment.
#' @param benefit Numeric predicted treatment effect for each participant, on the
#'   outcome difference scale (treatment minus control).
#' @param groups Positive integer number of quantile groups for calibration.
#'   Tied predictions are kept together, so fewer groups may result.
#' @param threshold Numeric treatment-effect threshold on the same scale as
#'   \code{benefit}. Treatment is recommended when predicted benefit exceeds it.
#' @param propensity Optional numeric vector of treatment probabilities, one per
#'   participant, strictly between 0 and 1. If supplied, inverse probability
#'   weighting is used; if NULL, unadjusted arm means are used.
#' @param higher.y Whether larger outcomes are desirable. If \code{FALSE}, the
#'   treatment rule recommends treatment when the predicted difference is less
#'   than \code{threshold}. Effects are still reported as treatment minus control.
#'
#' @return A list with \code{mean.bias} (observed treatment difference minus
#'   mean prediction), \code{calibration} (group sizes, arm counts, mean
#'   predictions, observed differences and standard errors), \code{calibration.fit}
#'   (group-level intercept, slope, R-squared and RMSE), and \code{policy}
#'   (rule treatment proportion and outcome differences for the rule versus
#'   its opposite, control for all, and treatment for all). Policy differences
#'   are oriented so positive values indicate better outcomes under the rule.
#'
#' @details Unadjusted estimators assume randomized treatment assignment with a
#' common, positive probability of each treatment. With propensity scores,
#' identification requires conditional exchangeability given the treatment-model
#' covariates, positivity and a correctly specified treatment model. Both
#' approaches require consistency and no interference for causal interpretation.
#' Predictions should be obtained from independent data or appropriate
#' cross-validation. Censoring and treatment costs are not adjusted for. A threshold changes the treatment rule, but policy results
#' remain outcome differences and do not subtract treatment costs.
#'
#' Calibration groups are formed from empirical quantile cut points; ties are
#' never split. Group summaries are differences in arm means, using normalized
#' inverse probability weights when propensity scores are provided. The fitted line and RMSE give equal weight to each group and describe
#' group-level approximations, not individual treatment effect prediction
#' error. A group lacking either treatment arm has an unavailable observed
#' difference; in that case the calibration fit is unavailable. The policy
#' estimates standardize arm-specific outcomes within the two rule categories
#' to their observed proportions, allowing unequal randomization probabilities.
#' A required empty treatment arm makes the corresponding policy result
#' unavailable. Unadjusted standard errors require at least two observations in
#' each arm. Weighted standard errors are unavailable: uncertainty in estimated
#' treatment models is not accounted for by the unadjusted formula. Truncated
#' or misspecified propensity scores can bias weighted estimates.
#'
#' @references Hoogland et al. (2021), \doi{10.1002/sim.9154}; Efthimiou et al.
#'   (2023), \doi{10.1002/sim.9665}.
#' @author Thomas Debray
#' @keywords internal
validate_benefit <- function(y, trt, benefit, groups = 5, threshold = 0,
                             higher.y = TRUE, propensity = NULL) {
  n <- length(y)
  if (!is.numeric(y) || !is.numeric(trt) || !is.numeric(benefit) ||
      n < 2L || length(trt) != n || length(benefit) != n ||
      any(!is.finite(y)) || any(!is.finite(trt)) || any(!is.finite(benefit)))
    stop("y, trt and benefit must be finite numeric vectors of equal length (at least two).")
  if (!all(trt %in% 0:1) || length(unique(trt)) != 2L)
    stop("trt must contain both treatment arms, coded 0 and 1.")
  if (length(groups) != 1L || !is.numeric(groups) || !is.finite(groups) ||
      groups < 1 || groups != floor(groups) || groups > n)
    stop("groups must be an integer between 1 and the sample size.")
  if (length(threshold) != 1L || !is.numeric(threshold) || !is.finite(threshold))
    stop("threshold must be a finite number.")
  if (length(higher.y) != 1L || !is.logical(higher.y) || is.na(higher.y))
    stop("higher.y must be TRUE or FALSE.")
  if (!is.null(propensity) && (!is.numeric(propensity) || length(propensity) != n ||
      any(!is.finite(propensity)) || any(propensity <= 0 | propensity >= 1)))
    stop("propensity must contain one finite probability strictly between 0 and 1 per observation.")
  warning("validate_benefit() is experimental and under development; results may change.",
          call. = FALSE)
  weights <- if (is.null(propensity)) rep(1, n) else
    ifelse(trt == 1, 1 / propensity, 1 / (1 - propensity))
  # Author: Thomas Debray. Normalize weights separately within each arm and stratum.
  arm_mean <- function(idx, arm) {
    selected <- idx & trt == arm
    if (!any(selected)) return(NA_real_)
    sum(weights[selected] * y[selected]) / sum(weights[selected])
  }
  cuts <- unique(as.numeric(stats::quantile(benefit, probs = seq(0, 1, length.out = groups + 1L), names = FALSE)))
  group <- if (length(cuts) == 1L) rep(1L, n) else
    as.integer(cut(benefit, breaks = cuts, include.lowest = TRUE, labels = FALSE))
  group <- match(group, sort(unique(group)))
  calibration <- do.call(rbind, lapply(sort(unique(group)), function(g) {
    idx <- group == g
    y0 <- y[idx & trt == 0]; y1 <- y[idx & trt == 1]
    data.frame(group = g, n = sum(idx), n.control = length(y0),
               n.treatment = length(y1), predicted = mean(benefit[idx]),
               observed = arm_mean(idx, 1) - arm_mean(idx, 0),
               se = if (is.null(propensity) && length(y0) > 1 && length(y1) > 1)
                 sqrt(stats::var(y0) / length(y0) + stats::var(y1) / length(y1)) else NA_real_)
  }))
  fit <- c(intercept = NA_real_, slope = NA_real_, r.squared = NA_real_, rmse = NA_real_)
  if (all(is.finite(calibration$observed))) {
    fit["rmse"] <- sqrt(mean((calibration$observed - calibration$predicted)^2))
    if (nrow(calibration) >= 2L && length(unique(calibration$predicted)) > 1L) {
      model <- stats::lm(observed ~ predicted, data = calibration)
      fit[c("intercept", "slope")] <- stats::coef(model)
      if (length(unique(calibration$observed)) > 1L)
        fit["r.squared"] <- summary(model)$r.squared
    }
  }
  recommend <- if (higher.y) benefit > threshold else benefit < threshold
  # Author: Thomas Debray. Standardize only over rule strata that exist.
  value <- function(rule) {
    result <- 0
    for (s in c(FALSE, TRUE)) {
      idx <- recommend == s
      if (!any(idx)) next
      arm <- rule(s)
      outcome <- arm_mean(idx, arm)
      if (!is.finite(outcome)) return(NA_real_)
      result <- result + mean(idx) * outcome
    }
    result
  }
  rule <- value(function(s) as.integer(s))
  opposite <- value(function(s) as.integer(!s))
  control <- value(function(s) 0L)
  treatment <- value(function(s) 1L)
  direction <- if (higher.y) 1 else -1
  list(mean.bias = arm_mean(rep(TRUE, n), 1) - arm_mean(rep(TRUE, n), 0) - mean(benefit),
       calibration = calibration, calibration.fit = fit,
       policy = c(proportion.treated = mean(recommend),
                  versus.opposite = direction * (rule - opposite),
                  versus.control = direction * (rule - control),
                  versus.treatment = direction * (rule - treatment)))
}
