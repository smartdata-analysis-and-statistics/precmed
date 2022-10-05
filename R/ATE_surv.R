#' Doubly robust estimator of the average treatment effect with Cox model for survival data
#'
#'
#' Doubly robust estimator of the average treatment effect between two treatments, which is the restricted mean time lost (RMTL) ratio
#' of treatment 1 over treatment 0 for survival outcomes.
#'
#' @param y Observed survival or censoring time; vector of size \code{n}.
#' @param d The event indicator, normally \code{1 = event, 0 = censored}; vector of size \code{n}.
#' @param x.cate Matrix of \code{p.cate} baseline covariates specified in the outcome model; dimension \code{n} by \code{p.cate}.
#' @param x.ps Matrix of \code{p.ps} baseline covariates specified in the propensity score model; dimension \code{n} by \code{p.ps}.
#' @param x.ipcw Matrix of \code{p.ipw} baseline covariate specified in inverse probability of censoring weighting; dimension \code{n} by \code{p.ipw}.
#' @param trt Treatment received; vector of size \code{n} with treatment coded as 0/1.
#' @param yf Follow-up time, interpreted as the potential censoring time; vector of size \code{n} if the potential censoring time is known.
#' @param tau0 The truncation time for defining restricted mean time lost.
#' @param surv.min Lower truncation limit for probability of being censored (positive and very close to 0).
#' @param ps.method A character value for the method to estimate the propensity score. Allowed values include one of:
#' \code{'glm'} for logistic regression with main effects only (default), or
#' \code{'lasso'} for a logistic regression with main effects and LASSO penalization on
#' two-way interactions (added to the model if interactions are not specified in \code{ps.model}).
#' Relevant only when \code{ps.model} has more than one variable.
#' @param minPS A numerical value (in [0, 1]) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in (0, 1]) above which estimated propensity scores should be
#' truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.
#' @param ipcw.method The censoring model. Allowed values are: \code{'breslow'} (Cox regression with Breslow estimator of the baseline survivor function),
#' \code{'aft (exponential)'}, \code{'aft (weibull)'}, \code{'aft (lognormal)'} or \code{'aft (loglogistic)'}. Default is \code{'breslow'}.
#'
#' @return Return a list of 4 elements:
#' \itemize{
#'   \item{\code{rmst1}: } A numeric value of the estimated restricted mean survival time n the group \code{trt = 1}.
#'   \item{\code{rmst0}: } A numeric value of the estimated restricted mean survival time n the group \code{trt = 0}.
#'   \item{\code{log.rmtl.ratio}: } A numeric value of the estimated log rmtl ratio.
#'   \item{\code{log.hazard.ratio}: } A numeric value of the estimated log hazard ratio.
#' }
#' @importFrom survival Surv coxph coxph.detail survreg

drsurv <- function(y, d, x.cate, x.ps, x.ipcw, trt, yf = NULL, tau0, surv.min = 0.025,
                   ps.method = "glm",
                   minPS = 0.01, maxPS = 0.99, ipcw.method = "breslow") {

  if (ps.method == "lasso") {
    ps <- glm.ps(trt = trt, x.ps = x.ps, xnew = NULL, minPS = minPS, maxPS = maxPS)
  } else {
    ps <- glm.simplereg.ps(trt = trt, x.ps = x.ps, xnew = NULL, minPS = minPS, maxPS = maxPS)
  }

  ynew <- pmin(y, tau0)
  dnew <- d
  dnew[ynew == tau0] <- 1

  fit <- coxph(Surv(ynew, dnew) ~ trt + x.cate)

  surv.new <- ipcw.surv(y = y, d = d, x.ipcw = x.ipcw, yf = yf, ipcw.method = ipcw.method, tau0 = tau0, surv.min = surv.min)

  weightc <- dnew / surv.new
  weightc <- weightc / mean(weightc)

  coxest1 <- cox.rmst(y = y[trt == 1], d = d[trt == 1], x.cate = as.matrix(x.cate[trt == 1, ]), xnew = x.cate, tau0 = tau0)

  coxest0 <- cox.rmst(y = y[trt == 0], d = d[trt == 0], x.cate = as.matrix(x.cate[trt == 0, ]), xnew = x.cate, tau0 = tau0)

  rmst1 <- mean((ynew - coxest1) * weightc * trt / ps) / mean(trt / ps) + mean(coxest1)
  rmst0 <- mean((ynew - coxest0) * weightc * (1 - trt) / (1 - ps)) / mean((1 - trt) / (1 - ps)) + mean(coxest0)

  if ((rmst1 > tau0 || rmst1 < 0) & (!is.na(rmst1))) rmst1 <- mean(coxest1)
  if ((rmst0 > tau0 || rmst0 < 0) & (!is.na(rmst0))) rmst0 <- mean(coxest0)

  delta.rmtl <- log((tau0 - rmst1) / (tau0 - rmst0))
  delta.hr <- as.vector(fit$coef[1])

  return(list(rmst1 = rmst1, rmst0 = rmst0, log.rmtl.ratio = delta.rmtl, log.hazard.ratio = delta.hr))
}


#' Estimate the ATE of the RMTL ratio and unadjusted hazard ratio in multiple bi-level subgroups defined by the proportions
#'
#' If only care about the higher subgroup (above cutoff), only need ate.rmtl.high and hr.high so set "onlyhigh" to be TRUE
#' Scores are adjusted to the opposite sign if \code{higher.y} == FALSE; scores stay the same if \code{higher.y} == FALSE;
#'  this is because estsurv() function always takes the subgroup of the top highest adjusted scores,
#'  and higher adjusted scores should always represent high responders of trt=1
#'
#' @param y Observed survival or censoring time; vector of size \code{n}.
#' @param d The event indicator, normally \code{1 = event, 0 = censored}; vector of size \code{n}.
#' @param x.cate Matrix of \code{p.cate} baseline covariates specified in the outcome model; dimension \code{n} by \code{p.cate}.
#' @param x.ps Matrix of \code{p.ps} baseline covariates specified in the propensity score model; dimension \code{n} by \code{p.ps}.
#' @param x.ipcw Matrix of \code{p.ipw} baseline covariate specified in inverse probability of censoring weighting; dimension \code{n} by \code{p.ipw}.
#' @param trt Treatment received; vector of size \code{n} with treatment coded as 0/1.
#' @param yf Follow-up time, interpreted as the potential censoring time; vector of size \code{n} if the potential censoring time is known.
#' @param tau0 The truncation time for defining restricted mean time lost.
#' @param score Estimated log CATE scores for all \code{n} observations from one of the five methods
#' (random forest, boosting, naive Poisson, two regressions, contrast regression); vector of size \code{n}.
#' @param higher.y A logical value indicating whether higher (\code{TRUE}) or lower (\code{FALSE})
#' @param prop Proportions corresponding to percentiles in the estimated log CATE scores that define subgroups to calculate ATE for;
#' vector of floats in (0,1] (if \code{onlyhigh=TRUE}) or in (0,1) (if \code{onlyhigh=FALSE}):
#'              Each element of \code{prop} represents the high/low cutoff in each bi-level subgroup and the length of \code{prop}
#'              is number of bi-level subgroups
#' @param onlyhigh Indicator of returning only the ATEs in the higher-than-cutoff category of the bi-level subgroups; boolean.
#' @param surv.min Lower truncation limit for probability of being censored (positive and very close to 0).
#' @param ps.method A character value for the method to estimate the propensity score. Allowed values include one of:
#' \code{'glm'} for logistic regression with main effects only (default), or
#' \code{'lasso'} for a logistic regression with main effects and LASSO penalization on
#' two-way interactions (added to the model if interactions are not specified in \code{ps.model}).
#' Relevant only when \code{ps.model} has more than one variable.
#' @param minPS A numerical value (in [0, 1]) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in (0, 1]) above which estimated propensity scores should be
#' truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.
#' @param ipcw.method The censoring model. Allowed values are: \code{'breslow'} (Cox regression with Breslow estimator of the baseline survivor function),
#' \code{'aft (exponential)'}, \code{'aft (weibull)'}, \code{'aft (lognormal)'} or \code{'aft (loglogistic)'}. Default is \code{'breslow'}.
#'
#' @return ate.rmtl.high: estimated ATEs (ratio of RMTL) in the multiple bi-level subgroups that are in the higher-than-cutoff category;
#' vector of size equal to the length of prop; always returned.
#'         ate.rmtl.low: estimated ATEs (ratio of RMTL) in the multiple bi-level subgroups that are in the lower-than-cutoff category;
#'         vector of size equal to the length of prop; returned only when \code{onlyhigh = TRUE}.
#'         hr.high: unadjusted hazard ratio in the multiple bi-level subgroups that are in the higher-than-cutoff category;
#'         vector of size equal to the length of prop; always returned.
#'         hr.low: unadjusted hazard ratio in the multiple bi-level subgroups that are in the lower-than-cutoff category;
#'         vector of size equal to the length of prop; returned only when \code{onlyhigh = TRUE}


estsurv.bilevel.subgroups <- function(y, d, x.cate, x.ps, x.ipcw, trt, yf, tau0 = tau0,
                                      score, higher.y, prop, onlyhigh,
                                      surv.min = 0.025, ps.method = "glm", minPS = 0.01, maxPS = 0.99, ipcw.method = "breslow") {

  if (higher.y == TRUE) score <- -score
  cut <- quantile(score, 1 - prop)
  n.subgroup <- length(prop)
  ate.rmtl.high <- hr.high <- rep(0, n.subgroup)
  warn.high <- err.high <- vector("list", length = n.subgroup)
  names(warn.high) <- names(err.high) <- paste("prop", round(prop, 2))
  if (onlyhigh == FALSE) {
    ate.rmtl.low <- hr.low <- rep(0, n.subgroup)
    warn.low <- err.low <- vector("list", length = n.subgroup)
    names(warn.low) <- names(err.low) <- paste("prop", round(1 - prop, 2))
  }

  for (b in seq(n.subgroup)) {
    idsub1 <- (score >= cut[b])
    fit.high <- survCatch(drsurv(y = y[idsub1], d = d[idsub1], x.cate = as.matrix(x.cate[idsub1, ]), x.ps = x.ps[idsub1, ], x.ipcw = as.matrix(x.ipcw[idsub1, ]),
                                 trt = trt[idsub1], yf = yf[idsub1],
                                 tau0 = tau0, surv.min = surv.min, ps.method = ps.method, minPS = minPS, maxPS = maxPS, ipcw.method = ipcw.method))
    ate.rmtl.high[b] <- fit.high$log.rmtl.ratio
    hr.high[b] <- fit.high$log.hazard.ratio
    warn.high[b] <- list(fit.high$warnings)
    err.high[b] <- list(fit.high$errors)

    if (onlyhigh == FALSE) {
      idsub0 <- (score < cut[b])
      fit.low <- survCatch(drsurv(y = y[idsub0], d = d[idsub0], x.cate = as.matrix(x.cate[idsub0, ]), x.ps = x.ps[idsub0, ], x.ipcw = as.matrix(x.ipcw[idsub0, ]),
                                  trt = trt[idsub0], yf = yf[idsub0],
                                  tau0 = tau0, surv.min = surv.min, ps.method = ps.method, minPS = minPS, maxPS = maxPS, ipcw.method = ipcw.method))
      ate.rmtl.low[b] <- fit.low$log.rmtl.ratio
      hr.low[b] <- fit.low$log.hazard.ratio
      warn.low[b] <- list(fit.low$warnings)
      err.low[b] <- list(fit.low$errors)
    }
  }

  if (onlyhigh == TRUE) {
    return(list(ate.rmtl.high = exp(ate.rmtl.high), hr.high = exp(hr.high), warn.high = warn.high, err.high = err.high))
  } else {
    return(list(ate.rmtl.high = exp(ate.rmtl.high), ate.rmtl.low = exp(ate.rmtl.low), hr.high = exp(hr.high), hr.low = exp(hr.low),
                warn.high = warn.high, warn.low = warn.low, err.high = err.high, err.low = err.low))
  }
}

#' Estimate the ATE of the RMTL ratio and unadjusted hazard ratio in one multilevel subgroup defined by the proportions
#'
#' Scores are adjusted to the opposite sign if \code{higher.y} == FALSE; scores stay the same if \code{higher.y} == FALSE;
#'  this is because estsurv function for multilevel subgroups start from the lowest to the highest adjusted scores,
#'  and higher adjusted scores should always represent high responders of trt=1
#'
#' @param y Observed survival or censoring time; vector of size \code{n}.
#' @param d The event indicator, normally \code{1 = event, 0 = censored}; vector of size \code{n}.
#' @param x.cate Matrix of \code{p.cate} baseline covariates specified in the outcome model; dimension \code{n} by \code{p.cate}.
#' @param x.ps Matrix of \code{p.ps} baseline covariates specified in the propensity score model; dimension \code{n} by \code{p.ps}.
#' @param x.ipcw Matrix of \code{p.ipw} baseline covariate specified in inverse probability of censoring weighting; dimension \code{n} by \code{p.ipw}.
#' @param trt Treatment received; vector of size \code{n} with treatment coded as 0/1.
#' @param yf Follow-up time, interpreted as the potential censoring time; vector of size \code{n} if the potential censoring time is known.
#' @param tau0 The truncation time for defining restricted mean time lost.
#' @param score Estimated log CATE scores for all \code{n} observations from one of the five methods
#' (random forest, boosting, naive Poisson, two regressions, contrast regression); vector of size \code{n}.
#' @param higher.y A logical value indicating whether higher (\code{TRUE}) or
#' lower (\code{FALSE}) values of the outcome are more desirable. Default is \code{TRUE}.
#' @param prop Proportions corresponding to percentiles in the estimated log CATE scores that define subgroups to calculate ATE for;
#' vector of floats in [0,1] always starting with 0 and ending with 1:
#'              Each element of \code{prop} represents inclusive cutoffs in the multilevel subgroup and the length of \code{prop}
#'              is number of categories in the multilevel subgroup
#' @param surv.min Lower truncation limit for probability of being censored (positive and very close to 0).
#' @param ps.method A character value for the method to estimate the propensity score. Allowed values include one of:
#' \code{'glm'} for logistic regression with main effects only (default), or
#' \code{'lasso'} for a logistic regression with main effects and LASSO penalization on
#' two-way interactions (added to the model if interactions are not specified in \code{ps.model}).
#' Relevant only when \code{ps.model} has more than one variable.
#' @param minPS A numerical value (in [0, 1]) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in (0, 1]) above which estimated propensity scores should be
#' truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.
#' @param ipcw.method The censoring model. Allowed values are: \code{'breslow'} (Cox regression with Breslow estimator of the baseline survivor function),
#' \code{'aft (exponential)'}, \code{'aft (weibull)'}, \code{'aft (lognormal)'} or \code{'aft (loglogistic)'}. Default is \code{'breslow'}.
#'
#' @return ate.rmtl: estimated ATEs (ratio of RMTL) of all categories in the one multilevel subgroup;
#' vector of size equal to the length of categories in the multilevel subgroup.
#'         ate.hr: unadjusted hazard ratio of all categories in the one multilevel subgroup;
#'         vector of size equal to the length of categories in the multilevel subgroup.

estsurv.multilevel.subgroups <- function(y, d, x.cate, x.ps, x.ipcw, trt, yf, tau0 = tau0,
                                         score, higher.y, prop, surv.min = 0.025,
                                         ps.method = 'glm', minPS = 0.01, maxPS = 0.99, ipcw.method = "breslow") {

  if (higher.y == TRUE) score <- -score
  cut <- quantile(score, prop)
  n.subgroup <- length(cut) - 1
  ate.rmtl <- hr <- rep(NA, n.subgroup)
  warn <- err <- vector("list", length = n.subgroup)
  names(warn) <- names(err) <- paste("prop",  round(prop[-1], 2))

  for (b in 1:n.subgroup) {
    idsub <- (score >= cut[b] & score <= cut[b + 1])
    fit <- survCatch(drsurv(y = y[idsub], d = d[idsub], x.cate = as.matrix(x.cate[idsub, ]), x.ps = x.ps[idsub, ], x.ipcw = as.matrix(x.ipcw[idsub, ]),
                            trt = trt[idsub], yf = yf[idsub],
                            tau0 = tau0, surv.min = surv.min, ps.method = ps.method, minPS = minPS, maxPS = maxPS, ipcw.method = ipcw.method))
    ate.rmtl[b] <- fit$log.rmtl.ratio
    hr[b] <- fit$log.hazard.ratio
    warn[b] <- list(fit$warnings)
    err[b] <- list(fit$errors)
  }
  return(list(ate.rmtl = exp(ate.rmtl), hr = exp(hr), warn = warn, err = err))
}
