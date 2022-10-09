#' Doubly robust estimator of and inference for the average treatment effect for survival data
#'
#' Doubly robust estimator of the average treatment effect between two treatments, which is the
#' restricted mean time lost ratio for survival outcomes. Bootstrap is used for inference.
#'
#'
#' @param cate.model A formula describing the outcome model to be fitted.
#' The outcome must appear on the left-hand side. For survival outcomes, a \code{Surv} object
#' must be used to describe the outcome.
#' @param ps.model A formula describing the propensity score (PS) model to be fitted. The treatment must
#' appear on the left-hand side. The treatment must be a numeric vector coded as 0/1.
#' If data are from a randomized controlled trial, specify \code{ps.model = ~1} as an intercept-only model.
#' @param data A data frame containing the variables in the outcome, propensity score, and inverse
#' probability of censoring models (if specified); a data frame with \code{n} rows (1 row per observation).
#' @param ipcw.model A formula describing the inverse probability of censoring weighting (IPCW)
#' model to be fitted. The left-hand side must be empty. Only applies for survival outcomes.
#' Default is \code{NULL}, which corresponds to specifying the IPCW with the same covariates
#' as the outcome model \code{cate.model}, plus the treatment.
#' @param followup.time A column name in \code{data} specifying the maximum follow-up time,
#' interpreted as the potential censoring time. Only applies for survival outcomes.
#' Default is \code{NULL}, which corresponds to unknown potential censoring time.
#' @param tau0 The truncation time for defining restricted mean time lost. Only applies for
#' survival outcomes. Default is \code{NULL}, which corresponds to setting the truncation time as the
#' maximum survival time in the data.
#' @param surv.min Lower truncation limit for the probability of being censored.
#' It must be a positive value and should be chosen close to 0. Only applies for survival
#' outcomes. Default is \code{0.025}.
#' @param ipcw.method A character value for the censoring model. Only applies for survival
#' outcomes. Allowed values are: \code{'breslow'} (Cox regression with Breslow estimator of t
#' he baseline survivor function), \code{'aft (exponential)'}, \code{'aft (weibull)'},
#' \code{'aft (lognormal)'} or \code{'aft (loglogistic)'} (accelerated failure time model
#' with different distributions for y variable). Default is \code{'breslow'}.
#' @param ps.method A character value for the method to estimate the propensity score.
#' Allowed values include one of:
#' \code{'glm'} for logistic regression with main effects only (default), or
#' \code{'lasso'} for a logistic regression with main effects and LASSO penalization on
#' two-way interactions (added to the model if interactions are not specified in \code{ps.model}).
#' Relevant only when \code{ps.model} has more than one variable.
#' @param minPS A numerical value (in [0, 1]) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in (0, 1]) above which estimated propensity scores should be
#' truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.
#' @param n.boot A numeric value indicating the number of bootstrap samples used. Default is \code{500}.
#' @param seed An optional integer specifying an initial randomization seed for reproducibility.
#' Default is \code{NULL}, corresponding to no seed.
#' @param verbose An integer value indicating whether intermediate progress messages should
#' be printed. \code{1} indicates messages are printed and \code{0} otherwise. Default is \code{0}.
#'
#' @return Return an object of class \code{atefit} with 6 elements:
#' \itemize{
#'   \item{\code{rmst1}: } A vector of numeric values of the estimated RMST, bootstrap standard error,
#'   lower and upper limits of 95\% confidence interval, and the p-value in the group \code{trt=1}.
#'   \item{\code{rmst0}: } A vector of numeric values of the estimated RMST, bootstrap standard error,
#'   lower and upper limits of 95\% confidence interval, and the p-value in the group \code{trt=0}.
#'   \item{\code{log.rmtl.ratio}: } A vector of numeric values of the estimated log RMTL ratio of
#'   \code{trt=1} over \code{trt=0}, bootstrap standard error, lower and upper limits of 95\% confidence
#'   interval, and the p-value.
#'   \item{\code{log.hazard.ratio}: } A vector of numeric values of the estimated adjusted log hazard ratio
#'   of \code{trt=1} over \code{trt=0}, bootstrap standard error, lower and upper limits of 95\% confidence
#'   interval, and the p-value.
#'   \item{\code{trt.boot}: } Estimates of \code{rmst1}, \code{rmst0},
#'   \code{log.rmtl.ratio} and \code{log.hazard.ratio} in each bootstrap sample.
#'   \item{\code{warning}: } A warning message produced if the treatment variable was not coded as 0/1.
#'   The key to map the original coding of the variable to a 0/1 key is displayed in the warning to facilitate
#'   the interpretation of the remaining of the output.
#' }
#'
#' @details This helper function estimates the average treatment effect (ATE) for survival data between two
#' treatment groups in a given dataset. The ATE is estimated with a doubly robust estimator that accounts for
#' imbalances in covariate distributions between the two treatment groups with inverse probability treatment and
#' censoring weighting. For survival outcomes, the estimated ATE is the estimated by RMTL ratio between treatment
#' 1 versus treatment 0. The log-transformed ATEs and log-transformed adjusted hazard ratios are returned, as well
#' as the estimated RMST in either treatment group. The variability of the estimated RMTL ratio is calculated
#' using bootstrap. Additional outputs include standard error of the log RMTL ratio, 95\% confidence interval,
#' p-value, and a histogram of the bootstrap estimates.
#'
#' @examples
#' library(survival)
#' tau0 <- with(survivalExample,
#'              min(quantile(y[trt == "drug1"], 0.95), quantile(y[trt == "drug0"], 0.95)))
#'
#' output <- atefitsurv(data = survivalExample,
#'                      cate.model = Surv(y, d) ~ age + female +
#'                                   previous_cost + previous_number_relapses,
#'                      ps.model = trt ~ age + previous_treatment,
#'                      tau0 = tau0,
#'                      n.boot = 50,
#'                      seed = 999,
#'                      verbose = 1)
#' output
#' plot(output)
#'
#' @export
#'
#' @importFrom dplyr mutate
#' @importFrom tidyr gather

atefitsurv <- function(data,
                       cate.model,
                       ps.model,
                       ps.method = "glm",
                       ipcw.model = NULL,
                       ipcw.method = "breslow",
                       minPS = 0.01,
                       maxPS = 0.99,
                       followup.time = NULL,
                       tau0 = NULL,
                       surv.min = 0.025,
                       n.boot = 500,
                       seed = NULL,
                       verbose = 0) {

  # Set seed once for reproducibility
  set.seed(seed)

  #### CHECK ARGUMENTS ####
  arg.checks(fun = "drinf", response = "survival", data = data, followup.time = followup.time, tau0 = tau0, surv.min = surv.min,
             ps.method = ps.method, minPS = minPS, maxPS = maxPS, ipcw.method = ipcw.method, n.boot = n.boot, plot.boot = FALSE)

  #### PRE-PROCESSING ####
  preproc <- data.preproc.surv(fun = "drinf", cate.model = cate.model, ps.model = ps.model, ipcw.model = ipcw.model, tau0 = tau0,
                               data = data, ps.method = ps.method, response = "survival")

  y <- preproc$y
  d <- preproc$d
  trt <- preproc$trt
  x.ps <- preproc$x.ps
  x.cate <- preproc$x.cate
  x.ipcw <- preproc$x.ipcw

  # Check if tau0 is large enough, i.e., if tau0 is larger than the 50% quantile of the observed or censoring time.
  if (tau0 < median(y)) warning("It is recommended to increase tau0 close to the largest observed or censored time.")

  if (is.null(followup.time)) {
    yf <- NULL
  } else {
    yf <- data[[followup.time]]
  }

  #### FUNCTION STARTS HERE ####
  # Point estimate
  est <- drsurv(y = y, d = d, x.cate = x.cate, x.ps = x.ps, x.ipcw = x.ipcw,
                trt = trt, yf = yf,
                tau0 = tau0,
                surv.min = surv.min,
                ps.method = ps.method,
                minPS = minPS,
                maxPS = maxPS,
                ipcw.method = ipcw.method)

  est <- unlist(est)

  # Apply bootstrap
  n <- length(y)
  save.boot <- matrix(NA, nrow = n.boot, ncol = length(est))
  colnames(save.boot) <- c("rmst1", "rmst0", "log.rmtl.ratio", "log.hazard.ratio")

  pb   <- txtProgressBar(min = 1,
                         max = n.boot,
                         style = 3)

  for (i in seq(n.boot)) {
    idsub.boot <- sample(n, size = n, replace  = TRUE)
    est.boot <- drsurv(y = y[idsub.boot],
                       d = d[idsub.boot],
                       x.cate = x.cate[idsub.boot, ],
                       x.ps = x.ps[idsub.boot, ],
                       x.ipcw = x.ipcw[idsub.boot, ],
                       trt = trt[idsub.boot],
                       yf = yf[idsub.boot],
                       tau0 = tau0, surv.min = surv.min, ps.method = ps.method, minPS = minPS, maxPS = maxPS, ipcw.method = ipcw.method)
    save.boot[i,] <- unlist(est.boot)

    if (verbose == 1) setTxtProgressBar(pb, i)
  }

  close(pb)

  se.est <- apply(save.boot, 2, sd)
  cil.est <- est - qnorm(0.975) * se.est
  ciu.est <- est + qnorm(0.975) * se.est
  p.est  <- 1 - pchisq(est^2 / se.est^2, df = 1)

  est.all <- data.frame(estimate = est,
                        SE = se.est,
                        CI.lower = cil.est,
                        CI.upper = ciu.est,
                        pvalue = p.est)

  out <- c()
  out$response <- "survival"
  out$rmst1 <- data.frame(est.all[1, ])
  out$rmst0 <- data.frame(est.all[2, ])
  out$log.rmtl.ratio <- data.frame(est.all[3, ])
  out$log.hazard.ratio <- data.frame(est.all[4, ])
  out$n.boot <- n.boot # Number of  bootstrap samples
  out$trt.boot <- save.boot
  out$warning <- preproc$warning

  class(out) <- "atefit"

  return(out)
}

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

#' Probability of being censored
#'
#' Probability of being censored which is used to correct the effect of right censoring.
#'
#' @param y Observed survival or censoring time; vector of size \code{n}.
#' @param d The event indicator, normally \code{1 = event, 0 = censored}; vector of size \code{n}.
#' @param x.ipcw Matrix of \code{p.ipw} baseline covariate specified in inverse probability of censoring weighting; dimension \code{n} by \code{p.ipw}.
#' @param yf Follow-up time, interpreted as the potential censoring time; vector of size \code{n} if the potential censoring time is known.
#' If unknown, set \code{yf == NULL} and \code{yf} will be taken as \code{y} in the function.
#' @param ipcw.method The censoring model. Allowed values are: \code{'breslow'} (Cox regression with Breslow estimator of the baseline survivor function),
#' \code{'aft (exponential)'}, \code{'aft (weibull)'}, \code{'aft (lognormal)'} or \code{'aft (loglogistic)'}. Default is \code{'breslow'}.
#' @param tau0 The truncation time for defining restricted mean time lost.
#' @param surv.min Lower truncation limit for probability of being censored (positive and very close to 0).
#'
#' @return A vector of size \code{n} with the estimated probabilities \code{Pr(C > min(y, tau0) | x.ipcw)}
#' @importFrom stats pnorm
#'
ipcw.surv <- function(y, d, x.ipcw, yf = NULL, ipcw.method = "breslow", tau0, surv.min = 0.025) {

  nf <- length(yf)
  if (nf > 0) {
    dfnew <- rep(1, nf)
  } else {
    yf <- y
    dfnew <- 1 - d
  }
  yfnew <- pmin(yf, tau0)
  dfnew[yfnew == tau0] <- 0

  x <- as.matrix(x.ipcw[, apply(x.ipcw, 2, var) > 0])
  xnew <- t(t(x) - colMeans(x))

  p <- length(xnew[1, ])
  n <- length(yfnew)
  surv.prob <- rep(1, n)

  if (ipcw.method == "breslow") {
    fit.cen <- coxph(Surv(yfnew, dfnew) ~ xnew)

    beta0 <- fit.cen$coef
    detail <- coxph.detail(fit.cen)
    hazard0 <- detail$hazard
    time0 <- detail$time
    score <- exp(as.vector(xnew[, is.na(beta0) == FALSE, drop = FALSE] %*% beta0[is.na(beta0) == FALSE]))

    for (i in seq(n)) {
      surv.prob[i] <- exp(-sum(hazard0[time0 <= min(tau0, y[i])]) * score[i])
    }
  } else {
    distribution <- str_extract(string = ipcw.method, pattern = "(?<=aft \\().*(?=\\))")
    fit.cen <- survreg(Surv(yfnew, dfnew) ~ xnew, dist = distribution)

    betas <- as.numeric(fit.cen$coefficients)
    xint <- as.matrix(cbind(1, xnew))
    linpred <- as.vector(xint[, is.na(betas) == FALSE, drop = FALSE] %*% betas[is.na(betas) == FALSE])
    score <- (log(yfnew) - linpred)/fit.cen$scale

    if (distribution %in% c("exponential", "weibull")) {
      surv.prob <- exp(-exp(score))
    } else if (distribution == "lognormal") {
      surv.prob <- 1 - pnorm(score)
    } else if (distribution == "loglogistic") {
      surv.prob <- 1 - exp(score) / (1 + exp(score))
    }
  }

  surv.prob <- pmax(surv.min, surv.prob)
  return(surv.prob)
}
