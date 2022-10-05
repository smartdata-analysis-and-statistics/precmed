#' Doubly robust estimator of the average treatment effect for count data
#'
#' Doubly robust estimator of the average treatment effect between two
#' treatments, which is the rate ratio of treatment 1 over treatment 0 for
#' count outcomes.
#'
#' @param y A numeric vector of size \code{n} with each element representing
#' the observed count outcome for each subject.
#' @param trt A numeric vector (in {0, 1}) of size \code{n} with each element
#' representing the treatment received for each subject.
#' @param x.cate A numeric matrix of dimension \code{n} by \code{p.cate} with
#' each column representing each baseline covariate specified in the outcome
#' model for all subjects.
#' @param x.ps A numeric matrix of dimension \code{n} by \code{p.ps + 1} with
#' a leading column of 1 as the intercept and each remaining column
#' representing each baseline covariate specified in the propensity score model
#' for all subjects.
#' @param time A numeric vector of size \code{n} with each element representing
#' the log-transformed person-years of follow-up for each subject.
#' @param ps.method A character value for the method to estimate the propensity
#' score. Allowed values include one of: \code{'glm'} for logistic regression
#' with main effects only (default), or \code{'lasso'} for a logistic regression
#' with main effects and LASSO penalization on two-way interactions (added to
#' the model if interactions are not specified in \code{ps.model}). Relevant
#' only when \code{ps.model} has more than one variable.
#' @param minPS A numerical value between 0 and 1 below which estimated
#' propensity scores should be truncated. Default is \code{0.01}.
#' @param maxPS A numerical value between 0 and 1 above which estimated
#' propensity scores should be truncated. Must be strictly greater than
#' \code{minPS}. Default is \code{0.99}.
#' @param interactions A logical value indicating whether the outcome model
#' should allow for treatment-covariate interaction by \code{x}. If \code{TRUE},
#' interactions will be assumed only if at least 10 patients received each
#' treatment option. Default is \code{TRUE}.
#'
#' @return Return a list of 4 elements:
#' \itemize{
#'   \item{\code{log.rate.ratio}: } A numeric value of the estimated log rate ratio.
#'   \item{\code{rate0}: } A numeric value of the estimated rate in the group trt=0.
#'   \item{\code{rate1}: } A numeric value of the estimated rate in the group trt=1.
#' }
#'

drcount <- function(y, trt, x.cate, x.ps, time,
                    ps.method = "glm", minPS = 0.01, maxPS = 0.99,
                    interactions = TRUE) {

  n <- length(y)

  withCallingHandlers({
    if (ps.method == "lasso") {
      ps <- glm.ps(trt = trt, x.ps = x.ps, xnew = NULL, minPS = minPS, maxPS = maxPS)
    } else {
      ps <- glm.simplereg.ps(trt = trt, x.ps = x.ps, xnew = NULL, minPS = minPS, maxPS = maxPS)
    }
  },
  warning = function(w) {
    if (grepl("numerically 0|1 occurred", conditionMessage(w)) | grepl("glm.fit: algorithm did not converge", conditionMessage(w)))
      invokeRestart("muffleWarning")  # suppress warning: "glm.fit: fitted probabilities/rates numerically 0 occurred." or "glm.fit: algorithm did not converge."
  })

  if (min(sum(trt), n - sum(trt)) > 10 & interactions == TRUE) {
    withCallingHandlers({
      # Calculate rate by treatment arm by fitting two separate models
      fit1 <- glm(y ~ x.cate + offset(time), family = "poisson", subset = (trt == 1))
      fit0 <- glm(y ~ x.cate + offset(time), family = "poisson", subset = (trt == 0))
    },
    warning = function(w) {
      if (grepl("numerically 0|1 occurred|algorithm did not converge", conditionMessage(w)))
        invokeRestart("muffleWarning")  # suppress warning: "glm.fit: fitted rates numerically 0 occurred." or "glm.fit: algorithm did not converge."
    })
    beta1 <- fit1$coef
    y1hat <- exp(cbind(1, x.cate)[, is.na(beta1) == FALSE, drop = FALSE] %*% beta1[is.na(beta1) == FALSE])
    beta0 <- fit0$coef
    y0hat <- exp(cbind(1, x.cate)[, is.na(beta0) == FALSE, drop = FALSE] %*% beta0[is.na(beta0) == FALSE])
  } else {
    withCallingHandlers({
      # Calculate rate by treatment arm by fitting a single model
      fit <- glm(y ~ x.cate + trt + offset(time), family = "poisson")
    },
    warning = function(w) {
      if (grepl("numerically 0|1 occurred|algorithm did not converge", conditionMessage(w)))
        invokeRestart("muffleWarning")  # suppress warning: "glm.fit: fitted rates numerically 0 occurred." or "glm.fit: algorithm did not converge."
    })
    beta <- fit$coef
    y1hat <- exp(cbind(1, x.cate, 1)[, is.na(beta) == FALSE, drop = FALSE] %*% beta[is.na(beta) == FALSE])
    y0hat <- exp(cbind(1, x.cate, 0)[, is.na(beta) == FALSE, drop = FALSE] %*% beta[is.na(beta) == FALSE])
  }
  mu1dr <- sum(y1hat + (y / exp(time) - y1hat) * trt / ps) / n
  mu0dr <- sum(y0hat + (y / exp(time) - y0hat) * (1 - trt) / (1 - ps)) / n

  if (is.na(mu1dr / mu0dr) == FALSE & mu1dr / mu0dr > 0) {
    logrr <- log(mu1dr / mu0dr)
  } else {
    logrr <- NA
  }

  return(list(log.rate.ratio = logrr, rate0 = mu0dr, rate1 = mu1dr))
}


#' Estimate the Average Treatment Effect of the log risk ratio in multiple
#' bi-level subgroups defined by the proportions
#'
#' If only care about the higher subgroup (above cutoff), only need trt.est.high so set \code{onlyhigh} to be TRUE
#' Scores are adjusted to the opposite sign if \code{higher.y} == FALSE; scores stay the same if \code{higher.y} == TRUE;
#'  this is because estcount.bilevel.subgroups() always takes the subgroup of the top highest adjusted scores,
#'  and higher adjusted scores should always represent high responders of trt=1
#'
#' @param y Observed outcome; vector of size \code{n} (observations)
#' @param trt Treatment received; vector of size \code{n} units with treatment coded as 0/1
#' @param x.cate Matrix of \code{p.cate} baseline covariates; dimension \code{n} by \code{p.cate} (covariates in the outcome model)
#' @param x.ps Matrix of \code{p.ps} baseline covariates (plus a leading column of 1 for the intercept);
#' dimension \code{n} by \code{p.ps + 1} (covariates in the propensity score model plus intercept)
#' @param time Log-transformed person-years of follow-up; vector of size \code{n}
#' @param score Estimated log CATE scores for all \code{n} observations from one of the four methods
#' (boosting, naive Poisson, two regressions, contrast regression); vector of size \code{n}
#' @param prop Proportions corresponding to percentiles in the estimated log CATE scores that define subgroups to calculate ATE for;
#' vector of floats in (0,1] (if onlyhigh=T) or in (0,1) (if onlyhigh=F):
#'              Each element of \code{prop} represents the high/low cutoff in each bi-level subgroup and the length of \code{prop}
#'              is number of bi-level subgroups
#' @param onlyhigh Indicator of returning only the ATEs in the higher-than-cutoff category of the bi-level subgroups; boolean
#' @param higher.y A logical value indicating whether higher (\code{TRUE}) or lower (\code{FALSE})
#' values of the outcome are more desirable. Default is \code{TRUE}.
#' @param ps.method A character value for the method to estimate the propensity score. Allowed values include one of:
#' \code{'glm'} for logistic regression with main effects only (default), or
#' \code{'lasso'} for a logistic regression with main effects and LASSO penalization on
#' two-way interactions (added to the model if interactions are not specified in \code{ps.model}).
#' Relevant only when \code{ps.model} has more than one variable.
#' @param minPS A numerical value (in [0, 1]) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in (0, 1]) above which estimated propensity scores should be
#' truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.
#'
#' @return ate.est.high: estimated ATEs in the multiple bi-level subgroups that are in the higher-than-cutoff category;
#' vector of size equal to the length of prop; always returned
#'         ate.est.low: estimated ATEs in the multiple bi-level subgroups that are in the lower-than-cutoff category;
#'         vector of size equal to the length of prop; returned only when \code{onlyhigh} == TRUE
estcount.bilevel.subgroups <- function(y, x.cate, x.ps, time, trt, score, higher.y,
                                       prop, onlyhigh,
                                       ps.method = "glm", minPS = 0.01, maxPS = 0.99) {

  if (higher.y == FALSE) score <- -score
  cut <- quantile(score, 1 - prop)
  n.subgroup <- length(prop)
  ate.est.high <- rep(0, n.subgroup)
  if (onlyhigh == FALSE) ate.est.low <- rep(0, n.subgroup)
  for (b in 1:n.subgroup) {
      idsub1 <- (score >= cut[b])
      ate.est.high[b] <- drcount(y = y[idsub1], x.cate = x.cate[idsub1, , drop = FALSE], x.ps = x.ps[idsub1, , drop = FALSE],
                                 trt = trt[idsub1], time = time[idsub1], ps.method = ps.method, minPS = minPS, maxPS = maxPS,
                                 interactions = FALSE)$log.rate.ratio
      if (onlyhigh == FALSE) {
        idsub0 <- (score < cut[b])
        ate.est.low[b] <- drcount(y = y[idsub0], x.cate = x.cate[idsub0, , drop = FALSE], x.ps = x.ps[idsub0, , drop = FALSE],
                                  trt = trt[idsub0], time = time[idsub0], ps.method = ps.method, minPS = minPS, maxPS = maxPS,
                                  interactions = FALSE)$log.rate.ratio
      }
  }
  if (onlyhigh == TRUE) {
    return(exp(ate.est.high))
  } else {
    return(list(ate.est.high = exp(ate.est.high), ate.est.low = exp(ate.est.low)))
  }
}


#' Estimate the ATE of the log RR ratio in one multilevel subgroup defined by the proportions
#'
#' Scores are adjusted to the opposite sign if \code{higher.y} == FALSE; scores stay the same if \code{higher.y} == TRUE;
#'  this is because subgroups defined in estcount.multilevel.subgroup() start from the lowest to the highest adjusted scores,
#'  and higher adjusted scores should always represent high responders of trt=1
#'
#' @param y Observed outcome; vector of size \code{n} (observations)
#' @param trt Treatment received; vector of size \code{n} units with treatment coded as 0/1
#' @param x.cate Matrix of \code{p.cate} baseline covariates; dimension \code{n} by \code{p.cate} (covariates in the outcome model)
#' @param x.ps Matrix of \code{p.ps} baseline covariates (plus a leading column of 1 for the intercept);
#' dimension \code{n} by \code{p.ps + 1} (covariates in the propensity score model plus intercept)
#' @param time Log-transformed person-years of follow-up; vector of size \code{n}
#' @param score Estimated log CATE scores for all \code{n} observations from one of the four methods
#' (boosting, naive Poisson, two regressions, contrast regression); vector of size \code{n}
#' @param prop Proportions corresponding to percentiles in the estimated log CATE scores that define subgroups to calculate ATE for;
#' vector of floats in [0,1] always starting with 0 and ending with 1:
#'              Each element of \code{prop} represents inclusive cutoffs in the multilevel subgroup and the length of \code{prop}
#'              is number of categories in the multilevel subgroup
#' @param higher.y A logical value indicating whether higher (\code{TRUE}) or lower (\code{FALSE})
#' values of the outcome are more desirable. Default is \code{TRUE}.
#' @param ps.method A character value for the method to estimate the propensity score. Allowed values include one of:
#' \code{'glm'} for logistic regression with main effects only (default), or
#' \code{'lasso'} for a logistic regression with main effects and LASSO penalization on
#' two-way interactions (added to the model if interactions are not specified in \code{ps.model}).
#' Relevant only when \code{ps.model} has more than one variable.
#' @param minPS A numerical value (in [0, 1]) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in (0, 1]) above which estimated propensity scores should be
#' truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.
#'
#' @return estimated ATEs of all categories in the one multilevel subgroup; vector of size equal to the length of categories in the multilevel subgroup
estcount.multilevel.subgroup <- function(y, x.cate, x.ps, time, trt, score, higher.y, prop,
                                         ps.method = "glm", minPS = 0.01, maxPS = 0.99){

  if (higher.y == FALSE) score <- -score
  cut <- quantile(score, prop)
  n.subgroup <- length(prop) - 1
  ate.est <- rep(NA, n.subgroup)
  for (b in seq(n.subgroup)) {
    idsub <- (score >= cut[b] & score <= cut[b + 1])
    ate.est[b] <- drcount(y = y[idsub],
                          x.cate = x.cate[idsub, , drop = FALSE],
                          x.ps = x.ps[idsub, , drop = FALSE],
                          trt = trt[idsub], time = time[idsub],
                          ps.method = ps.method,
                          minPS = minPS,
                          maxPS = maxPS,
                          interactions = FALSE)$log.rate.ratio
  }
  return(exp(ate.est))
}
