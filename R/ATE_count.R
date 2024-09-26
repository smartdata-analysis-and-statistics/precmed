#' Doubly robust estimator of and inference for the average treatment effect
#' for count data
#'
#' Doubly robust estimator of the average treatment effect between two
#' treatments, which is the rate ratio for count outcomes. Bootstrap is used for
#' inference.
#'
#' @param cate.model A formula describing the outcome model to be fitted.
#' The outcome must appear on the left-hand side.
#' @param ps.model A formula describing the propensity score (PS) model to be
#' fitted. The treatment must appear on the left-hand side. The treatment must
#' be a numeric vector coded as 0 or 1. If data are from a randomized controlled
#' trial, specify \code{ps.model = ~1} as an intercept-only model.
#' @param data A data frame containing the variables in the outcome, propensity
#' score, and inverse probability of censoring models (if specified); a data
#' frame with \code{n} rows (1 row per observation).
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
#' should assume treatment-covariate interaction by \code{x}. If \code{TRUE},
#' interactions will be assumed only if at least 10 patients received each
#' treatment option. Default is \code{TRUE}.
#' @param n.boot A numeric value indicating the number of bootstrap samples
#' used. Default is \code{500}.
#' @param seed An optional integer specifying an initial randomization seed for
#' reproducibility. Default is \code{NULL}, corresponding to no seed.
#' @param verbose An integer value indicating whether intermediate progress
#' messages should be printed. \code{1} indicates messages are printed and
#' \code{0} otherwise. Default is \code{0}.
#'
#' @return Return an item of the class \code{atefit} with the following
#' elements:
#' \itemize{
#'   \item{\code{log.rate.ratio}: } A vector of numeric values of the estimated
#'   ATE (expressed as a log rate ratio of \code{trt=1} over \code{trt=0}),
#'   the bootstrap standard error, the lower and upper limits of 95\% confidence
#'   interval, and the p-value.
#'   \item{\code{rate0}: } A numeric value of the estimated rate in the group
#'   \code{trt=0}.
#'   \item{\code{rate1}: } A numeric value of the estimated rate in the group
#'   \code{trt=1}.
#'   \item{\code{trt.boot}: } Estimated log rate ratios in each bootstrap
#'   sample.
#'   \item{\code{warning}: } A warning message produced if the treatment
#'   variable was not coded as 0 or 1. The key to map the original coding of the
#'   variable to a 0-1 coding is displayed in the warning to facilitate the
#'   interpretation of the remaining of the output.
#' }
#'
#' @details This helper function estimates the average treatment effect (ATE)
#' between two treatment groups in a given dataset. The ATE is estimated with a
#' doubly robust estimator that accounts for imbalances in covariate
#' distributions between the two treatment groups with inverse probability
#' treatment weighting. For count outcomes, the estimated ATE is the estimated
#' rate ratio between treatment 1 versus treatment 0.
#'
#' @examples
#' output <- atefitcount(data = countExample,
#'                       cate.model = y ~ age + female + previous_treatment +
#'                                    previous_cost + previous_number_relapses +
#'                                    offset(log(years)),
#'                       ps.model = trt ~ age + previous_treatment,
#'                       verbose = 1, n.boot = 50, seed = 999)
#' output
#' plot(output)
#' @export
#'
#' @importFrom dplyr mutate
#'

atefitcount <- function(data,
                        cate.model,
                        ps.model,
                        ps.method = "glm",
                        minPS = 0.01,
                        maxPS = 0.99,
                        interactions = TRUE,
                        n.boot = 500,
                        seed = NULL,
                        verbose = 0) {

  # Set seed once for reproducibility
  set.seed(seed)

  #### CHECK ARGUMENTS ####
  arg.checks(fun = "drinf",
             response = "count",
             data = data,
             ps.method = ps.method,
             minPS = minPS,
             maxPS = maxPS,
             interactions = interactions,
             n.boot = n.boot,
             plot.boot = FALSE)

  #### PRE-PROCESSING ####
  preproc <- data.preproc(fun = "drinf", cate.model = cate.model, ps.model = ps.model,
                          data = data, ps.method = ps.method)
  y <- preproc$y
  trt <- preproc$trt
  x.ps <- preproc$x.ps
  x.cate <- preproc$x.cate
  time <- preproc$time

  # Point estimate
  est <- drcount(y = y,
                 x.cate = x.cate,
                 x.ps = x.ps,
                 trt = trt,
                 time = time,
                 ps.method = ps.method,
                 minPS = minPS,
                 maxPS = maxPS,
                 interactions = interactions)
  logrr <- est$log.rate.ratio

  # Apply bootstrap
  n <- length(y)
  if (is.na(logrr)) {
    stop("Impossible to use bootstrap, log(RR)=NA.")
  }

  trt.boot <- rep(NA, n.boot)

  pb   <- txtProgressBar(min = 1,
                         max = n.boot,
                         style = 3)

  for (i.boot in seq(n.boot)) {
    idsub.boot <- sample(n, size = n, replace = TRUE)
    trt.boot[i.boot] <- drcount(y = y[idsub.boot],
                                x.cate = x.cate[idsub.boot, , drop = FALSE],
                                x.ps = x.ps[idsub.boot, , drop = FALSE],
                                trt = trt[idsub.boot],
                                time = time[idsub.boot],
                                ps.method = ps.method,
                                minPS = minPS,
                                maxPS = maxPS,
                                interactions = interactions)$log.rate.ratio

    if (verbose == 1) setTxtProgressBar(pb, i.boot)
  }
  close(pb)

  out <- c()
  out$response <- "count"
  se.est <- sd(trt.boot, na.rm = TRUE)
  out$log.rate.ratio <- data.frame(estimate = logrr,
                                   SE = se.est,
                                   CI.lower = logrr - qnorm(0.975) * se.est,
                                   CI.upper = logrr + qnorm(0.975) * se.est,
                                   pvalue = 1 - pchisq(logrr^2 / se.est^2, 1))
  rownames(out$log.rate.ratio) <- "log.rate.ratio"

  out$rate0 <- data.frame(estimate = est$rate0)
  rownames(out$rate0) <- "rate0"
  out$rate1 <- data.frame(estimate = est$rate1)
  rownames(out$rate1) <- "rate1"

  out$n.boot <- n.boot # Number of  bootstrap samples
  out$trt.boot <- trt.boot #bootstrap estimates of the treatment effect
  out$warning <- preproc$warning

  class(out) <- "atefit"


  return(out)
}

#' Doubly robust estimator of the average treatment effect for count data
#'
#' Doubly robust estimator of the average treatment effect between two
#' treatments, which is the rate ratio of treatment 1 over treatment 0 for
#' count outcomes.
#'
#' @param y A numeric vector of size \code{n} with each element representing
#' the observed count outcome for each subject.
#' @param trt A numeric vector (in \{0, 1\}) of size \code{n} with each element
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

  if (!higher.y) score <- -score
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
  if (onlyhigh) {
    return(exp(ate.est.high))
  }
  return(list(ate.est.high = exp(ate.est.high), ate.est.low = exp(ate.est.low)))
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
