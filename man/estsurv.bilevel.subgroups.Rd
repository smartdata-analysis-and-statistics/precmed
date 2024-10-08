% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ATE_surv.R
\name{estsurv.bilevel.subgroups}
\alias{estsurv.bilevel.subgroups}
\title{Estimate the ATE of the RMTL ratio and unadjusted hazard ratio in multiple bi-level subgroups defined by the proportions}
\usage{
estsurv.bilevel.subgroups(
  y,
  d,
  x.cate,
  x.ps,
  x.ipcw,
  trt,
  yf,
  tau0 = tau0,
  score,
  higher.y,
  prop,
  onlyhigh,
  surv.min = 0.025,
  ps.method = "glm",
  minPS = 0.01,
  maxPS = 0.99,
  ipcw.method = "breslow"
)
}
\arguments{
\item{y}{Observed survival or censoring time; vector of size \code{n}.}

\item{d}{The event indicator, normally \code{1 = event, 0 = censored}; vector of size \code{n}.}

\item{x.cate}{Matrix of \code{p.cate} baseline covariates specified in the outcome model; dimension \code{n} by \code{p.cate}.}

\item{x.ps}{Matrix of \code{p.ps} baseline covariates specified in the propensity score model; dimension \code{n} by \code{p.ps}.}

\item{x.ipcw}{Matrix of \code{p.ipw} baseline covariate specified in inverse probability of censoring weighting; dimension \code{n} by \code{p.ipw}.}

\item{trt}{Treatment received; vector of size \code{n} with treatment coded as 0/1.}

\item{yf}{Follow-up time, interpreted as the potential censoring time; vector of size \code{n} if the potential censoring time is known.}

\item{tau0}{The truncation time for defining restricted mean time lost.}

\item{score}{Estimated log CATE scores for all \code{n} observations from one of the five methods
(random forest, boosting, naive Poisson, two regressions, contrast regression); vector of size \code{n}.}

\item{higher.y}{A logical value indicating whether higher (\code{TRUE}) or lower (\code{FALSE})}

\item{prop}{Proportions corresponding to percentiles in the estimated log CATE scores that define subgroups to calculate ATE for;
vector of floats in `(0, 1]` (if \code{onlyhigh=TRUE}) or in `(0, 1)` (if \code{onlyhigh=FALSE}):
             Each element of \code{prop} represents the high/low cutoff in each bi-level subgroup and the length of \code{prop}
             is number of bi-level subgroups}

\item{onlyhigh}{Indicator of returning only the ATEs in the higher-than-cutoff category of the bi-level subgroups; boolean.}

\item{surv.min}{Lower truncation limit for probability of being censored (positive and very close to 0).}

\item{ps.method}{A character value for the method to estimate the propensity score. Allowed values include one of:
\code{'glm'} for logistic regression with main effects only (default), or
\code{'lasso'} for a logistic regression with main effects and LASSO penalization on
two-way interactions (added to the model if interactions are not specified in \code{ps.model}).
Relevant only when \code{ps.model} has more than one variable.}

\item{minPS}{A numerical value (in `[0, 1]`) below which estimated propensity scores should be
truncated. Default is \code{0.01}.}

\item{maxPS}{A numerical value (in `(0, 1]`) above which estimated propensity scores should be
truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.}

\item{ipcw.method}{The censoring model. Allowed values are: \code{'breslow'} (Cox regression with Breslow estimator of the baseline survivor function),
\code{'aft (exponential)'}, \code{'aft (weibull)'}, \code{'aft (lognormal)'} or \code{'aft (loglogistic)'}. Default is \code{'breslow'}.}
}
\value{
ate.rmtl.high: estimated ATEs (ratio of RMTL) in the multiple bi-level subgroups that are in the higher-than-cutoff category;
vector of size equal to the length of prop; always returned.
        ate.rmtl.low: estimated ATEs (ratio of RMTL) in the multiple bi-level subgroups that are in the lower-than-cutoff category;
        vector of size equal to the length of prop; returned only when \code{onlyhigh = TRUE}.
        hr.high: unadjusted hazard ratio in the multiple bi-level subgroups that are in the higher-than-cutoff category;
        vector of size equal to the length of prop; always returned.
        hr.low: unadjusted hazard ratio in the multiple bi-level subgroups that are in the lower-than-cutoff category;
        vector of size equal to the length of prop; returned only when \code{onlyhigh = TRUE}
}
\description{
If only care about the higher subgroup (above cutoff), only need ate.rmtl.high and hr.high so set "onlyhigh" to be TRUE
Scores are adjusted to the opposite sign if \code{higher.y} == FALSE; scores stay the same if \code{higher.y} == FALSE;
 this is because estsurv() function always takes the subgroup of the top highest adjusted scores,
 and higher adjusted scores should always represent high responders of trt=1
}
