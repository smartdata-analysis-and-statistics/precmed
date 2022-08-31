# ------------------------------------------------------------------
#
# Project: Precision Medicine MS - Comprehensive R package
#
# Purpose: Utility functions for Survival outcomes
#
# Platform: Windows
# R Version: 4.1.0
#


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

    for(i in 1:n) {
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



#' Estimate restricted mean survival time (RMST) based on Cox regression model
#'
#'
#' @param y Observed survival or censoring time; vector of size \code{n}.
#' @param d The event indicator, normally \code{1 = event, 0 = censored}; vector of size \code{n}.
#' @param x.cate Matrix of \code{p.cate} baseline covariates specified in the outcome model; dimension \code{n} by \code{p.cate}.
#' @param xnew Matrix of \code{p.cate} baseline covariates for which we want an estimate of the RMST; dimension \code{m} (observations in the new data set) by \code{p.cate}
#' @param tau0 The truncation time for defining restricted mean time lost.
#'
#' @return The estimated RMST for new subjects with covariates \code{xnew}; vector of size \code{m}.
#'
cox.rmst <- function(y, d, x.cate, xnew, tau0) {
  xnew <- as.matrix(xnew)
  x <- as.matrix(x.cate)

  vx <- apply(x, 2, var)
  xnew <- xnew[, vx > 0, drop = FALSE]
  x <- x[, vx > 0, drop = FALSE]

  mux <- colMeans(x)
  xnew <- t(t(xnew) - mux)
  x <- t(t(x) - mux)

  fit <- coxph(Surv(y, d) ~ x)
  beta <- fit$coef
  detail <- coxph.detail(fit)
  time <- detail$time
  surv0 <- exp(-cumsum(detail$hazard))

  n <- length(xnew[, 1])
  surv0 <- surv0[time <= tau0]
  time <- time[time <= tau0]
  hr <- exp(as.vector(xnew[, is.na(beta) == FALSE, drop = FALSE] %*% beta[is.na(beta) == FALSE, drop = FALSE]))
  surv <- crossprod(t(rep(1, n)), surv0)^hr

  gap <- c(time, tau0) - c(0, time)
  coxest <- colSums(t(cbind(1, surv)) * gap)

  return(coxest)
}




#' Split the given time-to-event dataset into balanced training and validation sets (within a pre-specified tolerance)
#' Balanced means 1) The ratio of treated and controls is maintained in the training and validation sets
#'                2) The covariate distributions are balanced between the training and validation sets
#'
#' @param y Observed survival or censoring time; vector of size \code{n}.
#' @param d The event indicator, normally \code{1 = event, 0 = censored}; vector of size \code{n}.
#' @param trt Treatment received; vector of size \code{n} with treatment coded as 0/1.
#' @param x.cate Matrix of \code{p.cate} baseline covariates specified in the outcome model; dimension \code{n} by \code{p.cate}.
#' @param x.ps Matrix of \code{p.ps} baseline covariates specified in the propensity score model; dimension \code{n} by \code{p.ps}.
#' @param x.ipcw Matrix of \code{p.ipw} baseline covariate specified in inverse probability of censoring weighting; dimension \code{n} by \code{p.ipw}.
#' @param yf Follow-up time, interpreted as the potential censoring time; vector of size \code{n} if the potential censoring time is known.
#' If unknown, set \code{yf == NULL} and \code{yf} will be taken as \code{y} in the function.
#' @param train.prop A numerical value (in (0, 1)) indicating the proportion of total data used
#' for training. Default is \code{3/4}.
#' @param error.max A numerical value > 0 indicating the tolerance (maximum value of error)
#' for the largest standardized absolute difference in the covariate distributions or in the
#' doubly robust estimated rate ratios between the training and validation sets. This is used
#' to define a balanced training-validation splitting. Default is \code{0.1}.
#' @param max.iter A positive integer value indicating the maximum number of iterations when
#' searching for a balanced training-validation split. Default is \code{5,000}.
#'
#' @return A list of 14 objects, 7training and 7 validation of y, trt, x.cate, x.ps, x.ipcw, time, yf:
#'             y.train          - observed survival or censoring time in the training set; vector of size \code{m} (observations in the training set)
#'             d.train          - event indicator in the training set; vector of size \code{m} coded as 0/1
#'             trt.train        - treatment received in the training set; vector of size \code{m} coded as 0/1
#'             x.cate.train     - baseline covariates for the outcome model in the training set; matrix of dimension \code{m} by \code{p.cate}
#'             x.ps.train       - baseline covariates (plus intercept) for the propensity score model in the training set; matrix of dimension \code{m} by \code{p.ps + 1}
#'             x.ipcw.train      - baseline covariates for inverse probability of censoring in the training set; matrix of dimension \code{m} by \code{p.ipw}
#'             yf.train         - follow-up time in the training set; if known, vector of size \code{m}; if unknown, \code{yf == NULL}
#'             y.valid          - observed survival or censoring time in the validation set; vector of size \code{n-m}
#'             d.valid          - event indicator in the validation set; vector of size \code{n-m} coded as 0/1
#'             trt.valid        - treatment received in the validation set; vector of size \code{n-m} coded as 0/1
#'             x.cate.valid     - baseline covariates for the outcome model in the validation set; matrix of dimension \code{n-m} by \code{p.cate}
#'             x.ps.valid       - baseline covariates (plus intercept) for the propensity score model in the validation set; matrix of dimension \code{n-m} by \code{p.ps + 1}
#'             x.ipcw.valid      - baseline covariates for inverse probability of censoring in the validation set; matrix of dimension \code{n-m} by \code{p.ipw}
#'             yf.valid         - follow-up time in the training set; if known, vector of size \code{n-m}; if unknown, \code{yf == NULL}
#'

balancesurv.split <- function(y, d, trt, x.cate, x.ps, x.ipcw, yf = NULL, train.prop = 3/4, error.max = 0.1, max.iter = 5000) {

  x <- cbind(x.cate, x.ps[, -1], x.ipcw)
  x <- x[, !duplicated(colnames(x)), drop = FALSE]
  sdx <- apply(x, 2, sd)

  n1 <- sum(trt)
  n0 <- sum(1 - trt)
  m1 <- round(n1 * train.prop, digits = 0)
  m0 <- round(n0 * train.prop, digits = 0)
  n <- n1 + n0
  m <- m1 + m0
  id1 <- (1:n)[trt == 1]
  id0 <- (1:n)[trt == 0]
  error <- errorbest <- Inf
  bestid.valid <- rep(NA, n - m)
  iter <- 0

  while ((error > error.max | is.na(error) == TRUE) & iter <= max.iter) {
    id.valid <- c(sample(id1, n1 - m1, replace = FALSE), sample(id0, n0 - m0, replace = FALSE))

    y.valid <- y[id.valid]
    d.valid <- d[id.valid]
    trt.valid <- trt[id.valid]
    x.valid  <- x[id.valid, , drop = FALSE]
    x.cate.valid <- x.cate[id.valid, , drop = FALSE]
    x.ps.valid <- x.ps[id.valid, , drop = FALSE]
    x.ipcw.valid <- x.ipcw[id.valid, , drop = FALSE]
    yf.valid <- yf[id.valid]

    y.train <- y[-id.valid]
    d.train <- d[-id.valid]
    trt.train <- trt[-id.valid]
    x.train <- x[-id.valid, , drop = FALSE]
    x.cate.train <- x.cate[-id.valid, , drop = FALSE]
    x.ps.train <- x.ps[-id.valid, , drop = FALSE]
    x.ipcw.train <- x.ipcw[-id.valid, , drop = FALSE]
    yf.train <- yf[-id.valid]

    diffx <- colMeans(x.train) - colMeans(x.valid)
    errorx <- max(abs(diffx / sdx))

    est.valid <- coxph(Surv(y.valid, d.valid) ~ trt.valid + x.valid)$coef[1]
    est.train <- coxph(Surv(y.train, d.train) ~ trt.train + x.train)$coef[1]

    error <- max(2 * abs(est.valid - est.train) / (abs(est.valid) + abs(est.train)), errorx)
    if (is.na(error) == FALSE & error < errorbest) {
      bestid.valid <- id.valid
      errorbest <- error
    }
    iter <- iter + 1
  }

  if (iter == max.iter + 1) {
    y.valid <- y[bestid.valid]
    d.valid <- d[bestid.valid]
    trt.valid <- trt[bestid.valid]
    x.cate.valid <- x.cate[bestid.valid, , drop = FALSE]
    x.ps.valid <- x.ps[bestid.valid, , drop = FALSE]
    x.ipcw.valid <- x.ipcw[bestid.valid, , drop = FALSE]

    y.train <- y[-bestid.valid]
    d.train <- d[-bestid.valid]
    trt.train <- trt[-bestid.valid]
    x.cate.train <- x.cate[-bestid.valid, , drop = FALSE]
    x.ps.train <- x.ps[-bestid.valid, , drop = FALSE]
    x.ipcw.train <- x.ipcw[-bestid.valid, , drop = FALSE]
    warning(paste("Maximum iteration reached and the SMD between training and validation set is still greater than error.max (error=", round(errorbest, 4), "). Consider increasing max.iter, decreasing error.max, or increasing sample size.", sep = ""))
  }

  return(list(y.train = y.train, d.train = d.train, trt.train = trt.train, x.cate.train = x.cate.train,
              x.ps.train = x.ps.train, x.ipcw.train = x.ipcw.train, yf.train = yf.train,
              y.valid = y.valid, d.valid = d.valid, trt.valid = trt.valid, x.cate.valid = x.cate.valid,
              x.ps.valid = x.ps.valid, x.ipcw.valid = x.ipcw.valid, yf.valid = yf.valid))
}


#' Check arguments that are common to all types of outcome
#' USed inside \code{arg.checks()}
#'
#' @param fun A function for which argument check is needed; "pm" for \code{pmcount()}, "cv" for \code{cvcount()},
#' and "drinf" for \code{drcount.inference()}. No default.
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
#' @param higher.y A logical value indicating whether higher (\code{TRUE}) or lower (\code{FALSE})
#' values of the outcome are more desirable. Default is \code{TRUE}.
#' @param abc A logical value indicating whether the area between curves (ABC) should be calculated
#' at each cross-validation iterations, for each \code{score.method}. Default is \code{TRUE}.
#' @param prop.cutoff A vector of numerical values (in (0, 1]) specifying percentiles of the
#' estimated log CATE scores to define nested subgroups. Each element represents the cutoff to
#' separate observations in nested subgroups (below vs above cutoff).
#' The length of \code{prop.cutoff} is the number of nested subgroups.
#' An equally-spaced sequence of proportions ending with 1 is recommended.
#' Default is \code{seq(0.5, 1, length = 6)}.
#' @param prop.multi A vector of numerical values (in [0, 1]) specifying percentiles of the
#' estimated log CATE scores to define mutually exclusive subgroups.
#' It should start with 0, end with 1, and be of \code{length(prop.multi) > 2}.
#' Each element represents the cutoff to separate the observations into
#' \code{length(prop.multi) - 1} mutually exclusive subgroups.
#' Default is \code{c(0, 1/3, 2/3, 1)}.
#' @param B A positive integer specifying the number of time cross-fitting is repeated in
#' \code{score.method = 'twoReg'} and \code{'contrastReg'}. Default is \code{3}.
#' @param Kfold A positive integer specifying the number of folds (parts) used in cross-fitting
#' to partition the data in \code{score.method = 'twoReg'} and \code{'contrastReg'}.
#' Default is \code{6}.
#' @param plot.gbmperf A logical value indicating whether to plot the performance measures in
#' boosting. Used only if \code{score.method = 'boosting'} or if \code{score.method = 'twoReg'}
#' or \code{'contrastReg'} and \code{initial.predictor.method = 'boosting'}. Default is \code{TRUE}.
#' @param tree.depth A positive integer specifying the depth of individual trees in boosting
#' (usually 2-3). Used only if \code{score.method = 'boosting'} or
#' if \code{score.method = 'twoReg'} or \code{'contrastReg'} and
#' \code{initial.predictor.method = 'boosting'}. Default is \code{2}.
#' @param n.trees.boosting A positive integer specifying the maximum number of trees in boosting
#' (usually 100-1000). Used only if \code{score.method = 'boosting'} or
#' if \code{score.method = 'twoReg'} or \code{'contrastReg'} and
#' \code{initial.predictor.method = 'boosting'}. Default is \code{200}.
#' @param error.maxNR A numerical value > 0 indicating the minimum value of the mean absolute
#' error in Newton Raphson algorithm. Used only if \code{score.method = 'contrastReg'}.
#' Default is \code{0.001}.
#' @param max.iterNR A positive integer indicating the maximum number of iterations in the
#' Newton Raphson algorithm. Used only if \code{score.method = 'contrastReg'}.
#' Default is \code{150}.
#' @param tune A vector of 2 numerical values > 0 specifying tuning parameters for the
#' Newton Raphson algorithm. \code{tune[1]} is the step size, \code{tune[2]} specifies a
#' quantity to be added to diagonal of the slope matrix to prevent singularity.
#' Used only if \code{score.method = 'contrastReg'}. Default is \code{c(0.5, 2)}.
#' @param train.prop A numerical value (in (0, 1)) indicating the proportion of total data used
#' for training. Default is \code{3/4}.
#' @param cv.n A positive integer value indicating the number of cross-validation iterations.
#' Default is \code{10}.
#' @param error.max A numerical value > 0 indicating the tolerance (maximum value of error)
#' for the largest standardized absolute difference in the covariate distributions or in the
#' doubly robust estimated rate ratios between the training and validation sets. This is used
#' to define a balanced training-validation splitting. Default is \code{0.1}.
#' @param max.iter A positive integer value indicating the maximum number of iterations when
#' searching for a balanced training-validation split. Default is \code{5,000}.
#' @param n.boot A numeric value indicating the number of bootstrap samples used. This is only relevant
#' if \code{inference = TRUE}. Default is \code{500}.
#' @param plot.boot A logic value indicating whether histograms of the bootstrapped log(rate ratio) should
#' be produced at every \code{n.boot/10}-th iteration and whether the final histogram should be outputted.
#' Default is \code{FALSE}.
#'
#' @return Nothing. Will stop if arguments are incorrect.

arg.checks.common <- function(fun,
                              ps.method, minPS, maxPS,
                              higher.y = NULL, abc = NULL,
                              prop.cutoff = NULL, prop.multi = NULL,
                              B = NULL, Kfold = NULL, plot.gbmperf = NULL,
                              tree.depth = NULL, n.trees.boosting = NULL,
                              error.maxNR = NULL, max.iterNR = NULL, tune = NULL,
                              train.prop = NULL, cv.n = NULL,
                              error.max = NULL, max.iter = NULL,
                              n.boot = NULL, plot.boot = NULL) {


  #### Check common arguments ####
  # Check values of ps.method
  if (!(ps.method %in% c("glm", "lasso"))) stop("ps.method must be either 'glm' or 'lasso'.")

  # Check values of minPS, maxPS
  if (minPS > 1 | minPS < 0) stop("minPS must be a number between 0 and 1, inclusive.")
  if (maxPS > 1 | maxPS < 0) stop("maxPS must be a number between 0 and 1, inclusive.")
  if (minPS >= maxPS) stop("minPS must be smaller than maxPS.")


  if (fun %in% c("cv", "pm")) {
    # Check values of higher.y
    if (!(higher.y %in% c(TRUE, FALSE))) stop("higher.y has to be boolean.")

    # Check values of prop
    if (any(prop.cutoff > 1) | any(prop.cutoff < 0)) stop("prop.cutoff must be > 0 and <= 1.")

    # Check cross-fitting arguments
    if (any(c(B, Kfold, tree.depth, n.trees.boosting) <= 0)) stop("B, Kfold, tree.depth, or n.trees.boosting must be > 0.")
    if (!(plot.gbmperf %in% c(TRUE, FALSE))) stop("plot.gbmperf has to be boolean.")

    # Check control arguments
    if (length(tune) != 2) stop("tune must be a vector with two elements.")
    if (any(c(error.maxNR, max.iterNR, tune[1], tune[2]) <= 0)) stop("error.maxNR, max.iterNR and tune must be > 0.")

    if (fun == "cv") {
      # Check values of abc
      if (!(abc %in% c(TRUE, FALSE))) stop("abc has to be boolean.")

      # Check values of prop.multi
      if (any(prop.multi > 1) | any(prop.multi < 0)) stop("prop.multi must be >= 0 and <= 1.")
      if (length(prop.multi) == 2 & sum(prop.multi == c(0, 1)) == 2) stop("prop.multi must defined at least 2 subgroups.")

      # Check values of CV
      if (train.prop >= 1| train.prop <= 0) stop("train.prop must be a number between 0 and 1, exclusive.")
      if (cv.n <= 0) stop("cv.n must be > 0.")
      if (cv.n %% 1 != 0) stop("cv.n must be an integer.")

      # Check control values for balance.split
      if (any(c(error.max, max.iter) <= 0)) stop("error.max and max.iterNR must be > 0.")
    }
  } else if (fun == "drinf") {
    # Check n.boot positive
    if (n.boot <= 0) stop("n.boot must be > 0.")

    # Check plot.boot boolean
    if (!(plot.boot %in% c(TRUE, FALSE))) stop("plot.boot has to be boolean.")
  }
}


#' Check arguments
#' Catered to all types of outcome
#' Apply at the beginning of \code{pmcount()}, \code{cvcount()}, \code{drcount.inference()}, \code{pmsurv()}, \code{cvsurv()}, and \code{drsurv.inference()}
#'
#' @param fun A function for which argument check is needed; "pm" for \code{pmcount()} and \code{pmsurv()}, "cv" for \code{cvcount()} and \code{cvsurv()},
#' and "drinf" for \code{drcount.inference()} and \code{drsurv.inference()}. No default.
#' @param response The type of response. Always 'survival' for this function.
#' @param data A data frame containing the variables in the outcome and propensity score models;
#' a data frame with \code{n} rows (1 row per observation).
#' @param followup.time Follow-up time, interpreted as the potential censoring time. If the potential censoring time is known,
#' followup.time is the name of a corresponding column in the data. Otherwise, set \code{followup.time == NULL}.
#' @param tau0 The truncation time for defining restricted mean time lost.
#' @param surv.min Lower truncation limit for probability of being censored (positive and very close to 0).
#' @param ipcw.method The censoring model. Allowed values are: \code{'breslow'} (Cox regression with Breslow estimator of the baseline survivor function),
#' \code{'aft (exponential)'}, \code{'aft (weibull)'}, \code{'aft (lognormal)'} or \code{'aft (loglogistic)'}. Default is \code{'breslow'}.
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
#' @param higher.y A logical value indicating whether higher (\code{TRUE}) or lower (\code{FALSE})
#' values of the outcome are more desirable. Default is \code{TRUE}.
#' @param score.method A vector of one or multiple methods to estimate the CATE score.
#' Allowed values are: \code{'boosting'}, \code{'poisson'}, \code{'twoReg'}, \code{'contrastReg'},
#' \code{'negBin'}. Default specifies all 5 methods.
#' @param abc A logical value indicating whether the area between curves (ABC) should be calculated
#' at each cross-validation iterations, for each \code{score.method}. Default is \code{TRUE}.
#' @param prop.cutoff A vector of numerical values (in (0, 1]) specifying percentiles of the
#' estimated log CATE scores to define nested subgroups. Each element represents the cutoff to
#' separate observations in nested subgroups (below vs above cutoff).
#' The length of \code{prop.cutoff} is the number of nested subgroups.
#' An equally-spaced sequence of proportions ending with 1 is recommended.
#' Default is \code{seq(0.5, 1, length = 6)}.
#' @param prop.multi A vector of numerical values (in [0, 1]) specifying percentiles of the
#' estimated log CATE scores to define mutually exclusive subgroups.
#' It should start with 0, end with 1, and be of \code{length(prop.multi) > 2}.
#' Each element represents the cutoff to separate the observations into
#' \code{length(prop.multi) - 1} mutually exclusive subgroups.
#' Default is \code{c(0, 1/3, 2/3, 1)}.
#' @param train.prop A numerical value (in (0, 1)) indicating the proportion of total data used
#' for training. Default is \code{3/4}.
#' @param cv.n A positive integer value indicating the number of cross-validation iterations.
#' Default is \code{10}.
#' @param error.max A numerical value > 0 indicating the tolerance (maximum value of error)
#' for the largest standardized absolute difference in the covariate distributions or in the
#' doubly robust estimated rate ratios between the training and validation sets. This is used
#' to define a balanced training-validation splitting. Default is \code{0.1}.
#' @param max.iter A positive integer value indicating the maximum number of iterations when
#' searching for a balanced training-validation split. Default is \code{5,000}.
#' @param initial.predictor.method A character vector for the method used to get initial
#' outcome predictions conditional on the covariates in \code{cate.model}
#' in \code{score.method = 'twoReg'} and \code{'contrastReg'}. Allowed values include
#' one of \code{'randomForest'}, \code{'boosting'} and \code{'logistic'} (fastest). Default is \code{'randomForest'}.
#' @param tree.depth A positive integer specifying the depth of individual trees in boosting
#' (usually 2-3). Used only if \code{score.method = 'boosting'} or
#' if \code{score.method = 'twoReg'} or \code{'contrastReg'} and
#' \code{initial.predictor.method = 'boosting'}. Default is \code{2}.
#' @param n.trees.rf A positive integer specifying the maximum number of trees in random forest.
#' Used if \code{score.method = 'ranfomForest'} or
#' if \code{initial.predictor.method = 'randomForest'} with
#' \code{score.method = 'twoReg'} or \code{'contrastReg'}. Default is \code{1000}.
#' @param n.trees.boosting A positive integer specifying the maximum number of trees in boosting
#' (usually 100-1000). Used if \code{score.method = 'boosting'} or
#' if \code{initial.predictor.method = 'boosting'} with
#' \code{score.method = 'twoReg'} or \code{'contrastReg'}. Default is \code{150}.
#' @param B A positive integer specifying the number of time cross-fitting is repeated in
#' \code{score.method = 'twoReg'} and \code{'contrastReg'}. Default is \code{3}.
#' @param Kfold A positive integer specifying the number of folds (parts) used in cross-fitting
#' to partition the data in \code{score.method = 'twoReg'} and \code{'contrastReg'}.
#' Default is \code{6}.
#' @param plot.gbmperf A logical value indicating whether to plot the performance measures in
#' boosting. Used only if \code{score.method = 'boosting'} or if \code{score.method = 'twoReg'}
#' or \code{'contrastReg'} and \code{initial.predictor.method = 'boosting'}. Default is \code{TRUE}.
#' @param error.maxNR A numerical value > 0 indicating the minimum value of the mean absolute
#' error in Newton Raphson algorithm. Used only if \code{score.method = 'contrastReg'}.
#' Default is \code{0.001}.
#' @param max.iterNR A positive integer indicating the maximum number of iterations in the
#' Newton Raphson algorithm. Used only if \code{score.method = 'contrastReg'}.
#' Default is \code{150}.
#' @param tune A vector of 2 numerical values > 0 specifying tuning parameters for the
#' Newton Raphson algorithm. \code{tune[1]} is the step size, \code{tune[2]} specifies a
#' quantity to be added to diagonal of the slope matrix to prevent singularity.
#' Used only if \code{score.method = 'contrastReg'}. Default is \code{c(0.5, 2)}.
#' @param n.boot A numeric value indicating the number of bootstrap samples used. This is only relevant
#' if \code{inference = TRUE}. Default is \code{500}.
#' @param plot.boot A logic value indicating whether histograms of the bootstrapped log(rate ratio) should
#' be produced at every \code{n.boot/10}-th iteration and whether the final histogram should be outputted.
#' Default is \code{FALSE}.
#' @param interactions A logical value indicating whether the outcome model should assume interactions
#' between \code{x} and \code{trt}. If \code{TRUE}, interactions will be assumed only if at least 10 patients
#' received each treatment option. Default is \code{TRUE}.
#'
#' @return Nothing. Will stop if arguments are incorrect.

arg.checks <- function(fun, response, data,
                       followup.time = NULL, tau0 = NULL, surv.min = NULL,
                       ipcw.method = NULL,
                       ps.method, minPS, maxPS,
                       higher.y = NULL, score.method = NULL, abc = NULL,
                       prop.cutoff = NULL, prop.multi = NULL,
                       train.prop = NULL, cv.n = NULL,
                       error.max = NULL, max.iter = NULL,
                       initial.predictor.method = NULL,
                       tree.depth = NULL, n.trees.rf = NULL, n.trees.boosting = NULL,
                       B = NULL, Kfold = NULL, plot.gbmperf = NULL,
                       error.maxNR = NULL, max.iterNR = NULL, tune = NULL,
                       n.boot = NULL, plot.boot = NULL, interactions = NULL){

  #### Check common arguments ####
  arg.checks.common(fun = fun,
                    ps.method = ps.method, minPS = minPS, maxPS = maxPS,
                    higher.y = higher.y,
                    prop.cutoff = prop.cutoff, prop.multi = prop.multi,
                    B = B, Kfold = Kfold, plot.gbmperf = plot.gbmperf,
                    tree.depth = tree.depth, n.trees.boosting = n.trees.boosting,
                    error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune,
                    abc = abc,
                    train.prop = train.prop, cv.n = cv.n,
                    error.max = error.max, max.iter = max.iter,
                    n.boot = n.boot, plot.boot = plot.boot)

  #### Check argument only used for count outcome ####
  if (response == "count") {
    # Check initial predictor method
    if (is.null(initial.predictor.method) == FALSE) {
      if (!(initial.predictor.method %in% c("poisson", "boosting", "gam"))) stop("initial.predictor.method must be 'poisson', 'boosting' or 'gam'.")
    }

    # Check values of score.method
    if (any(!(score.method %in% c("boosting", "poisson", "twoReg", "contrastReg", "negBin")))) stop("Elements of score.method must come from: 'boosting', 'poisson', 'twoReg', 'contrastReg', 'negBin'.")

    if (fun == "drinf"){
      # Check interactions boolean
      if (!(interactions %in% c(TRUE, FALSE))) stop("interactions has to be boolean.")
    }
  }

  #### Check argument only used for survival outcome ####
  if (response == "survival") {
    # Check if followup.time is either NULL or a column name of the data
    if (!is.null(followup.time)) {
      if(!as.character(followup.time) %in% names(data)) stop("followup.time must be either NULL or one of the column names of data.")
    }

    # Check if tau0 > 0
    if (is.null(tau0) == FALSE) {
      if (tau0 <= 0) stop("tau0 must be positive.")
    }

    # Check if surv.min is positive and very close to 0
    if (!(surv.min > 0 & surv.min < 0.1)) stop("surv.min must be positive and close to 0.")

    # Check n.trees.rf
    if (any(c(n.trees.rf) <= 0)) stop("n.trees.rf must be > 0.")

    # Check values for ipcw.method
    if (!(ipcw.method %in% c("breslow", "aft (exponential)", "aft (weibull)", "aft (lognormal)", "aft (loglogistic)"))) stop("ipcw.method must be either 'breslow', 'aft (exponential)', 'aft (weibull)', 'aft (lognormal)', 'aft (loglogistic)'.")

    # Check initial predictor method
    if (is.null(initial.predictor.method) == FALSE) {
      if (!(initial.predictor.method %in% c("randomForest", "boosting", "logistic"))) stop("initial.predictor.method must be 'randomForest', 'boosting' or 'logistic'.")
    }

    # Check values of score.method
    if (any(!(score.method %in% c("boosting", "poisson", "twoReg", "contrastReg", "randomForest")))) stop("Elements of score.method must come from: 'boosting', 'poisson', 'twoReg', 'contrastReg', 'randomForest'.")

  }
}


#' Data preprocessing
#' Apply at the beginning of \code{pmcount()}, \code{cvcount()}, \code{pmsurv()}, and \code{cvsurv()}, after \code{arg.checks()}
#'
#' @param fun A function for which argument check is needed; "pm" for \code{pmcount()} and \code{pmsurv()}, "cv" for \code{cvcount()} and \code{cvsurv()},
#' and "drinf" for \code{drcount.inference()} and \code{drsurv.inference()}. No default.
#' @param cate.model A formula describing the outcome model to be fitted.
#' The outcome must appear on the left-hand side.
#' @param ps.model A formula describing the propensity score model to be fitted.
#' The treatment must appear on the left-hand side. The treatment must be a numeric vector
#' coded as 0/1. If data are from a RCT, specify \code{ps.model} as an intercept-only model.
#' @param ipcw.model A formula describing inverse probability of censoring weighting(IPCW) model to be fitted.
#' If covariates are the same as outcome model, set \code{ipcw.model = NULL}.
#' Otherwise, the left-hand side must be empty and the right-hand side is a covariates model.
#' @param tau0 The truncation time for defining restricted mean time lost. Default is \code{NULL},
#' which corresponds to setting the truncation time as the maximum survival time in the data
#' @param data A data frame containing the variables in the outcome, propensity score, and IPCW models;
#' a data frame with \code{n} rows (1 row per observation).
#' @param prop.cutoff A vector of numerical values (in (0, 1]) specifying percentiles of the
#' estimated log CATE scores to define nested subgroups. Each element represents the cutoff to
#' separate observations in nested subgroups (below vs above cutoff).
#' The length of \code{prop.cutoff} is the number of nested subgroups.
#' An equally-spaced sequence of proportions ending with 1 is recommended.
#' Default is \code{seq(0.5, 1, length = 6)}.
#' @param prop.multi A vector of numerical values (in [0, 1]) specifying percentiles of the
#' estimated log CATE scores to define mutually exclusive subgroups.
#' It should start with 0, end with 1, and be of \code{length(prop.multi) > 2}.
#' Each element represents the cutoff to separate the observations into
#' \code{length(prop.multi) - 1} mutually exclusive subgroups.
#' Default is \code{c(0, 1/3, 2/3, 1)}.
#' @param ps.method A character value for the method to estimate the propensity score.
#' Allowed values include one of:
#' \code{'glm'} for logistic regression with main effects only (default), or
#' \code{'lasso'} for a logistic regression with main effects and LASSO penalization on
#' two-way interactions (added to the model if interactions are not specified in \code{ps.model}).
#' Relevant only when \code{ps.model} has more than one variable.
#' @param initial.predictor.method A character vector for the method used to get initial
#' outcome predictions conditional on the covariates. Only applies when \code{score.method}
#' includes \code{'twoReg'} or \code{'contrastReg'}. Allowed values include one of
#' \code{'randomForest'} (survival outcomes only), \code{'boosting'}, \code{'logistic'}
#' (survival outcomes only, fast), \code{'poisson'} (count outcomes only, fast), and
#' \code{'gam'} (count outcomes only). Default is \code{NULL}, which assigns \code{'boosting'}
#' for count outcomes and \code{'randomForest'} for survival outcomes.
#' @param response The type of response variables; \code{count} (default) or \code{survival}.
#'
#' @return A list of elements:
#'            - y: outcome; vector of length \code{n} (observations)
#'            - d : the event indicator; vector of length \code{n}; only if \code{respone = "survival"}
#'            - trt: binary treatment; vector of length \code{n}
#'            - x.ps: matrix of \code{p.ps} baseline covariates specified in the propensity score model (plus intercept); dimension \code{n} by \code{p.ps + 1}
#'            - x.cate: matrix of \code{p.cate} baseline covariates specified in the outcome model; dimension \code{n} by \code{p.cate}
#'            - x.ipcw: matrix of \code{p.ipw} baseline covarites specified in inverse probability of censoring weighting model; dimension \code{n} by \code{p.ipw}
#'            - time: offset; vector of length \code{n}; only if \code{response = "count"}
#'            - if \code{fun = "pm"}:
#'                - prop: formatted \code{prop.cutoff}
#'                - prop.no1: formatted \code{prop.cutoff} with 1 removed if applicable; otherwise prop.no1 is the same as prop
#'            - if \code{fun = "cv"}
#'                - prop.onlyhigh: formatted \code{prop.cutoff} with 0 removed if applicable
#'                - prop.bi; formatted \code{prop.cutoff} with 0 and 1 removed if applicable
#'                - prop.multi: formatted \code{prop.multi}, starting with 0 and ending with 1
#'


data.preproc.surv <- function(fun, cate.model, ps.model, ipcw.model = NULL, tau0 = NULL, data, prop.cutoff = NULL, prop.multi = NULL,
                              ps.method, initial.predictor.method = NULL, response = "count") {

  ## cate.model and ps.model as formulas
  cate.model <- as.formula(cate.model)
  ps.model <- as.formula(ps.model)

  ## Extraction step
  ### If count data: extract y, time, trt, x.cate, x.ps in matrix form from cate.model and ps.model
  ### If survival data: extract y, d, trt, x.cate, x.ps, x.ipcw in matrix form from cate.model, ps.model and ipcw.model
  cate.mat <- model.frame(cate.model, data, na.action = 'na.pass')
  ps.mat <- model.frame(ps.model, data, na.action = 'na.pass')
  resp <- model.response(cate.mat)
  if (is.null(resp)) {
    stop("Outcome must be supplied on the left-hand side of the cate.model formula.")
  } else {
    if(response == "count") {
      y <- resp
      time <- model.offset(cate.mat)
      if (is.null(time) == TRUE) { # Eventually some outcomes will not need an offset
        time <- rep(0, nrow(data))
        warning("No offset supplied. Offset set to 0.")
      }
    } else if (response == "survival") {
      y <- resp[, 1]
      d <- resp[, 2]
    }
  }

  # Assign default to tau0 if NULL
  if (response == "survival" & is.null(tau0) == TRUE) {
    tau0 <- max(y[d == 1])
    warning("No value supplied for tau0. Default sets tau0 to the maximum survival time.")
  }

  # trt
  trt <- model.response(ps.mat)
  if (is.null(trt) == TRUE) stop("Treatment must be supplied on the left-hand side of the ps.model formula.")

  # Covariate matrices
  x.cate <- model.matrix(cate.model, cate.mat)[, -1, drop = FALSE]
  x.ps <- model.matrix(ps.model, ps.mat)
  if (ncol(x.ps) == 1 & ps.method == "lasso") stop("LASSO penalization irrelevant when ps.model specified as a function of an intercept only. Consider setting ps.method='glm'.")

  # Covariate matrices for IPCW
  if (response == "survival") {
    if (is.null(ipcw.model) == TRUE) {
      x.ipcw <- cbind(x.cate, trt)
    } else {
      ipcw.model <- as.formula(ipcw.model)
      ipcw.mat <- model.frame(ipcw.model, data, na.action = 'na.pass')
      x.ipcw <- model.matrix(ipcw.model, ipcw.mat)[, -1, drop = FALSE]
    }
  }

  ## Check negative count or survival
  if (any(y < 0)) stop(paste0("Negative y values not allowed for ", response, " outcomes."))

  ## Check missing data
  if (response == "count") {
    if (any(is.na(y)) | any(is.na(time)) | any(is.na(trt)) | any(is.na(x.cate)) | any(is.na(x.ps))) stop("Missing data not allowed in cate.model, ps.model or offset.")
  } else if (response == "survival") {
    if (any(is.na(y)) | any(is.na(trt)) | any(is.na(x.cate)) | any(is.na(x.ps)) | any(is.na(x.ipcw))) stop("Missing data not allowed in cate.model or ps.model")
  }

  ## Check treatment binary numeric and coded as 0/1
  cat <- warn.dr <- NULL
  if (as.numeric(length(unique(trt))) > 2) {
    stop("trt must be binary.")
  } else if (!all(unique(trt) %in% c(0, 1))) {
    cat <- cat <- sort(unique(trt))
    trt <- ifelse(trt == cat[1], 0, 1)
    warn.dr <- warning(paste0("Variable trt was recoded to 0/1 with ", cat[1], "->0 and ", cat[2], "->1.\n"))
  } else if (is.factor(trt) == TRUE) {
    trt <- as.numeric(trt == 1)
  }

  ## Check event binary numeric and coded as 0/1 for survival outcomes
  if (response == "survival") {
    if (sum(is.na(d)) == length(d)) {
      stop("The event indicator must be supplied as 0/1.")
    } else if (any(is.na(d))) {
      stop("Missing data in event indicator not allowed.")
    } else if (as.numeric(length(unique(d))) > 2) {
      stop("The event indicator must be binary.")
    }
  }

  ## Assign default to initial.predictor.method if NULL
  if (is.null(initial.predictor.method) == TRUE & fun %in% c("pm", "cv")) {
    if (response == "count") {
      initial.predictor.method <- "boosting"
    } else if (response == "survival") {
      initial.predictor.method <- "randomForest"
    }
  }

  ## Output data
  if (response == "count") {
    out.data <- list(y = y, time = time, trt = trt, x.cate = x.cate, x.ps = x.ps, cat.trt = cat, tau0 = tau0, initial.predictor.method = initial.predictor.method)
  } else if (response == "survival") {
    out.data <- list(y = y, d = d, trt = trt, x.cate = x.cate, x.ps = x.ps, x.ipcw = x.ipcw, cat.trt = cat, tau0 = tau0, initial.predictor.method = initial.predictor.method)
  }

  if (fun == "pm") {
    ## Check values of prop
    prop <- sort(prop.cutoff) # sort the proportions from small to large
    if (prop[1] == 0) {
      prop <- prop[-1] # if first element is 0, remove it because this means we leave out 0% of individuals
      warning("The first element of prop.cutoff cannot be 0 and has been removed.")
    }
    out.data$prop.no1 <- prop
    if (prop[length(prop)] == 1){
      out.data$prop.no1 <- prop[-length(prop)]
    }
    out.data$prop <- prop
    return(out.data)

  } else if (fun == "cv") {
    ## Check values of prop.cutoff and prop.multi
    prop.onlyhigh <- prop.bi <- sort(prop.cutoff)
    prop.multi <- sort(prop.multi)
    if (prop.onlyhigh[1] == 0) {
      prop.onlyhigh <- prop.onlyhigh[-1]
      warning("The first element of prop.cutoff cannot be 0 and has been removed.") # only need to show once so no warning in the prop.bi check below
    }
    if (prop.bi[1] == 0) {
      prop.bi <- prop.bi[-1]
    }
    if (prop.bi[length(prop.bi)] == 1) {
      prop.bi <- prop.bi[-length(prop.bi)]
    }
    if (prop.multi[1] != 0) {
      prop.multi <- c(0, prop.multi)
      warning("The first element of prop.multi must be 0 and has been added.")
    }
    if (prop.multi[length(prop.multi)] != 1) {
      prop.multi <- c(prop.multi, 1)
      warning("The last element of prop.multi must be 1 and has been added.")
    }

    out.data$prop.onlyhigh <- prop.onlyhigh
    out.data$prop.bi <- prop.bi
    out.data$prop.multi <- prop.multi

    return(out.data)

  } else if (fun == "drinf") {
    out.data$warning <- warn.dr
    return(out.data)
  }
}



#' Catch errors and warnings when estimating the ATEs in the nested subgroup
#'
#' Storing the errors and warnings that occurred when estimating the ATEs in the nested subgroups.
#' If there are no errors and no warnings, the estimated log.rmtl.ratio and log.hazard.ratio are provided.
#' If there are warnings but no errors, the estimated log.rmtl.ratio and log.hazard.ratio are provided with a warning attribute set.
#' If there are errors, the NA values are returned for log.rmtl.ratio and log.hazard.ratio. A error attribute set is also provided.
#'
#' @param fun The drsurv function...
#' @return A list containing

survCatch <- function(fun) {
  err.msg <- err.call <- NULL
  warn.msg <- warn.call <- NULL
  fit <- withCallingHandlers(
    tryCatch(expr = fun,
             error = function(e) {
               err.msg <<- conditionMessage(e)
               err.call <<- conditionCall(e)
               list(log.rmtl.ratio = NA, log.hazard.ratio = NA)
             }),
    warning = function(w) {
      warn.msg <<- append(warn.msg, conditionMessage(w))
      warn.call <<- append(warn.call, conditionCall(w))
      invokeRestart("muffleWarning")
    })

  if (!is.null(err.msg)) {
    err.call <- toString(as.expression(err.call))
    errors <- paste0("Error in ", err.call, " : ", err.msg)
  } else {
    errors <- NULL
  }

  if (!is.null(warn.msg)) {
    warn.call <- toString(as.expression(warn.call))
    warnings <- paste0("Warning in ", warn.call, " : ", warn.msg)
  } else {
    warnings <- NULL
  }

  out <- c()
  out$log.rmtl.ratio <- fit$log.rmtl.ratio
  out$log.hazard.ratio <- fit$log.hazard.ratio
  out$warnings <- warnings
  out$errors <- errors
  return(out)
}
