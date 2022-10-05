# ------------------------------------------------------------------
#
# Project: Precision Medicine MS (precmed) - Comprehensive R package
#
# Purpose: Utility functions for Count outcomes
#
# Platform: Windows
# R Version: 4.1.0
#


#' Propensity score estimation with a linear model
#'
#' Propensity score based on a multivariate logistic regression with main effects only
#'
#' @param trt Treatment received; vector of size \code{n} (observations) with treatment coded as 0/1
#' @param x.ps A matrix of \code{p.ps} baseline covariates (plus a leading column of 1 for the intercept);
#' dimension \code{n} by \code{p.ps + 1} (covariates in the propensity score model plus intercept)
#' @param xnew A matrix of \code{p.ps} baseline covariates (plus a leading column of 1 for the intercept)
#' for which we want PS predictions; dimension \code{m} (observations in the new data set) by \code{p.ps + 1}
#' @param minPS A numerical value (in [0, 1]) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in (0, 1]) above which estimated propensity scores should be
#' truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.
#'
#' @return The estimated propensity score for each unit; vector of size \code{n} (if \code{xnew} is NULL) or \code{m}

glm.simplereg.ps <- function(trt, x.ps, xnew = NULL, minPS = 0.01, maxPS = 0.99) {
  fit <- glm(trt ~ -1 + x.ps, family = "binomial")
  beta <- fit$coef
  if (is.null(xnew) == TRUE) {
    x <- x.ps[, is.na(beta) == FALSE, drop = FALSE]
    prob <- as.vector(1 / (1 + exp(-x %*% beta[is.na(beta) == FALSE])))
  } else {
    x <- xnew[, is.na(beta) == FALSE, drop = FALSE]
    prob <- as.vector(1 / (1 + exp(-x %*% beta[is.na(beta) == FALSE])))
  }
  prob <- pmin(prob, maxPS)
  prob <- pmax(prob, minPS)
  return(prob)
}


#' Propensity score estimation with LASSO
#'
#' Propensity score based on a multivariate logistic regression with LASSO penalization on the two-way interactions
#'
#' @param trt Treatment received; vector of size \code{n} (observations) with treatment coded as 0/1
#' @param x.ps Matrix of \code{p.ps} baseline covariates (plus a leading column of 1 for the intercept);
#' dimension \code{n} by \code{p.ps + 1} (covariates in the propensity score model plus intercept)
#' @param xnew Matrix of \code{p.ps} baseline covariates (plus a leading column of 1 for the intercept)
#' for which we want propensity scores predictions; dimension \code{m} (observations in the new data set) by \code{p.ps + 1}
#' @param minPS A numerical value (in [0, 1]) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in (0, 1]) above which estimated propensity scores should be
#' truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.
#'
#' @return The trimmed propensity score for each unit; vector of size \code{n} (if \code{xnew} is NULL) or \code{m}

glm.ps <- function(trt, x.ps, xnew = NULL, minPS = 0.01, maxPS = 0.99) {
  x <- x.ps[, -1, drop = FALSE]
  p <- length(x[1, ])
  xmat <- matrix(NA, nrow = nrow(x), ncol = p * (p - 1) / 2)
  if (is.null(xnew) == FALSE) {
    xnew <- as.matrix(xnew)
    xnew <- xnew[, -1, drop = FALSE]
    xmatnew <- matrix(NA, nrow = nrow(xnew), ncol = p * (p - 1) / 2)
  }

  # Identify two-way interactions with SMD > 0.1 between treatment groups
  if (p > 1) {
    for (i in 1:(p - 1)) {
      a <- 1 + ((p + 1 - 2) + (p + 1 - i)) * (i - 1) / 2
      xmat[, a:(a + p - 1 - i)] <- x[, (i + 1):p] * x[, i]
      if (is.null(xnew) == FALSE) xmatnew[, a:(a + p - 1 - i)] <- xnew[, (i + 1):p] * xnew[, i]
    }
    mu1 <- apply(xmat[trt == 1, ], 2, mean)
    mu0 <- apply(xmat[trt == 0, ], 2, mean)
    sigma1 <- apply(xmat[trt == 1, ], 2, var)
    sigma0 <- apply(xmat[trt == 0, ], 2, var)
    diff <- sqrt(2) * abs(mu1 - mu0) / sqrt(sigma1 + sigma0 + 1e-6)
    xmat <- xmat[, (diff > 0.1 & sigma1 + sigma0 > 0)]
    if (is.null(xnew) == FALSE) xmatnew <- xmatnew[, (diff > 0.1 & sigma1 + sigma0 > 0)]
  }

  xmat <- cbind(x, xmat)
  q <- dim(xmat)[2]
  if (is.null(xnew) == FALSE) xmatnew <- cbind(xnew, xmatnew)

  if (q > p) {
    # Estimate PS with LASSO with penalty term only for interactions
    pf0 <- c(rep(0, p), rep(1, q - p))
    fit <- glmnet(x = xmat, y = trt, family = "binomial", penalty.factor = pf0)
    fit.cv <- cv.glmnet(x = xmat, y = trt, family = "binomial", penalty.factor = pf0)
    beta <- coef(fit, s = fit.cv$lambda.min)
    if (is.null(xnew) == TRUE) {
      prob <- as.vector(1 / (1 + exp(-cbind(1, xmat) %*% beta)))
    } else {
      prob <- as.vector(1 / (1 + exp(-cbind(1, xmatnew) %*% beta)))
    }
  } else if (q == p) {
    # Estimate PS with logistic regression
    fit <- glm(trt ~ x, family = "binomial")
    beta <- fit$coef
    if (is.null(xnew) == TRUE) {
      prob <- as.vector(1 / (1 + exp(-cbind(1, x)[, is.na(beta) == FALSE, drop = FALSE] %*% beta[is.na(beta) == FALSE])))
    } else {
      prob <- as.vector(1 / (1 + exp(-cbind(1, xnew)[, is.na(beta) == FALSE, drop = FALSE] %*% beta[is.na(beta) == FALSE])))
    }
  }

  prob <- pmin(prob, maxPS)
  prob <- pmax(prob, minPS)
  return(prob)
}


#' Split the given dataset into balanced training and validation sets
#' (within a pre-specified tolerance)
#' Balanced means 1) The ratio of treated and controls is maintained in the training and validation sets
#'                2) The covariate distributions are balanced between the training and validation sets
#'
#' @param y Observed outcome; vector of size \code{n} (observations)
#' @param trt Treatment received; vector of size \code{n} with treatment coded as 0/1
#' @param x.cate Matrix of \code{p.cate} baseline covariates; dimension \code{n} by \code{p.cate}
#' (covariates in the outcome model)
#' @param x.ps Matrix of \code{p.ps} baseline covariates (plus a leading column of 1 for the intercept);
#' dimension \code{n} by \code{p.ps + 1} (covariates in the propensity score model plus intercept)
#' @param time Log-transformed person-years of follow-up; vector of size \code{n}
#' @param minPS A numerical value (in [0, 1]) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in (0, 1]) above which estimated propensity scores should be
#' truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.
#' @param train.prop A numerical value (in (0, 1)) indicating the proportion of total data used
#' for training. Default is \code{3/4}.
#' @param error.max A numerical value > 0 indicating the tolerance (maximum value of error)
#' for the largest standardized absolute difference in the covariate distributions or in the
#' doubly robust estimated rate ratios between the training and validation sets. This is used
#' to define a balanced training-validation splitting. Default is \code{0.1}.
#' @param max.iter A positive integer value indicating the maximum number of iterations when
#' searching for a balanced training-validation split. Default is \code{5,000}.
#'
#' @return A list of 10 objects, 5 training and 5 validation of y, trt, x.cate, x.ps, time:
#'             y.train          - observed outcome in the training set; vector of size \code{m} (observations in the training set)
#'             trt.train        - treatment received in the training set; vector of size \code{m} coded as 0/1
#'             x.cate.train     - baseline covariates for the outcome model in the training set; matrix of dimension \code{m} by \code{p.cate}
#'             x.ps.train       - baseline covariates (plus intercept) for the propensity score model in the training set; matrix of dimension \code{m} by \code{p.ps + 1}
#'             time.train       - log-transformed person-years of follow-up in the training set; vector of size \code{m}
#'             y.valid          - observed outcome in the validation set; vector of size \code{n-m}
#'             trt.valid        - treatment received in the validation set; vector of size \code{n-m} coded as 0/1
#'             x.cate.valid     - baseline covariates for the outcome model in the validation set; matrix of dimension \code{n-m} by \code{p.cate}
#'             x.ps.valid       - baseline covariates (plus intercept) for the propensity score model in the validation set; matrix of dimension \code{n-m} by \code{p.ps + 1}
#'             time.valid       - log-transformed person-years of follow-up in the validation set; vector of size \code{n-m}
balance.split <- function(y, trt, x.cate, x.ps, time,
                          minPS = 0.01, maxPS = 0.99,
                          train.prop = 3/4, error.max = 0.1, max.iter = 5000) {

  x <- cbind(x.cate, x.ps[, -1])
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
  error <-  errorbest <- Inf
  bestid.valid <- rep(NA, n - m)
  iter <- 0

  while ((error > error.max | is.na(error) == TRUE) & iter <= max.iter) {
    id.valid <- c(sample(x = id1, size = n1 - m1, replace = FALSE), sample(x = id0, size = n0 - m0, replace = FALSE))

    y.valid <- y[id.valid]
    trt.valid <- trt[id.valid]
    x.valid <- x[id.valid, , drop = FALSE]
    x.cate.valid <- x.cate[id.valid, , drop = FALSE]
    x.ps.valid <- x.ps[id.valid, , drop = FALSE]
    time.valid <- time[id.valid]

    y.train <- y[-id.valid]
    trt.train <- trt[-id.valid]
    x.train <- x[-id.valid, , drop = FALSE]
    x.cate.train <- x.cate[-id.valid, , drop = FALSE]
    x.ps.train <- x.ps[-id.valid, , drop = FALSE]
    time.train <- time[-id.valid]

    diffx <- colMeans(x.train) - colMeans(x.valid)
    errorx <- max(abs(diffx / sdx))

    est.valid  <- drcount(y = y.valid, x.cate = x.cate.valid, x.ps = x.ps.valid, trt = trt.valid, time = time.valid,
                          ps.method = "glm", minPS = minPS, maxPS = maxPS, interactions = TRUE)$log.rate.ratio
    est.train <- drcount(y = y.train, x.cate = x.train, x.ps = x.ps.train, trt = trt.train, time = time.train,
                         ps.method = "glm", minPS = minPS, maxPS = maxPS, interactions = TRUE)$log.rate.ratio

    error <- max(2 * abs(est.valid - est.train) / (abs(est.valid) + abs(est.train)), errorx)
    if (is.na(error) == FALSE & error < errorbest) {
      bestid.valid <- id.valid
      errorbest <- error
    }
    iter <- iter + 1
  }

  if (all(is.na(bestid.valid)) == TRUE) stop("No balanced data split found.")

  if (iter == max.iter + 1) {
    y.valid <- y[bestid.valid]
    trt.valid <- trt[bestid.valid]
    x.cate.valid <- x.cate[bestid.valid, , drop = FALSE]
    x.ps.valid <- x.ps[bestid.valid, , drop = FALSE]
    time.valid <- time[bestid.valid]

    y.train <- y[-bestid.valid]
    trt.train <- trt[-bestid.valid]
    x.cate.train <- x.cate[-bestid.valid, , drop = FALSE]
    x.ps.train <- x.ps[-bestid.valid, , drop = FALSE]
    time.train <- time[-bestid.valid]
    warning(paste("Maximum iteration reached and the SMD between training and validation set is still greater than error.max (error=", round(errorbest, 4), "). Consider increasing max.iter, decreasing error.max, or increasing sample size.", sep = ""))
  }

  return(list(y.train = y.train, trt.train = trt.train, x.cate.train = x.cate.train, x.ps.train = x.ps.train, time.train = time.train,
              y.valid = y.valid, trt.valid = trt.valid, x.cate.valid = x.cate.valid, x.ps.valid = x.ps.valid, time.valid = time.valid))
}



#' Data preprocessing
#' Apply at the beginning of \code{pmcount()} and \code{cvcount()}, after \code{arg.checks()}
#'
#' @param fun A function for which argument check is needed; "pm" for \code{pmcount()}, "cv" for \code{cvcount()},
#' and "drinf" for \code{drcount.inference()}. No default.
#' @param cate.model A formula describing the outcome model to be fitted.
#' The outcome must appear on the left-hand side.
#' @param ps.model A formula describing the propensity score model to be fitted.
#' The treatment must appear on the left-hand side. The treatment must be a numeric vector
#' coded as 0/1. If data are from a RCT, specify \code{ps.model} as an intercept-only model.
#' @param data A data frame containing the variables in the outcome and propensity score models;
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
#' \code{'boosting'}, \code{'poisson'} (fast), and \code{'gam'}. Default is \code{NULL}, which assigns
#' \code{'boosting'} for count outcomes.
#'
#' @return A list of 6 elements:
#'            - y: outcome; vector of length \code{n} (observations)
#'            - trt: binary treatment; vector of length \code{n}
#'            - x.ps: matrix of \code{p.ps} baseline covariates (plus intercept); dimension \code{n} by \code{p.ps + 1}
#'            - x.cate: matrix of \code{p.cate} baseline covariates; dimension \code{n} by \code{p.cate}
#'            - time: offset; vector of length \code{n}
#'            - if \code{fun = "pm"}:
#'                - prop: formatted \code{prop.cutoff}
#'            - if \code{fun = "cv"}
#'                - prop.onlyhigh: formatted \code{prop.cutoff} with 0 removed if applicable
#'                - prop.bi; formatted \code{prop.cutoff} with 0 and 1 removed if applicable
#'                - prop.multi: formatted \code{prop.multi}, starting with 0 and ending with 1

data.preproc <- function(fun, cate.model, ps.model, data, prop.cutoff = NULL, prop.multi = NULL, ps.method, initial.predictor.method = NULL) {

  ## cate.model and ps.model as formulas
  cate.model <- as.formula(cate.model)
  ps.model <- as.formula(ps.model)

  ## Extraction step: extract y, trt, x.cate, x.ps, time in matrix form from cate.model and ps.model
  cate.mat <- model.frame(cate.model, data, na.action = 'na.pass')
  ps.mat <- model.frame(ps.model, data, na.action = 'na.pass')

  # y
  y <- model.response(cate.mat)
  if (is.null(y) == TRUE) stop("Outcome must be supplied on the left-hand side of the cate.model formula.")

  # trt
  trt <- model.response(ps.mat)
  if (is.null(trt) == TRUE) stop("Treatment must be supplied on the left-hand side of the ps.model formula.")

  # Covariate matrices
  x.cate <- model.matrix(cate.model, cate.mat)[, -1, drop = FALSE]
  x.ps <- model.matrix(ps.model, ps.mat)
  if (ncol(x.ps) == 1 & ps.method == "lasso") stop("LASSO penalization irrelevant when ps.model specified as a function of an intercept only. Consider setting ps.method='glm'.")

  # time
  time <- model.offset(cate.mat)
  if (is.null(time) == TRUE) { # Eventually some outcomes will not need an offset
    time <- rep(0, nrow(data))
    warning("No offset supplied. Offset set to 0.")
  }

  ## Check missing data
  if (any(is.na(y)) | any(is.na(trt)) | any(is.na(x.cate)) | any(is.na(x.ps)) | any(is.na(time))) stop("Missing data not allowed in cate.model, ps.model or offset.")

  ## Check negative count
  if (any(y < 0)) stop("Negative y values not allowed for count outcomes.")

  ## Check treatment binary numeric and coded as 0/1
  cat <- warn.dr <- NULL
  if (as.numeric(length(unique(trt))) != 2) {
    stop("trt should describe two distinct treatments.")
  } else if (!all(unique(trt) %in% c(0, 1))) {
    cat <- sort(unique(trt))
    trt <- ifelse(trt == cat[1], 0, 1)
    warn.dr <- warning(paste0("Variable trt was recoded to 0/1 with ", cat[1], "->0 and ", cat[2], "->1.\n"))
  } else if (is.factor(trt) == TRUE) {
    trt <- as.numeric(trt == 1)
  }

  ## Assign default to initial.predictor.method if NULL
  if (is.null(initial.predictor.method) == TRUE & fun %in% c("pm", "cv")) initial.predictor.method <- "boosting"

  if (fun == "pm") {
    ## Check values of prop
    prop <- sort(prop.cutoff) # sort the proportions from small to large
    if (prop[1] == 0) {
      prop <- prop[-1] # if first element is 0, remove it because this means we leave out 0% of individuals
      warning("The first element of prop.cutoff cannot be 0 and has been removed.")
    }

    return(list(y = y, trt = trt, x.cate = x.cate, x.ps = x.ps, time = time, prop = prop))
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

    return(list(y = y, trt = trt, x.cate = x.cate, x.ps = x.ps, time = time,
                prop.onlyhigh = prop.onlyhigh, prop.bi = prop.bi, prop.multi = prop.multi, cat.trt = cat, initial.predictor.method = initial.predictor.method))
  } else if (fun == "drinf") {
    return(list(y = y, trt = trt, x.cate = x.cate, x.ps = x.ps, time = time, warning = warn.dr))
  }
}
