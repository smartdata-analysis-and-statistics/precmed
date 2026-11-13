#' Propensity score estimation with a linear model
#'
#' Propensity score based on a multivariate logistic regression with main effects only
#'
#' @param trt Treatment received; vector of size \code{n} (observations) with treatment coded as 0/1
#' @param x.ps A matrix of \code{p.ps} baseline covariates (plus a leading column of 1 for the intercept);
#' dimension \code{n} by \code{p.ps + 1} (covariates in the propensity score model plus intercept)
#' @param xnew A matrix of \code{p.ps} baseline covariates (plus a leading column of 1 for the intercept)
#' for which we want PS predictions; dimension \code{m} (observations in the new data set) by \code{p.ps + 1}
#' @param minPS A numerical value (in `[0, 1]`) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in `(0, 1]`) above which estimated propensity scores should be
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
#' @param minPS A numerical value (in `[0, 1]`) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in `(0, 1]`) above which estimated propensity scores should be
#' truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.
#'
#' @param solver Interaction-screening backend: \code{"C"} (default) or
#' \code{"R"}. C screens pairs before allocating retained training and prediction
#' columns. Both use the strict SMD > 0.1 rule and positive sample variance.
#'
#' @return The trimmed propensity score for each unit; vector of size \code{n} (if \code{xnew} is NULL) or \code{m}

glm.ps <- function(trt, x.ps, xnew = NULL, minPS = 0.01, maxPS = 0.99,
                   solver = c("C", "R")) {
  solver <- match.arg(solver)
  x <- x.ps[, -1, drop = FALSE]
  p <- length(x[1, ])
  if (!is.null(xnew)) {
    xnew <- as.matrix(xnew)[, -1, drop = FALSE]
  }
  interactions <- ps_interactions(x, trt, xnew, solver)
  xmat <- interactions$train
  if (!is.null(xnew)) xmatnew <- interactions$prediction

  xmat <- cbind(x, xmat)
  q <- dim(xmat)[2]
  if (is.null(xnew) == FALSE) xmatnew <- cbind(xnew, xmatnew)

  if (q > p) {
    # Estimate PS with LASSO with penalty term only for interactions
    pf0 <- c(rep(0, p), rep(1, q - p))
    fit.cv <- cv.glmnet(x = xmat, y = trt, family = "binomial", penalty.factor = pf0)
    # cv.glmnet already fitted this same full-data regularization path.
    beta <- coef(fit.cv$glmnet.fit, s = fit.cv$lambda.min)
    if (is.null(xnew) == TRUE) {
      prob <- as.vector(1 / (1 + exp(-cbind(1, xmat) %*% beta)))
    } else {
      prob <- as.vector(1 / (1 + exp(-cbind(1, xmatnew) %*% beta)))
    }
  } else if (q == p) {
    # Estimate PS with logistic regression
    fit <- if (p == 0) glm(trt ~ 1, family = "binomial") else
      glm(trt ~ x, family = "binomial")
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







