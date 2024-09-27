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







