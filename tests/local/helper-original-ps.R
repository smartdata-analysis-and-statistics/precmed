# Frozen PS implementation before native interaction screening.
original_glm_ps <- function(trt, x.ps, xnew = NULL, minPS = 0.01, maxPS = 0.99) {
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
