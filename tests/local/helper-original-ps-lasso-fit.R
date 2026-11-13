# Frozen PS implementation immediately before reusing the CV full-data fit.
original_ps_lasso_fit <- function(trt, x.ps, xnew = NULL, minPS = 0.01, maxPS = 0.99,
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
