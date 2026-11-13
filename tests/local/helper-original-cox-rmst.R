# Frozen pre-native RMST implementation.
original_cox_rmst <- function(y, d, x.cate, xnew, tau0) {
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
