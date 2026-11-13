# Frozen R survival contrast solver before adding the C backend.
original_contrast_survival <- function(ynew, dnew, trt, x.cate, tau0, weightsurv,
                          ps, f1.predictor, f0.predictor, error.maxNR = 1e-3, max.iterNR = 100, tune = c(0.5, 2)) {
  x.aug <- cbind(1, x.cate)
  p.aug <- length(x.aug[1,])

  weightc <- weightsurv / mean(weightsurv)

  y <- 1 - ynew / tau0
  f1.predictor <- 1 - f1.predictor / tau0
  f0.predictor <- 1 - f0.predictor / tau0

  # Estimate coefficients with NR algorithm
  beta <- rep(0, p.aug)
  mae <- Inf
  iter <- 0
  eta.max <- 0
  while (mae > error.maxNR && iter <= max.iterNR && eta.max < 1000 * tau0) {
    eta <- as.numeric(exp(x.aug %*% beta))
    eta.max <- max(eta)
    if (eta.max == Inf) break
    error <- (trt * (y - eta * f0.predictor / 2 - f1.predictor / 2) * (1 - ps) - (1 - trt) * (y * eta - f0.predictor * eta / 2 - f1.predictor / 2) * ps) / (eta * ps + (1 - ps))
    score <- colSums(x.aug * weightc * error)
    slopewt <- (y + f0.predictor * (trt / ps - 1) / 2 + f1.predictor * ((1 - trt) / (1 - ps) - 1) / 2) * eta * ps * (1 - ps) / (eta * ps + (1 - ps))^2
    slope <- crossprod(x.aug * weightc * slopewt, x.aug)
    if (iter == 0) {
      epsilon.min <- eigen(slope, only.values = TRUE)$value[p.aug]
      epsilon0 <- epsilon.min + epsilon.min * (epsilon.min < 0)
    }
    beta <- beta + solve(slope + diag(tune[2] * abs(epsilon0), p.aug, p.aug)) %*% score * tune[1]
    mae <- sum(abs(score))
    iter <- iter + 1
  }

  converge1 <- 1 * (mae <= error.maxNR)
  converge2 <- 0

  # If NP did not converge, solve for minimizing the L2-norm (sum of squares) of the score function to avoid inversing the slope inverse (generally slower than NP)
  if (converge1 == 0) {
    lossf <- function(beta) {
      eta <- as.numeric(exp(x.aug %*% beta))
      error <- (trt * (y - eta * f0.predictor / 2 - f1.predictor / 2) * (1 - ps) - (1 - trt) * (eta * y - eta * f0.predictor / 2 - f1.predictor / 2) * ps) / (eta * ps + (1 - ps))
      score <- colSums(x.aug * weightc * error)
      return(sum(score^2))
    }

  initial.value <- lossf(rep(0, p.aug))
  fit <- optim(rep(0, p.aug), fn = lossf)
  beta <- fit$par
  converge2 <- 1 * (fit$value < initial.value / 100)
  }

  return(list(coef = beta, converge = (converge1 + converge2 > 0)))
}
