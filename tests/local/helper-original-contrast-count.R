# Frozen R solver before introducing the optional C backend.
original_contrast_count <- function(y, x.cate, time, trt,
                              ps, f1.predictor, f0.predictor,
                              error.maxNR = 1e-3, max.iterNR = 150, tune = c(0.5, 2.0)) {

  y <- y * exp(-time)
  x.aug <- cbind(1, x.cate)
  p.aug <- length(x.aug[1, ])

  # Estimate coefficients with NR algorithm
  beta <- rep(0, p.aug)
  mae <- Inf
  iter <- 0
  epsilon.min <- Inf
  while (mae > error.maxNR && iter <= max.iterNR && epsilon.min > 0) {
    eta <- as.numeric(exp(x.aug %*% beta))
    error <- (trt * (y - eta * f0.predictor / 2 - f1.predictor / 2) * (1 - ps) - (1 - trt) * (eta * y - eta * f0.predictor / 2 - f1.predictor / 2) * ps) / (eta * ps + (1 - ps))
    score <- colSums(x.aug * error)
    slopewt <- (y + f0.predictor * (trt / ps - 1) / 2 + f1.predictor * ((1 - trt) / (1 - ps) - 1) / 2) * eta * ps * (1 - ps) / (eta * ps + (1 - ps))^2
    slope <- t(x.aug * slopewt) %*% x.aug
    epsilon.min <- eigen(slope, only.values = TRUE)$value[p.aug]
    if (iter == 0) epsilon0 <- epsilon.min + epsilon.min * (epsilon.min < 0) # fixed to all iterations
    beta <- beta + solve(slope + diag(tune[2] * abs(epsilon0), p.aug, p.aug)) %*% score * tune[1] # adding the diagonal matrix to slove potential singualrity issues
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
      score <- colSums(x.aug * error)
      return(sum(score^2))
    }

    initial.value <- lossf(rep(0, p.aug))
    fit <- optim(rep(0, p.aug), fn = lossf) #, control=list(trace=T))
    beta <- fit$par
    converge2 <- 1 * (fit$value < initial.value / 100)
  }

  beta <- as.vector(beta)
  eta <- as.numeric(exp(x.aug %*% beta))
  error <- (trt * (y - eta * f0.predictor / 2 - f1.predictor / 2) * (1 - ps) - (1 - trt) * (eta * y - eta * f0.predictor / 2 - f1.predictor / 2) * ps) / (eta * ps + (1 - ps))
  score <- colSums(x.aug * error)
  slopewt <- (y + f0.predictor * (trt / ps - 1) / 2 + f1.predictor * ((1 - trt) / (1 - ps) - 1) / 2) * eta * ps * (1 - ps) / (eta * ps + (1 - ps))^2
  slope <- t(x.aug * slopewt) %*% x.aug
  sigma <- solve(slope) %*% (t(x.aug * error^2) %*% x.aug) %*% solve(slope)

  return(list(coef = beta, vcov = sigma, converge = (converge1 + converge2 > 0)))
}
