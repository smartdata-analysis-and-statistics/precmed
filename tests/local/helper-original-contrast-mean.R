# Frozen continuous contrast regression before matrix fitting/inverse reuse.
original_contrast_mean <- function(y, x.cate, trt, ps, f1.predictor, f0.predictor){

  #ps=resultcf$ps
  #f1.predictor=resultcf$f1.predictor
  #f0.predictor=resultcf$f0.predictor
  fbar.predictor <- (f1.predictor + f0.predictor) / 2
  resid <- y - fbar.predictor

  x.aug <- cbind(1, x.cate)
  xaug.star <- x.aug * (trt + ps - 2 * trt * ps) / 2
  # (1-ps)/2 when trt = 0, ps/2 when trt = 1

  outcome <- resid * (trt - ps) # (1-ps) when trt = 1, -ps when trt = 0
  fit <- lm(outcome ~ xaug.star - 1) # omitting intercept
  beta <- fit$coef

  ###### estimate the variance of the weights in ITR score 4  sandwich estimator
  error <- fit$res
  slope <- t(xaug.star) %*% xaug.star
  #sigma=solve(slope)%*%(t(xaug.star*error^2)%*%xaug.star)%*%solve(slope)
  sigma <- MASS::ginv(slope) %*% (t(xaug.star*error^2) %*% xaug.star) %*% MASS::ginv(slope)

  return(list(coef = beta, vcov = sigma))
}
