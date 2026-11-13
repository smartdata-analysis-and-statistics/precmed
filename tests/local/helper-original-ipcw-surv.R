# Frozen pre-optimization censoring probability implementation.
original_ipcw_surv <- function(y, d, x.ipcw, yf = NULL, ipcw.method = "breslow", tau0, surv.min = 0.025) {

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

    for (i in seq(n)) {
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
