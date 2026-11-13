# Frozen continuous-outcome balanced splitting before matrix reuse.
original_balancemean_split <- function(y, trt, x.cate, x.ps,
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

    y.train <- y[-id.valid]
    trt.train <- trt[-id.valid]
    x.train <- x[-id.valid, , drop = FALSE]
    x.cate.train <- x.cate[-id.valid, , drop = FALSE]
    x.ps.train <- x.ps[-id.valid, , drop = FALSE]

    diffx <- colMeans(x.train) - colMeans(x.valid)
    errorx <- max(abs(diffx / sdx))

    #Replace drcount to drmean
    est.valid  <- drmean(y = y.valid, x.cate = x.cate.valid, x.ps = x.ps.valid, trt = trt.valid,
                         ps.method = "glm", minPS = minPS, maxPS = maxPS, interactions = TRUE)$mean.diff
    est.train <- drmean(y = y.train, x.cate = x.train, x.ps = x.ps.train, trt = trt.train,
                        ps.method = "glm", minPS = minPS, maxPS = maxPS, interactions = TRUE)$mean.diff

    #Change error?
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

    y.train <- y[-bestid.valid]
    trt.train <- trt[-bestid.valid]
    x.cate.train <- x.cate[-bestid.valid, , drop = FALSE]
    x.ps.train <- x.ps[-bestid.valid, , drop = FALSE]

    warning(paste("Maximum iteration reached and the SMD between training and validation set is still greater than error.max (error=", round(errorbest, 4), "). Consider increasing max.iter, decreasing error.max, or increasing sample size.", sep = ""))
  }

  return(list(y.train = y.train, trt.train = trt.train, x.cate.train = x.cate.train, x.ps.train = x.ps.train,
              y.valid = y.valid, trt.valid = trt.valid, x.cate.valid = x.cate.valid, x.ps.valid = x.ps.valid, bestid.valid = bestid.valid))
}
