# Frozen pre-optimization implementation for the local regression test.
original_intxcount <- function(y, trt, x.cate, x.ps, time,
                      score.method = c("boosting", "poisson", "twoReg", "contrastReg", "negBin"),
                      ps.method = "glm", minPS = 0.01, maxPS = 0.99,
                      initial.predictor.method = "boosting",
                      xvar.smooth = NULL,
                      tree.depth = 2, n.trees.boosting = 200, B = 3, Kfold = 6, plot.gbmperf = TRUE,
                      error.maxNR = 1e-3, max.iterNR = 150, tune = c(0.5, 2.0), ...) {

  result <- vector("list", length(score.method) + 1)
  names(result) <- c(paste0("result.", score.method), "best.iter")

  # Calculate the total number of patients who received the treatment (N1)
  # and the total number of patients who did not receive the treatment (N0)
  N1 <- sum(trt)               # N1: number of treated patients (trt = 1)
  N0 <- length(trt) - N1        # N0: number of untreated patients (trt = 0)

  # Total number of patients
  N <- length(trt)              # N is simply the total number of patients

  # Calculate p.aug: the number of covariates in x.cate plus 1 for the intercept
  p.aug <- ncol(x.cate) + 1

  datatot <- data.frame("y" = y, x.cate, "time" = time)

  ######### cross-fitting

  # Generate indices for N1 (normal order) and N0 (reverse order)
  index1 <- generate_kfold_indices(N1, Kfold, reverse = FALSE)
  index0 <- generate_kfold_indices(N0, Kfold, reverse = TRUE)


  delta.twoReg.mat <- delta.contrastReg.mat <- matrix(NA, B, p.aug)
  sigma.contrastReg.mat <- matrix(0, p.aug, p.aug)
  converge <- rep(NA, B)
  best.iter <- 0
  fgam <- NULL

  if (any(c("twoReg", "contrastReg") %in% score.method)) {
    for (bb in seq(B)) {
      index1cv <- sample(x = index1, size = N1, replace = FALSE)
      index0cv <- sample(x = index0, size = N0, replace = FALSE)
      index <- rep(NA, N)
      index[trt == 1] <- index1cv
      index[trt == 0] <- index0cv

      f1.predictcv <- f0.predictcv <- pscv <- rep(NA, N)
      for (k in seq(Kfold)) {
        datatot_train <- datatot[index != k, ]
        x_ps_train <- x.ps[index != k, , drop = FALSE]
        trt_train <- trt[index != k]

        datatot_valid <- datatot[index == k, ]
        x_ps_valid <- x.ps[index == k, , drop = FALSE]

        data1 <- datatot_train[trt_train == 1, ]
        data0 <- datatot_train[trt_train == 0, ]

        if (initial.predictor.method == "boosting") {
          # if model has a single predictor, GBM must have cv.folds = 0 https://github.com/zoonproject/zoon/issues/130
          cate.cvfold <- ifelse(ncol(x.cate) == 1, 0, 5)
          fit1.boosting <- gbm(y ~ . - time + offset(time),
                               data = data1,
                               distribution = "poisson",
                               interaction.depth = tree.depth,
                               n.trees = n.trees.boosting,
                               cv.folds = cate.cvfold, ...)
          best1.iter <- max(10, gbm.perf(fit1.boosting, method = "cv", plot.it = plot.gbmperf))


          fit1.boosting <- gbm(y ~ . - time + offset(time), data = data1, distribution = "poisson",
                               interaction.depth = tree.depth, n.trees = n.trees.boosting, cv.folds = cate.cvfold, ...)
          best1.iter <- max(10, gbm.perf(fit1.boosting, method = "cv", plot.it = plot.gbmperf))
          withCallingHandlers({
            f1.predictcv[index == k] <- predict(object = fit1.boosting, newdata = datatot_valid, n.trees = best1.iter, type = "response")
          },
          warning = function(w) {
            if (grepl("does not add the offset", conditionMessage(w)))
              invokeRestart("muffleWarning")  # suppress warning: "predict.gbm does not add the offset to the predicted values."
          })

          fit0.boosting <- gbm(y ~ . - time + offset(time), data = data0, distribution = "poisson",
                               interaction.depth = tree.depth, n.trees = n.trees.boosting, cv.folds = cate.cvfold, ...)
          best0.iter <- max(10, gbm.perf(fit0.boosting, method = "cv", plot.it = plot.gbmperf))
          withCallingHandlers({
            f0.predictcv[index == k] <- predict(object = fit0.boosting, newdata = datatot_valid, n.trees = best0.iter, type = "response")
          },
          warning = function(w) {
            if (grepl("does not add the offset", conditionMessage(w)))
              invokeRestart("muffleWarning")  # suppress warning: "predict.gbm does not add the offset to the predicted values."
          })

          best.iter <- max(best.iter, best1.iter, best0.iter)


        } else if (initial.predictor.method == "poisson") {

          fit1.pois <- glm(y ~ . - time + offset(time), data = data1, family = "poisson")
          f1.predictcv[index == k] <- predict(object = fit1.pois, newdata = datatot_valid, type = "response")

          fit0.pois <- glm(y ~ . - time + offset(time), data = data0, family = "poisson")
          f0.predictcv[index == k] <- predict(object = fit0.pois, newdata = datatot_valid, type = "response")

        } else if (initial.predictor.method == "gam") {

          xvars <- colnames(x.cate)
          if (is.null(xvar.smooth)) {
            fgam <- paste0("y ~ ", paste0("s(", xvars, ")", collapse = "+"))
          } else {
            xvar.smooth2 <- xvars[str_detect(xvars, paste(paste0(xvar.smooth, "$"), collapse = "|"))] # conver to the preprocessed names
            xvar.linear <- setdiff(xvars, xvar.smooth2) # the remaining xvars in x.cate but not in xvar.smooth are linear predictors
            fgam <- paste0("y ~ ", paste0(xvar.linear, collapse = "+"), "+", paste0("s(", xvar.smooth2, ")", collapse = "+"))
          }
          fit1.gam <- mgcv::gam(as.formula(fgam), offset = time, data = data1, family = "poisson")
          f1.predictcv[index == k] <- predict(object = fit1.gam, newdata = datatot_valid, type = "response")

          fit0.gam <- mgcv::gam(as.formula(fgam), offset = time, data = data0, family = "poisson")
          f0.predictcv[index == k] <- predict(object = fit0.gam, newdata = datatot_valid, type = "response")
        }

        if (ps.method == "glm") {
          pscv[index == k] <- glm.simplereg.ps(x.ps = x_ps_train, trt = trt_train, xnew = x_ps_valid, minPS = minPS, maxPS = maxPS)
        } else {
          pscv[index == k] <- glm.ps(x.ps = x_ps_train, trt = trt_train, xnew = x_ps_valid, minPS = minPS, maxPS = maxPS)
        }

      }#end of Kfold loops

      if ("twoReg" %in% score.method) {
        ## bb-th cross fitting two regression estimator
        beta1.final <- onearmglmcount.dr(y = y, x.cate = x.cate, time = time, trt = trt, ps = pscv, f.predictor = f1.predictcv)
        beta0.final <- onearmglmcount.dr(y = y, x.cate = x.cate, time = time, trt = 1 - trt, ps = 1 - pscv, f.predictor = f0.predictcv)
        delta.twoReg.mat[bb, ] <- as.vector(beta1.final - beta0.final)
      }#end of if ("twoReg" %in% score.method)

      if ("contrastReg" %in% score.method) {
        ## bb-th cross fitting contrast regression estimator
        fit_two <- twoarmglmcount.dr(y = y, x.cate = x.cate, time = time, trt = trt, ps = pscv,
                                     f1.predictor = f1.predictcv, f0.predictor = f0.predictcv,
                                     error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune)
        delta.contrastReg.mat[bb, ] <- fit_two$coef
        converge[bb] <- fit_two$converge
        if (converge[bb]) sigma.contrastReg.mat <- sigma.contrastReg.mat + fit_two$vcov
      }#end of if ("contrastReg" %in% score.method)
    }#end of B loops
  }#end of if c("twoReg", "contrastReg") %in% score.method

  if ("boosting" %in% score.method) {
    ## Boosting method based on the entire data (score 1)
    data1 <- datatot[trt == 1, ]
    # if model has a single predictor, GBM must have cv.folds = 0 https://github.com/zoonproject/zoon/issues/130
    cate.cvfold <- ifelse(ncol(x.cate) == 1, 0, 5)
    cate.gbmmethod <- ifelse(ncol(x.cate) == 1, "OOB", "cv")
    if (cate.gbmmethod == "OOB") warning("If the model has a single predictor, GBM must use OOB")

    fit1.boosting <- gbm(y ~ . - time + offset(time), data = data1, distribution = "poisson",
                           interaction.depth = tree.depth, n.trees = n.trees.boosting, cv.folds = cate.cvfold, ...)
    best1.iter <- max(10, gbm.perf(fit1.boosting, method = cate.gbmmethod, plot.it = plot.gbmperf))

    data0 <- datatot[trt == 0, ]
    fit0.boosting <- gbm(y ~ . - time + offset(time), data = data0, distribution = "poisson",
                           interaction.depth = tree.depth, n.trees = n.trees.boosting, cv.folds = cate.cvfold, ...)
    best0.iter <- max(10, gbm.perf(fit0.boosting, method = cate.gbmmethod, plot.it = plot.gbmperf))

    result$result.boosting <- list(fit0.boosting = fit0.boosting, best0.iter = best0.iter, fit1.boosting = fit1.boosting, best1.iter = best1.iter)
    best.iter <- max(best.iter, best1.iter, best0.iter)
  }

  if ("poisson" %in% score.method) {
    ## Naive Poisson regression method  (score 2)
    beta1.ini <- glm(y ~ x.cate + offset(time), family = "poisson", subset = (trt == 1))$coef
    beta0.ini <- glm(y ~ x.cate + offset(time), family = "poisson", subset = (trt == 0))$coef
    delta.poisson <- beta1.ini - beta0.ini
    names(delta.poisson) <- c("(Intercept)", colnames(x.cate))
    result$result.poisson <- delta.poisson
  }

  if ("twoReg" %in% score.method) {
    ## Final two regression estimator (score 3)
    delta.twoReg <- colMeans(delta.twoReg.mat)
    names(delta.twoReg) <- c("(Intercept)", colnames(x.cate))
    result$result.twoReg <- delta.twoReg
  }

  if ("contrastReg" %in% score.method) {
    ## Final contrast regression estimator (score 4)
    converge.contrastReg <- (sum(converge) > 0)
    if (converge.contrastReg) {
      delta.contrastReg <- colMeans(delta.contrastReg.mat[converge == TRUE, , drop = FALSE])
      sigma.contrastReg <- sigma.contrastReg.mat/sum(converge)
    } else {
      delta.contrastReg <- colMeans(delta.contrastReg.mat)
      sigma.contrastReg <- sigma.contrastReg.mat
    }
    names(delta.contrastReg) <- colnames(sigma.contrastReg) <- rownames(sigma.contrastReg) <- c("(Intercept)", colnames(x.cate))
    result$result.contrastReg <- list(delta.contrastReg = delta.contrastReg,
                                      sigma.contrastReg = sigma.contrastReg,
                                      converge.contrastReg = converge.contrastReg)
  }


  if ("negBin" %in% score.method) {
    ## Naive negative binomial regression method  (score 5)
    beta1.ini <- glm.nb(y ~ x.cate + offset(time), subset = (trt == 1), maxit = 500)$coef
    beta0.ini <- glm.nb(y ~ x.cate + offset(time), subset = (trt == 0), maxit = 500)$coef
    delta.negBin <- beta1.ini - beta0.ini
    names(delta.negBin) <- c("(Intercept)", colnames(x.cate))
    result$result.negBin <- delta.negBin
  }

  result$best.iter <- best.iter
  result$fgam <- fgam

  return(result)
}
