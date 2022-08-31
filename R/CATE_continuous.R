# ------------------------------------------------------------------
#
# Project: Precision Medicine MS - Comprehensive R package
#
# Purpose: CATE functions
#
# Platform: Windows
# R Version: 4.1.0
#



#' Doubly robust estimators of the coefficients in the two regression
#'
#' @param y Observed outcome; vector of size \code{n}
#' @param x.cate Matrix of \code{p} baseline covariates; dimension \code{n} by \code{p}
#' @param trt Treatment received; vector of size \code{n} units with treatment coded as 0/1
#' @param ps Estimated propensity scores for all observations; vector of size \code{n}
#' @param f.predictor Initial prediction of the outcome (expected number of relapses for one unit of exposure time) conditioned
#' on the covariates \code{x} for one treatment group \code{r}; \code{mu_r(x)}, step 1 in the two regression; vector of size \code{n}
#'
#' @return Doubly robust estimators of the regression coefficients \code{beta_r} in the doubly robust estimating equation
#' where \code{r = 0, 1} is treatment received; vector of size \code{p} + 1 (intercept included)
#' @importFrom stats lm
##onearmglmmean.dr <- function(y, x.cate, trt, ps, f.predictor) {
##  f.predictor <- as.vector(f.predictor)

##  x <- as.matrix(cbind(1, (f.predictor), x.cate))

##  withCallingHandlers({
##    fit <- glm(y ~ (f.predictor) + x.cate, family = "gaussian", weights = trt / ps)

##    beta <- fit$coef
##    yhat <- (as.matrix(x[, is.na(beta) == FALSE, drop = FALSE]) %*% beta[is.na(beta) == FALSE])

##    fit2 <- glm(yhat ~ x.cate, family = "gaussian")

##  },
##  warning = function(w) { # don't change the = to <- in withCallingHandlers
##    if (grepl("non-integer", conditionMessage(w)))
##      invokeRestart("muffleWarning") # suppress warnings in glm(): "In dpois(y, mu, log = TRUE) : non-integer x = 0.557886."
##  })

##  return(fit2$coef)
##}

onearmglmmean.dr <- function(y, x.cate, trt, ps, f.predictor){

  f.predictor <- as.vector(f.predictor)

  # first, solve the weighted estimating equation in twin regression (second estimating equation with the weights). Because the weights are trt/ps, only patients who are treated are used (or untreated if trt=1-trt)
  fit <- lm(y ~ f.predictor+x.cate, weights = trt/ps)

  # second, derive the “calibrated” prediction
  coef.one <- fit$coef

  x.aug.one <- cbind(1, f.predictor, x.cate)
  yhat <- x.aug.one[, is.na(coef.one) == FALSE] %*% coef.one[is.na(coef.one) == FALSE]

  # solve the twin regression estimating equation (without weight)
  fit2  <- lm(yhat ~ x.cate)
  beta <- fit2$coef

  return(beta)
}

#' Doubly robust estimators of the coefficients in the contrast regression
#'  as well as their covariance matrix
#'
#' Solving the estimating equation \code{bar S_n (delta) = 0}
#'
#' @param y Observed outcome; vector of size \code{n}
#' @param x.cate Matrix of \code{p.cate} baseline covariates; dimension \code{n} by \code{p.cate}
#' @param trt Treatment received; vector of size \code{n} units with treatment coded as 0/1
#' @param ps Estimated propensity scores for all observations; vector of size \code{n}
#' @param f1.predictor Initial predictions of the outcome (expected number of relapses for one unit of exposure time)
#' conditioned on the covariates \code{x} for treatment group trt = 1; \code{mu_1(x)}, step 1 in the two regression; vector of size \code{n}
#' @param f0.predictor  Initial predictions of the outcome (expected number of relapses for one unit of exposure time)
#' conditioned on the covariates \code{x} for treatment group trt = 0; \code{mu_0(x)}, step 1 in the two regression; vector of size \code{n}
#' @return coef: Doubly robust estimators of the regression coefficients \code{delta_0}; vector of size \code{p} + 1 (intercept included)
#'         vcov: Variance-covariance matrix of the estimated coefficient \code{delta_0}; matrix of size \code{p} + 1 by \code{p} + 1
#' @importFrom stats lm

twoarmglmmean.dr  <- function(y, x.cate, trt, ps, f1.predictor, f0.predictor){

  #ps=resultcf$ps
  #f1.predictor=resultcf$f1.predictor
  #f0.predictor=resultcf$f0.predictor
  fbar.predictor <- (f1.predictor + f0.predictor) / 2
  resid <- y-fbar.predictor

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



#' Estimate the CATE model using specified scoring methods
#'
#' Coefficients of the CATE estimated with boosting, linear regression, two regression, contrast regression, random forest, generalized additive model
#'
#' @param y Observed outcome; vector of size \code{n} (observations)
#' @param x.cate Matrix of \code{p.cate} baseline covariates; dimension \code{n} by \code{p.cate} (covariates in the outcome model)
#' @param x.init Matrix of \code{p.init} baseline covariates; dimension \code{n} by \code{p.init}
#' It must be specified when \code{score.method = contrastReg} or \code{twoReg}.
#' @param x.ps Matrix of \code{p.ps} baseline covariates (plus a leading column of 1 for the intercept);
#' dimension \code{n} by \code{p.ps + 1} (covariates in the propensity score model plus intercept)
#' @param trt Treatment received; vector of size \code{n} units with treatment coded as 0/1
#' @param score.method A vector of one or multiple methods to estimate the CATE score.
#' Allowed values are: \code{'boosting'}, \code{'gaussian'}, \code{'twoReg'}, \code{'contrastReg'},
#' \code{'randomForest'}, \code{'gam'}. Default specifies all 6 methods.
#' @param ps.method A character value for the method to estimate the propensity score.
#' Allowed values include one of:
#' \code{'glm'} for logistic regression with main effects only (default), or
#' \code{'lasso'} for a logistic regression with main effects and LASSO penalization on
#' two-way interactions (added to the model if interactions are not specified in \code{ps.model}).
#' Relevant only when \code{ps.model} has more than one variable.
#' @param minPS A numerical value (in [0, 1]) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A number above which estimated propensity scores should be trimmed; scalar
#' @param initial.predictor.method A character vector for the method used to get initial
#' outcome predictions conditional on the covariates in \code{cate.model}
#' in \code{score.method = 'twoReg'} and \code{'contrastReg'}. Allowed values include
#' one of \code{'gaussian'} (fastest), \code{'boosting'} (default) and \code{'gam'}.
#' @param xvar.smooth.init A vector of characters indicating the name of the variables used as
#' the smooth terms if \code{initial.predictor.method = 'gam'}. The variables must be selected
#' from the variables listed in \code{init.model}.
#' Default is \code{NULL}, which uses all variables in \code{init.model}.
#' @param xvar.smooth.score A vector of characters indicating the name of the variables used as
#' the smooth terms if \code{score.method = 'gam'}. The variables must be selected
#' from the variables listed in \code{cate.model}.
#' Default is \code{NULL}, which uses all variables in \code{cate.model}.
#' @param tree.depth A positive integer specifying the depth of individual trees in boosting
#' (usually 2-3). Used only if \code{score.method = 'boosting'} or
#' if \code{score.method = 'twoReg'} or \code{'contrastReg'} and
#' \code{initial.predictor.method = 'boosting'}. Default is \code{2}.
#' @param n.trees.rf A positive integer specifying the number of trees. Used only if
#' \code{score.method = 'randomForest'}. Default is \code{1000}.
#' @param n.trees.boosting A positive integer specifying the maximum number of trees in boosting
#' (usually 100-1000). Used only if \code{score.method = 'boosting'} or
#' if \code{score.method = 'twoReg'} or \code{'contrastReg'} and
#' \code{initial.predictor.method = 'boosting'}. Default is \code{200}.
#' @param B A positive integer specifying the number of time cross-fitting is repeated in
#' \code{score.method = 'twoReg'} and \code{'contrastReg'}. Default is \code{3}.
#' @param Kfold A positive integer specifying the number of folds (parts) used in cross-fitting
#' to partition the data in \code{score.method = 'twoReg'} and \code{'contrastReg'}.
#' Default is \code{6}.
#' @param plot.gbmperf A logical value indicating whether to plot the performance measures in
#' boosting. Used only if \code{score.method = 'boosting'} or if \code{score.method = 'twoReg'}
#' or \code{'contrastReg'} and \code{initial.predictor.method = 'boosting'}. Default is \code{TRUE}.
#' @param ... Additional arguments for \code{gbm()}
#'
#' @return Depending on what score.method is, the outputs is a combination of the following:
#'           result.boosting: Results of boosting fit and best iteration, for trt = 0 and trt = 1 separately
#'           result.gaussian: Linear regression estimator (beta1 - beta0); vector of length \code{p.cate} + 1
#'           result.twoReg: Two regression estimator (beta1 - beta0); vector of length \code{p.cate} + 1
#'           result.contrastReg: A list of the contrast regression results with 3 elements:
#'               $delta.contrastReg: Contrast regression DR estimator; vector of length \code{p.cate} + 1
#'               $sigma.contrastReg: Variance covariance matrix for delta.contrastReg; matrix of size \code{p.cate} + 1 by \code{p.cate} + 1
#'           result.randomForest: Results of random forest fit and best iteration, for trt = 0 and trt = 1 separately
#'           result.gam: Results of generalized additive model fit and best iteration, for trt = 0 and trt = 1 separately
#'           best.iter: Largest best iterations for boosting (if used)
#'           fgam: Formula applied in GAM when \code{initial.predictor.method = 'gam'}
#'           warn.fit: Warnings occurred when fitting \code{score.method}
#'           err.fit:: Errors occurred when fitting \code{score.method}

intxmean <- function(y, trt, x.cate, x.init, x.ps,
                     score.method = c("boosting", "gaussian", "twoReg", "contrastReg", "gam", "randomForest"),
                     ps.method = "glm", minPS = 0.01, maxPS = 0.99,
                     initial.predictor.method = "boosting",
                     xvar.smooth.init, xvar.smooth.score,
                     tree.depth = 2, n.trees.rf = 1000, n.trees.boosting = 200, B = 1, Kfold = 2, plot.gbmperf = TRUE, ...) {

  result <- vector("list", length(score.method) + 1)
  names(result) <- c(paste0("result.", score.method), "best.iter")

  N1 <- sum(trt)
  N0 <- sum(1 - trt)
  N <- N1 + N0
  p.aug <- ncol(x.cate) + 1

  #datatot <- data.frame(y, x.cate, time)
  #colnames(datatot) <- c("y", colnames(x.cate), "time")

  datatot <- data.frame(y, x.cate)
  colnames(datatot) <- c("y", colnames(x.cate))

  datatot.init <- data.frame(y, x.init)
  colnames(datatot.init) <- c("y", colnames(x.init))

  ######### cross-fitting  ---------------------------------------------------------------

  index1 <- rep(1:Kfold, floor(N1 / Kfold))
  if (N1 > Kfold * floor(N1 / Kfold)) index1 <- c(index1, 1:(N1 - Kfold * floor(N1 / Kfold)))

  index0 <- rep(1:Kfold, floor(N0 / Kfold))
  if (N0 > Kfold * floor(N0 / Kfold)) index0 <- c(index0, Kfold + 1 - 1:(N0 - Kfold * floor(N0 / Kfold)))

  delta.twoReg.mat <- delta.contrastReg.mat <- matrix(NA, B, p.aug)
  sigma.contrastReg.mat <- matrix(0, p.aug, p.aug)
  converge <- rep(NA, B)
  best.iter <- 0
  fgam.init <- NULL

  #warn.high <- err.high <- vector("list", length = n.subgroup)
  #names(warn.high) <- names(err.high) <- paste("prop", round(prop, 2))

  warn.fit <- err.fit <- list()

  if (any(c("twoReg", "contrastReg") %in% score.method)) {
    for (bb in 1:B) {
      index1cv <- sample(x = index1, size = N1, replace = FALSE)
      index0cv <- sample(x = index0, size = N0, replace = FALSE)
      index <- rep(NA, N)
      index[trt == 1] <- index1cv
      index[trt == 0] <- index0cv

      f1.predictcv <- f0.predictcv <- pscv <- rep(NA, N)
      for (k in 1:Kfold) {
        datatot_train <- datatot[index != k, ]
        datatot_train_init <- datatot.init[index != k, ]
        x_ps_train <- x.ps[index != k, , drop = FALSE]
        trt_train <- trt[index != k]

        datatot_valid <- datatot[index == k, ]
        datatot_valid_init <- datatot.init[index == k,]
        x_ps_valid <- x.ps[index == k, , drop = FALSE]

        data1 <- datatot_train[trt_train == 1, ]
        data0 <- datatot_train[trt_train == 0, ]

        data1.init <- datatot_train_init[trt_train == 1, ]
        data0.init <- datatot_train_init[trt_train == 0, ]

        if (initial.predictor.method == "boosting") {

          # TODO: if model has a single predictor, GBM must have cv.folds = 0 https://github.com/zoonproject/zoon/issues/130
          ## Removed offset
          fit1.boosting <- gbm(y ~ ., data = data1.init, distribution = "gaussian",
                               interaction.depth = tree.depth, n.trees = n.trees.boosting, cv.folds = 5)
          best1.iter <- max(10, gbm.perf(fit1.boosting, method = "cv", plot.it = plot.gbmperf))
          withCallingHandlers({
            f1.predictcv[index == k] <- predict(object = fit1.boosting, newdata = datatot_valid_init, n.trees = best1.iter, type = "response")
          },
          ## Comment out here.
          warning = function(w) {
            if (grepl("does not add the offset", conditionMessage(w)))
              invokeRestart("muffleWarning")  # suppress warning: "predict.gbm does not add the offset to the predicted values."
          })

          fit0.boosting <- gbm(y ~ ., data = data0.init, distribution = "gaussian",
                               interaction.depth = tree.depth, n.trees = n.trees.boosting, cv.folds = 5)
          best0.iter <- max(10, gbm.perf(fit0.boosting, method = "cv", plot.it = plot.gbmperf))
          withCallingHandlers({
            f0.predictcv[index == k] <- predict(object = fit0.boosting, newdata = datatot_valid_init, n.trees = best0.iter, type = "response")
          },
          warning = function(w) {
            if (grepl("does not add the offset", conditionMessage(w)))
              invokeRestart("muffleWarning")  # suppress warning: "predict.gbm does not add the offset to the predicted values."
          })

          best.iter <- max(best.iter, best1.iter, best0.iter)

          ##Changed poisson to gaussian
        } else if (initial.predictor.method == "gaussian") {

          fit1.gaus <- glm(y ~ ., data = data1.init, family = "gaussian")

          datatot_valid_new <- datatot_valid_init
          datatot_valid_new[, "y"] <- rep(1, nrow(datatot_valid_new))
          names(datatot_valid_new[which(names(datatot_valid_new) == "y")]) <- "intercept"

          f1.predictcv[index == k] <- as.matrix(datatot_valid_new[, is.na(fit1.gaus$coefficients) == FALSE]) %*% fit1.gaus$coefficients[is.na(fit1.gaus$coefficients) == FALSE]

          fit0.gaus <- glm(y ~ ., data = data0.init, family = "gaussian")

          datatot_valid_new <- datatot_valid_init
          datatot_valid_new[, "y"] <- rep(1, nrow(datatot_valid_new))
          names(datatot_valid_new[which(names(datatot_valid_new) == "y")]) <- "intercept"

          f0.predictcv[index == k] <- as.matrix(datatot_valid_new[, is.na(fit0.gaus$coefficients) == FALSE]) %*% fit0.gaus$coefficients[is.na(fit0.gaus$coefficients) == FALSE]


        } else if (initial.predictor.method == "gam") {

          xvars <- colnames(x.init)
          if (is.null(xvar.smooth.init) == TRUE){
            fgam.init <- paste0("y ~ ", paste0("s(", xvars, ")", collapse = "+"))
          } else {
            xvar.smooth2.init <- xvars[stringr::str_detect(xvars, paste(paste0(xvar.smooth.init, "$"), collapse = "|"))] # conver to the preprocessed names
            xvar.linear.init <- setdiff(xvars, xvar.smooth2.init) # the remaining xvars in x.cate but not in xvar.smooth are linear predictors
            fgam.init <- paste0("y ~ ", paste0(xvar.linear.init, collapse = "+"), "+", paste0("s(", xvar.smooth2.init, ")", collapse = "+"))
          }
          ## changed poisson to gaussian, removed offset option
          fit1.gam.init <- mgcv::gam(as.formula(fgam.init), data = data1.init, family = "gaussian")
          f1.predictcv[index == k] <- predict(object = fit1.gam.init, newdata = datatot_valid, type = "response")

          fit0.gam.init <- mgcv::gam(as.formula(fgam.init), data = data0.init, family = "gaussian")
          f0.predictcv[index == k] <- predict(object = fit0.gam.init, newdata = datatot_valid, type = "response")
        }

        if (ps.method == "glm") {
          pscv[index == k] <- glm.simplereg.ps(x.ps = x_ps_train, trt = trt_train, xnew = x_ps_valid, minPS = minPS, maxPS = maxPS)
        } else {
          pscv[index == k] <- glm.ps(x.ps = x_ps_train, trt = trt_train, xnew = x_ps_valid, minPS = minPS, maxPS = maxPS)
        }

      }#end of Kfold loops

      if ("twoReg" %in% score.method) {
        ## bb-th cross fitting two regression estimator
        ## Replace count function to mean
        beta1.final <- onearmglmmean.dr(y = y, x.cate = x.cate, trt = trt, ps = pscv, f.predictor = f1.predictcv)
        beta0.final <- onearmglmmean.dr(y = y, x.cate = x.cate, trt = 1 - trt, ps = 1 - pscv, f.predictor = f0.predictcv)
        delta.twoReg.mat[bb, ] <- as.vector(beta1.final - beta0.final)
      }#end of if ("twoReg" %in% score.method)

      ##      if ("contrastReg" %in% score.method) {
      ##        ## bb-th cross fitting contrast regression estimator
      ##        fit_two <- twoarmglmmean.dr(y = y, x.cate = x.cate, trt = trt, ps = pscv,
      ##                                    f1.predictor = f1.predictcv, f0.predictor = f0.predictcv,
      ##                                    error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune)
      ##        delta.contrastReg.mat[bb, ] <- fit_two$coef
      ##        converge[bb] <- fit_two$converge
      ##        if(converge[bb] == TRUE) sigma.contrastReg.mat <- sigma.contrastReg.mat + fit_two$vcov
      ##      }#end of if ("contrastReg" %in% score.method)


      if ("contrastReg" %in% score.method) {
        ## bb-th cross fitting contrast regression estimator
        fit_two <- twoarmglmmean.dr(y = y, x.cate = x.cate, trt = trt, ps = pscv,
                                    f1.predictor = f1.predictcv, f0.predictor = f0.predictcv)
        delta.contrastReg.mat[bb, ] <- fit_two$coef
        #        converge[bb] <- fit_two$converge
        #        if(converge[bb] == TRUE) sigma.contrastReg.mat <- sigma.contrastReg.mat + fit_two$vcov
        sigma.contrastReg.mat <- sigma.contrastReg.mat + fit_two$vcov
      }#end of if ("contrastReg" %in% score.method)


    }#end of B loops
  }#end of if c("twoReg", "contrastReg") %in% score.method

  if ("boosting" %in% score.method) {
    ## Boosting method based on the entire data (score 1)
    data1.boosting <- datatot[trt == 1, ]
    # TODO: if model has a single predictor, GBM must have cv.folds = 0 https://github.com/zoonproject/zoon/issues/130
    ##Changed poisson to gaussian, removed offset term

    fit1.boosting <- meanCatch(gbm(y ~ ., data = data1.boosting, distribution = "gaussian",
                                   interaction.depth = tree.depth, n.trees = n.trees.boosting, cv.folds = 5))
    if (length(fit1.boosting$errors) > 0){
      best1.iter <- NA
    } else{
      best1.iter <- max(10, gbm.perf(fit1.boosting$fit, method = "cv", plot.it = plot.gbmperf))
    }

    data0.boosting <- datatot[trt == 0, ]
    fit0.boosting <- meanCatch( gbm(y ~ ., data = data0.boosting, distribution = "gaussian",
                                    interaction.depth = tree.depth, n.trees = n.trees.boosting, cv.folds = 5))

    if (length(fit0.boosting$errors) > 0){
      best0.iter <- NA
    } else{
      best0.iter <- max(10, gbm.perf(fit0.boosting$fit, method = "cv", plot.it = plot.gbmperf))
    }

    result$result.boosting <- list(fit1.boosting = fit1.boosting$fit, fit0.boosting = fit0.boosting$fit)

    temp.iter <- c(best.iter, best1.iter, best0.iter)
    temp.iter <- temp.iter[!is.na(temp.iter)]
    if (length(temp.iter) == 0){best.iter <- NA} else {best.iter <- max(temp.iter)}

    warn.fit[["boosting"]] <- append(fit1.boosting$warnings, fit0.boosting$warnings)
    err.fit[["boosting"]] <- append(fit1.boosting$errors, fit0.boosting$errors)


  }
  ##changed poisson to gaussian
  if ("gaussian" %in% score.method) {
    ## Naive Poisson regression method  (score 2)
    beta1.ini <- glm(y ~ x.cate, family = "gaussian", subset = (trt == 1))$coef
    beta0.ini <- glm(y ~ x.cate, family = "gaussian", subset = (trt == 0))$coef
    delta.gaussian <- beta1.ini - beta0.ini
    names(delta.gaussian) <- c("(Intercept)", colnames(x.cate))
    result$result.gaussian <- delta.gaussian
  }


  ##TODO: Try to give message when there are not enough observations
  if ("gam" %in% score.method) {

    xvars <- colnames(x.cate)
    if (is.null(xvar.smooth.score) == TRUE){
      fgam.score <- paste0("y ~ ", paste0("s(", xvars, ")", collapse = "+"))
    } else {
      xvar.smooth2.score <- xvars[stringr::str_detect(xvars, paste(paste0(xvar.smooth.score, "$"), collapse = "|"))] # conver to the preprocessed names
      xvar.linear.score <- setdiff(xvars, xvar.smooth2.score) # the remaining xvars in x.cate but not in xvar.smooth are linear predictors
      fgam.score <- paste0("y ~ ", paste0(xvar.linear.score, collapse = "+"), "+", paste0("s(", xvar.smooth2.score, ")", collapse = "+"))
    }

    data1.gam.cate<- datatot[trt == 1, ]
    data0.gam.cate <- datatot[trt == 0, ]

    fit1.gam.cate <- meanCatch( mgcv::gam(as.formula(fgam.score), data = data1.gam.cate, family = "gaussian"))
    fit0.gam.cate <- meanCatch( mgcv::gam(as.formula(fgam.score), data = data0.gam.cate, family = "gaussian"))

    result$result.gam <- list(fit1.gam = fit1.gam.cate$fit, fit0.gam = fit0.gam.cate$fit)

    warn.fit[["gam"]] <- append(fit1.gam.cate$warnings, fit0.gam.cate$warnings)
    err.fit[["gam"]] <- append(fit1.gam.cate$errors, fit0.gam.cate$errors)

  }


  if ("randomForest" %in% score.method) {

    data1.rf <- datatot[trt == 1, ]
    data0.rf <- datatot[trt == 0, ]

    fit1.rf <- meanCatch(randomForestSRC::rfsrc(y ~ ., data = data1.rf, ntree = n.trees.rf))
    fit0.rf <- meanCatch(randomForestSRC::rfsrc(y ~ ., data = data0.rf, ntree = n.trees.rf))

    result$result.randomForest <- list(fit1.rf = fit1.rf$fit, fit0.rf = fit0.rf$fit)

    warn.fit[["randomForest"]] <- append(fit1.rf$warnings, fit0.rf$warnings)
    err.fit[["randomForest"]] <- append(fit1.rf$errors, fit0.rf$errors)
  }

  if ("twoReg" %in% score.method) {
    ## Final two regression estimator (score 3)
    delta.twoReg <- colMeans(delta.twoReg.mat)
    names(delta.twoReg) <- c("(Intercept)", colnames(x.cate))
    result$result.twoReg <- delta.twoReg
  }

  if ("contrastReg" %in% score.method) {
    ## Final contrast regression estimator (score 4)
    #    converge.contrastReg <- (sum(converge) > 0)
    #    if(converge.contrastReg == TRUE){
    #      delta.contrastReg <- colMeans(delta.contrastReg.mat[converge == TRUE, , drop = FALSE])
    #      sigma.contrastReg <- sigma.contrastReg.mat/sum(converge)
    #    } else {
    #      delta.contrastReg <- colMeans(delta.contrastReg.mat)
    #      sigma.contrastReg <- sigma.contrastReg.mat
    #    }
    delta.contrastReg <- colMeans(delta.contrastReg.mat[, , drop = FALSE])
    sigma.contrastReg <- sigma.contrastReg.mat/B

    names(delta.contrastReg) <- colnames(sigma.contrastReg) <- rownames(sigma.contrastReg) <- c("(Intercept)", colnames(x.cate))
    result$result.contrastReg <- list(delta.contrastReg = delta.contrastReg,
                                      sigma.contrastReg = sigma.contrastReg
                                      #                                      #converge.contrastReg = converge.contrastReg
    )
  }

  result$best.iter <- best.iter
  result$fgam <- fgam.init
  result$warn.fit <- warn.fit
  result$err.fit <- err.fit

  return(result)
}



#' Calculate the CATE score given the baseline covariates for specified scoring method methods
#'
#' Based on intxmean results of the CATE coefficients estimated with boosting, linear regression, two regression, contrast regression, random forest, generalized additive model
#'
#' @param fit List of objects generated from intxmean: outputs of boosting, linear regression, two regression, contrast regression, random forest, generalized additive model
#' @param x.cate Matrix of \code{p.cate} baseline covariates; dimension \code{n} (observations) by \code{p.cate} (covariates in the outcome model)
#' @param score.method A vector of one or multiple methods to estimate the CATE score.
#' Allowed values are: \code{'boosting'}, \code{'gaussian'}, \code{'twoReg'}, \code{'contrastReg'}, \code{'randomForest'}, \code{'gam'}.
#' Default specifies all 6 methods.
#'
#' @return score.boosting: Estimated CATE score for all \code{n} observations with the boosting method; vector of size \code{n}
#'         score.gaussian: Estimated CATE score for all \code{n} observations with the linear regression method; vector of size \code{n}
#'         score.twoReg: Estimated CATE score for all \code{n} observations with the two regression method; vector of size \code{n}
#'         score.contrastReg: Estimated CATE score for all \code{n} observations with the contrast regression method; vector of size \code{n}
#'         score.randomForest: Estimated CATE score for all \code{n} observations with the random forest method; vector of size \code{n}
#'         score.gam: Estimated CATE score for all \code{n} observations with the generalized additive model; vector of size \code{n}
#'         score = NA if the corresponding method is not called
scoremean <- function(fit, x.cate,
                      score.method = c("boosting", "gaussian", "twoReg", "contrastReg", "randomForest", "gam")) {

  result <- vector("list", length(score.method))
  names(result) <- paste0("score.", score.method)
  x.aug <- cbind(1, x.cate)

  if ("boosting" %in% score.method) {
    fit0.boosting <- fit$result.boosting$fit0.boosting
    best0.iter <- fit$result.boosting$best0.iter
    fit1.boosting <- fit$result.boosting$fit1.boosting
    best1.iter <- fit$result.boosting$best1.iter

    datanew <- data.frame(x = x.cate)
    colnames(datanew) <- c(colnames(x.cate))
    suppressMessages({
      predict0 <- predict(object = fit0.boosting, newdata = datanew, n.trees = best0.iter)
      predict1 <- predict(object = fit1.boosting, newdata = datanew, n.trees = best1.iter)
    },
    classes = "message")
    result$score.boosting <- predict1 - predict0
  }

  if ("gaussian" %in% score.method) {
    delta.gaussian <- fit$result.gaussian
    result$score.gaussian <- as.numeric(as.matrix(x.aug[,is.na(delta.gaussian) == FALSE]) %*% delta.gaussian[is.na(delta.gaussian) == FALSE])
  }

  if ("twoReg" %in% score.method) {
    delta.twoReg <- fit$result.twoReg
    result$score.twoReg <- as.numeric(as.matrix(x.aug[,is.na(delta.twoReg) == FALSE]) %*% delta.twoReg[is.na(delta.twoReg) == FALSE])
  }

  if ("contrastReg" %in% score.method) {
    delta.contrastReg <- fit$result.contrastReg$delta.contrastReg
    result$score.contrastReg <- as.numeric(as.matrix(x.aug[,is.na(delta.contrastReg) == FALSE]) %*% delta.contrastReg[is.na(delta.contrastReg) == FALSE])
  }

  if ("gam" %in% score.method) {

    fit0.gam <- fit$result.gam$fit0.gam
    fit1.gam <- fit$result.gam$fit1.gam

    datanew <- data.frame(x = x.cate)
    colnames(datanew) <- c(colnames(x.cate))

    predict0 <- predict(object = fit0.gam, newdata = datanew)
    predict1 <- predict(object = fit1.gam, newdata = datanew)

    result$score.gam <- predict1 - predict0
  }

  if ("randomForest" %in% score.method) {

    fit0.rf <- fit$result.randomForest$fit0.rf
    fit1.rf <- fit$result.randomForest$fit1.rf

    datanew <- data.frame(x = x.cate)
    colnames(datanew) <- c(colnames(x.cate))

    predict0 <- predict(object = fit0.rf, newdata = datanew)$predicted
    predict1 <- predict(object = fit1.rf, newdata = datanew)$predicted

    result$score.randomForest <- predict1 - predict0
  }

  return(result)
}
