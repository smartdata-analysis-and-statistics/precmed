# ------------------------------------------------------------------
#
# Project: Precision Medicine MS - Comprehensive R package
#
# Purpose: CATE functions for Survival outcomes
#
# Platform: Windows
# R Version: 4.1.0
#
#   Modifications:
#
#   Date			By			Description
# --------		--------	-----------------------------
#  14Jun2021  dl      Start the script
#  02Jul2021  dl      Add onearmsurv.dr() and twoarmsurv.dr()
#  07Jul2021  dl      Add intxsurv() and scoresurv()
#  08JUL2021  gs      Reviewed onearmsurv() and twoarmsurv()
#  09JUL2021  dl      Minor changes to twoarmsurv.dr() and intxsurv()
#  13JUL2021  dl      Changed x to x.cate in onearmsurv.dr(), twoarmsurv.dr(), and scoresurv();
#                     Changed argument x to x.cate and x.ps in intxsurv();
#  13JUL2021  gs      Reviewed intxsurv() and scoresurv() with minor edits
#  14JUL2021  dl      Add argument x.ipw to intxsurv()
#  23JUL2021  dl      Modify documentation
#  08AUG2021  gs      Recover from NA coefficients in fit in onearmsurv.dr()
#                     Move two and contrast reg empty vector definition inside loop
#  13AUG2021  dl      Separate ntrees.rf and ntrees.gbm; minor changes
#  09SEP2021  gs      Add best.iter output in intxsurv when initial.predictor = "boosting"
#  14SEP2021  gs      Rename x.ipw to x.ipcw
#  29SEP2021  pj      Remove argument seed.cf in intxsurv
#  05OCT2021  gs      Add ipcw.method argument to survprd -> inherited in intxsurv
#  06OCT2021  pj      Remove set.seed functions
#                     Specify gam's gam function in intxsurv not to confuse with mgcv's gam
#  13OCT2021  pj      Specify gam's s function in intxsurv
#  15OCT2021  gs      Lower case ipcw.method argument names, change survprd to ipcw.surv
#  22OCT2021  gs      Add initial.predictor logistic regression, wrap non-integer warning messages
#  03NOV2021  pj      Change gam::gam to mgcv::gam and revise formula to have smoothed term on f.predict.ini and all other remaining covariates in data
#  04FEB2022  pj      Add missing argument in documentation and change argument weight to weights
#  04MAR2022  gs      Change n.trees to n.trees.boosting
# ------------------------------------------------------------------

#' Doubly robust estimators of the coefficients in the two regression
#'
#' @param ynew Truncated survival or censoring time; vector of size \code{n}.
#' @param dnew The event indicator after truncation, \code{1 = event or censored after truncation, 0 = censored before truncation};
#' vector of size \code{n}.
#' @param trt Treatment received; vector of size \code{n} with treatment coded as 0/1.
#' @param x.cate Matrix of \code{p.cate} baseline covariates specified in the outcome model; dimension \code{n} by \code{p.cate}.
#' @param tau0 The truncation time for defining restricted mean time lost.
#' @param weightsurv Estimated inverse probability of censoring weights with truncation for all observations; vector of size \code{n}.
#' @param ps Estimated propensity scores for all observations; vector of size \code{n}
#' @param f.predictor Initial prediction of the outcome (restricted mean time loss) conditioned on the covariates \code{x.cate} for one treatment group \code{r};
#' \code{mu_r(x.cate)}, step 1 in the two regression; vector of size \code{n}
#'
#' @return Doubly robust estimators of the two regression coefficients \code{beta_r} where \code{r = 0, 1} is treatment received; vector of size \code{p.cate} + 1 (intercept included)
#'

onearmsurv.dr <- function(ynew, dnew, trt, x.cate, tau0, weightsurv, ps, f.predictor) {
  weightc <- weightsurv / mean(weightsurv)

  ynew <- tau0 - ynew
  f.predictor <- tau0 - f.predictor
  x <- as.matrix(cbind(1, log(f.predictor), x.cate))

  withCallingHandlers({
    fit <- glm(ynew ~ log(f.predictor) + x.cate, family = "poisson", weights = weightc * trt / ps)
    beta <- fit$coef
    yhat <- exp(as.matrix(x[, is.na(beta) == FALSE, drop = FALSE]) %*% beta[is.na(beta) == FALSE])
    fit2 <- glm(yhat ~ x.cate, family = "poisson")
  },
  warning = function(w) { # don't change the = to <- in withCallingHandlers
    if (grepl("non-integer", conditionMessage(w)))
      invokeRestart("muffleWarning") # suppress warnings in glm(): "In dpois(y, mu, log = TRUE) : non-integer x = 0.557886."
  })
  return(fit2$coef)
}




#' Doubly robust estimators of the coefficients in the contrast regression
#'  as well as their covariance matrix and convergence information
#'
#' Newton-Raphson algorithm is used to solve the estimating equation \code{bar S_n (delta) = 0}
#'
#' @param ynew Truncated survival time; vector of size \code{n}
#' @param dnew Event indicator after truncation; vector of size \code{n}
#' @param trt Treatment received; vector of size \code{n} with treatment coded as 0/1.
#' @param x.cate Matrix of \code{p.cate} baseline covariates specified in the outcome model; dimension \code{n} by \code{p.cate}.
#' @param tau0 The truncation time for defining restricted mean time lost.
#' @param weightsurv Estimated inverse probability of censoring weights with truncation for all observations; vector of size \code{n}.
#' @param ps Estimated propensity scores for all observations; vector of size \code{n}
#' @param f1.predictor Initial predictions of the outcome (restricted mean time loss) conditioned on the covariates \code{x.cate} for treatment group trt = 1;
#' \code{mu_1(x.cate)}, step 1 in the two regression; vector of size \code{n}
#' @param f0.predictor  Initial predictions of the outcome (restricted mean time loss) conditioned on the covariates \code{x.cate} for treatment group trt = 0;
#' \code{mu_0(x.cate)}, step 1 in the two regression; vector of size \code{n}
#' @param error.maxNR A numerical value > 0 indicating the minimum value of the mean absolute
#' error in Newton Raphson algorithm. Used only if \code{score.method = 'contrastReg'}.
#' Default is \code{0.001}.
#' @param max.iterNR A positive integer indicating the maximum number of iterations in the
#' Newton Raphson algorithm. Used only if \code{score.method = 'contrastReg'}.
#' Default is \code{100}.
#' @param tune A vector of 2 numerical values > 0 specifying tuning parameters for the
#' Newton Raphson algorithm. \code{tune[1]} is the step size, \code{tune[2]} specifies a
#' quantity to be added to diagonal of the slope matrix to prevent singularity.
#' Used only if \code{score.method = 'contrastReg'}. Default is \code{c(0.5, 2)}.
#'
#' @return coef: Doubly robust estimators of the contrast regression coefficients \code{delta_0}; vector of size \code{p.cate} + 1 (intercept included)
#'         converge: Indicator that the Newton Raphson algorithm converged for \code{delta_0}; boolean
#'
#' @importFrom stats optim


twoarmsurv.dr <- function(ynew, dnew, trt, x.cate, tau0, weightsurv,
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
  while(mae > error.maxNR && iter <= max.iterNR && eta.max < 1000 * tau0) {
    eta <- as.numeric(exp(x.aug %*% beta))
    eta.max <- max(eta)
    if(eta.max == Inf) break
    error <- (trt * (y - eta * f0.predictor / 2 - f1.predictor / 2) * (1 - ps) - (1 - trt) * (y * eta - f0.predictor * eta / 2 - f1.predictor / 2) * ps) / (eta * ps + (1 - ps))
    score <- colSums(x.aug * weightc * error)
    slopewt <- (y + f0.predictor * (trt / ps - 1) / 2 + f1.predictor * ((1 - trt) / (1 - ps) - 1) / 2) * eta * ps * (1 - ps) / (eta * ps + (1 - ps))^2
    slope <- crossprod(x.aug * weightc * slopewt, x.aug)
    if(iter == 0) {
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
  if(converge1 == 0) {
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


#' Estimate the CATE model using specified scoring methods for survival outcomes
#'
#' Coefficients of the CATE estimated with random forest, boosting, naive Poisson, two regression, and contrast regression
#'
#' @param y Observed survival or censoring time; vector of size \code{n}.
#' @param d The event indicator, normally \code{1 = event, 0 = censored}; vector of size \code{n}.
#' @param trt Treatment received; vector of size \code{n} with treatment coded as 0/1.
#' @param x.cate Matrix of \code{p.cate} baseline covariates specified in the outcome model; dimension \code{n} by \code{p.cate}.
#' @param x.ps Matrix of \code{p.ps} baseline covariates specified in the propensity score model; dimension \code{n} by \code{p.ps}.
#' @param x.ipcw Matrix of \code{p.ipw} baseline covariate specified in inverse probability of censoring weighting; dimension \code{n} by \code{p.ipw}.
#' @param yf Follow-up time, interpreted as the potential censoring time; vector of size \code{n} if the potential censoring time is known.
#' @param tau0 The truncation time for defining restricted mean time lost.
#' @param surv.min Lower truncation limit for probability of being censored (positive and very close to 0).
#' @param score.method A vector of one or multiple methods to estimate the CATE score.
#' Allowed values are:  \code{'randomForest'}, \code{'boosting'}, \code{'poisson'}, \code{'twoReg'},
#' \code{'contrastReg'}. Default specifies all 5 methods.
#' @param ps.method A character vector for the method to estimate the propensity score.
#' Allowed values include one of:
#' \code{'glm'} for logistic regression with main effects only (default), or
#' \code{'lasso'} for a logistic regression with main effects and LASSO penalization on
#' two-way interactions (added to the model if interactions are not specified in \code{ps.model}).
#' Relevant only when \code{ps.model} has more than one variable.
#' @param minPS A numerical value (in [0, 1]) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A number above which estimated propensity scores should be trimmed; scalar
#' @param ipcw.method The censoring model. Allowed values are: \code{'breslow'} (Cox regression with Breslow estimator of the baseline survivor function),
#' \code{'aft (exponential)'}, \code{'aft (weibull)'}, \code{'aft (lognormal)'} or \code{'aft (loglogistic)'}. Default is \code{'breslow'}.
#' @param initial.predictor.method A character vector for the method used to get initial
#' outcome predictions conditional on the covariates in \code{cate.model}
#' in \code{score.method = 'twoReg'} and \code{'contrastReg'}. Allowed values include
#' one of \code{'randomForest'}, \code{'boosting'} and \code{'logistic'} (fastest). Default is \code{'randomForest'}.
#' @param tree.depth A positive integer specifying the depth of individual trees in boosting
#' (usually 2-3). Used only if \code{score.method = 'boosting'} or
#' if \code{score.method = 'twoReg'} or \code{'contrastReg'} and
#' \code{initial.predictor.method = 'boosting'}. Default is \code{3}.
#' @param n.trees.rf A positive integer specifying the maximum number of trees in random forest.
#' Used if \code{score.method = 'ranfomForest'} or
#' if \code{initial.predictor.method = 'randomForest'} with
#' \code{score.method = 'twoReg'} or \code{'contrastReg'}. Default is \code{1000}.
#' @param n.trees.boosting A positive integer specifying the maximum number of trees in boosting
#' (usually 100-1000). Used if \code{score.method = 'boosting'} or
#' if \code{initial.predictor.method = 'boosting'} with
#' \code{score.method = 'twoReg'} or \code{'contrastReg'}. Default is \code{150}.
#' @param B A positive integer specifying the number of time cross-fitting is repeated in
#' \code{score.method = 'twoReg'} and \code{'contrastReg'}. Default is \code{3}.
#' @param Kfold A positive integer specifying the number of folds (parts) used in cross-fitting
#' to partition the data in \code{score.method = 'twoReg'} and \code{'contrastReg'}.
#' Default is \code{5}.
#' @param plot.gbmperf A logical value indicating whether to plot the performance measures in
#' boosting. Used only if \code{score.method = 'boosting'} or if \code{score.method = 'twoReg'}
#' or \code{'contrastReg'} and \code{initial.predictor.method = 'boosting'}. Default is \code{TRUE}.
#' @param error.maxNR A numerical value > 0 indicating the minimum value of the mean absolute
#' error in Newton Raphson algorithm. Used only if \code{score.method = 'contrastReg'}.
#' Default is \code{0.001}.
#' @param max.iterNR A positive integer indicating the maximum number of iterations in the
#' Newton Raphson algorithm. Used only if \code{score.method = 'contrastReg'}.
#' Default is \code{100}.
#' @param tune A vector of 2 numerical values > 0 specifying tuning parameters for the
#' Newton Raphson algorithm. \code{tune[1]} is the step size, \code{tune[2]} specifies a
#' quantity to be added to diagonal of the slope matrix to prevent singularity.
#' Used only if \code{score.method = 'contrastReg'}. Default is \code{c(0.5, 2)}.
#'
#' @return Depending on what score.method is, the outputs is a combination of the following:
#'           result.randomForest: Results of random forest fit, for trt = 0 and trt = 1 separately
#'           result.boosting: Results of boosting fit, for trt = 0 and trt = 1 separately
#'           result.poisson: Naive Poisson estimator (beta1 - beta0); vector of length \code{p.cate} + 1
#'           result.twoReg: Two regression estimator (beta1 - beta0); vector of length \code{p.cate} + 1
#'           result.contrastReg: A list of the contrast regression results with 2 elements:
#'               $delta.contrastReg: Contrast regression DR estimator; vector of length \code{p.cate} + 1
#'               $converge.contrastReg: Indicator that the Newton Raphson algorithm converged for \code{delta_0}; boolean

intxsurv <- function(y, d, trt, x.cate, x.ps, x.ipcw, yf = NULL, tau0, surv.min = 0.025,
                     score.method = c("randomForest", "boosting", "poisson", "twoReg", "contrastReg"),
                     ps.method = "glm", minPS = 0.01, maxPS = 0.99, ipcw.method = "breslow",
                     initial.predictor.method = "randomForest",
                     tree.depth = 3, n.trees.rf = 1000, n.trees.boosting = 150,  B = 3, Kfold = 5, plot.gbmperf = TRUE,
                     error.maxNR = 1e-3, max.iterNR = 100, tune = c(0.5, 2)) {

  result <- vector("list", length(score.method))
  names(result) <- c(paste0("result.", score.method))

  N1 <- sum(trt)
  N0 <- sum(1 - trt)
  N <- N1 + N0
  p.aug <- ncol(x.cate) + 1

  ynew <- pmin(y, tau0)
  dnew <- d
  dnew[ynew == tau0] <- 1

  surv.new <- ipcw.surv(y = y, d = d, x.ipcw = x.ipcw, yf = yf, ipcw.method = ipcw.method, tau0 = tau0, surv.min = surv.min)
  weightc <- dnew / surv.new
  weightc <- weightc / mean(weightc)

  datatotrf <- data.frame(y = ynew, d = dnew, x = x.cate)
  colnames(datatotrf) <- c("y", "d", colnames(x.cate))
  datatot <- data.frame(y = ynew / tau0, x = x.cate)
  colnames(datatot) <- c("y", colnames(x.cate))

  ######## Cross fitting ------------------------------------------------------------------

  index1 <- rep(1:Kfold, floor(N1 / Kfold))
  if(N1 > Kfold * floor(N1 / Kfold)) index1 <- c(index1, 1:(N1 - Kfold * floor(N1 / Kfold)))

  index0 <- rep(1:Kfold, floor(N0 / Kfold))
  if(N0 > Kfold * floor(N0 / Kfold)) index0 <- c(index0, Kfold + 1 - 1:(N0 - Kfold * floor(N0 / Kfold)))

  best.iter <- 0

  if (any(c("twoReg", "contrastReg") %in% score.method)) {

    delta.twoReg.mat <- delta.contrastReg.mat <- matrix(NA, B, p.aug)
    converge <- rep(NA, B)

    for(bb in 1:B) {
      index1cv <- sample(x = index1, size = N1, replace = FALSE)
      index0cv <- sample(x = index0, size = N0, replace = FALSE)
      index <- rep(NA, N)
      index[trt == 1] <- index1cv
      index[trt == 0] <- index0cv

      f1.predictcv <- f0.predictcv <- pscv <- rep(NA, N)
      for (k in 1:Kfold) {
        if (initial.predictor.method == "randomForest"){

          datatot_train <- datatotrf[index != k, ]
          x_ps_train <- x.ps[index != k, , drop = FALSE]
          trt_train=trt[index != k]

          datatot_valid <- datatotrf[index == k, ]
          x_ps_valid <- x.ps[index == k, , drop = FALSE]

          data1rf <- datatot_train[trt_train == 1,]
          fit1 <- rfsrc(Surv(y, d) ~ ., data = data1rf, ntree = n.trees.rf)
          surv1.prd <- predict(object = fit1, newdata = datatot_valid)$survival
          time1 <- fit1$time
          m1 <- length(time1)
          timegap1 <- time1 - c(0, time1[-m1])
          f1.predictcv[index == k] <- colSums(timegap1 * rbind(1, t(surv1.prd[, - m1]))) / tau0

          data0rf <- datatot_train[trt_train == 0,]
          fit0 <- rfsrc(Surv(y, d) ~ ., data = data0rf, ntree = n.trees.rf)
          surv0.prd <- predict(object = fit0, newdata = datatot_valid)$survival
          time0 <- fit0$time
          m0 <- length(time0)
          timegap0 <- time0 - c(0, time0[-m0])
          f0.predictcv[index == k] <- colSums(timegap0 * rbind(1, t(surv0.prd[, - m0]))) / tau0

        } else if (initial.predictor.method == "boosting") {

          datatot_train <- datatot[index != k, ]
          x_ps_train <- x.ps[index != k, ]
          trt_train <- trt[index != k]
          weightnew_train <- weightc[index != k]

          datatot_valid <- datatot[index == k, ]
          x_ps_valid <- x.ps[index == k, ]

          data1 <- datatot_train[trt_train == 1, ]
          data0 <- datatot_train[trt_train == 0, ]
          weight1 <- weightnew_train[trt_train == 1]
          weight0 <- weightnew_train[trt_train == 0]

          fit1.boosting <- gbm(y ~ ., data = data1, distribution = "gaussian", interaction.depth = tree.depth, n.trees = n.trees.boosting,
                               cv.folds = 5, weights = weight1)
          best1.iter <- max(gbm.perf(object = fit1.boosting, plot.it = plot.gbmperf, method = "cv"), 10)
          f.predict.ini <- predict.gbm(object = fit1.boosting, newdata = data1, n.trees = best1.iter)
          f.predict.ini <- jitter(f.predict.ini, factor = 0.15)
          withCallingHandlers({
            fit1.gam <- mgcv::gam(as.formula(paste0("y ~ s(f.predict.ini) + ", paste0(setdiff(colnames(data1), "y"), collapse = "+"))),
                                  family = "binomial", data = data1, weight = weight1)
          },
          warning = function(w) {
            if (grepl("non-integer", conditionMessage(w)))
              invokeRestart("muffleWarning") # suppress warnings in gam(): "In eval(family$initialize) : non-integer #successes in a binomial glm!"
          })
          f1.predict.ini <- predict(fit1.boosting, newdata = datatot_valid, best1.iter)
          newdata1 <- data.frame(f.predict.ini = f1.predict.ini, data.frame(datatot_valid))
          colnames(newdata1) <- c("f.predict.ini", colnames(data1))
          f1.predictcv[index == k] <- predict(fit1.gam, newdata = newdata1, type = "response")

          fit0.boosting <- gbm(y ~ ., data = data0, distribution = "gaussian", interaction.depth = tree.depth, n.trees = n.trees.boosting,
                               cv.folds = 5, weights = weight0)
          best0.iter <- max(gbm.perf(object = fit0.boosting, plot.it = plot.gbmperf, method = "cv"), 10)
          f.predict.ini <- predict.gbm(object = fit0.boosting, newdata = data0, n.trees = best0.iter)
          f.predict.ini <- jitter(f.predict.ini, factor = 0.15)
          withCallingHandlers({
            fit0.gam <- mgcv::gam(as.formula(paste0("y ~ s(f.predict.ini) + ", paste0(setdiff(colnames(data0), "y"), collapse = "+"))),
                                  family = "binomial", data = data0, weight = weight0)
          },
          warning = function(w) {
            if (grepl("non-integer", conditionMessage(w)))
              invokeRestart("muffleWarning") # suppress warnings in gam(): "In eval(family$initialize) : non-integer #successes in a binomial glm!"
          })
          f0.predict.ini <- predict.gbm(fit0.boosting, newdata = datatot_valid, best0.iter)
          newdata0 <- data.frame(f.predict.ini = f0.predict.ini, data.frame(datatot_valid))
          colnames(newdata0) <- c("f.predict.ini", colnames(data0))
          f0.predictcv[index == k] <- predict(fit0.gam, newdata = newdata0, type = "response")

          best.iter <- max(best.iter, best1.iter, best0.iter)

        } else if (initial.predictor.method == "logistic") {

          datatot_train <- datatot[index != k, ]
          x_ps_train <- x.ps[index != k, ]
          trt_train <- trt[index != k]
          weightnew_train <- weightc[index != k]

          datatot_valid <- datatot[index == k, ]
          x_ps_valid <- x.ps[index == k, ]

          data1 <- datatot_train[trt_train == 1, ]
          data0 <- datatot_train[trt_train == 0, ]
          weight1 <- weightnew_train[trt_train == 1]
          weight0 <- weightnew_train[trt_train == 0]

          withCallingHandlers({
            fit1.bin <- glm(y ~ ., data = data1, family = "binomial", weights = weight1)
            f1.predictcv[index == k] <- predict(fit1.bin, newdata = datatot_valid, type = "response")
          },
          warning = function(w) { # don't change the = to <- in withCallingHandlers
            if (grepl("non-integer", conditionMessage(w)))
              invokeRestart("muffleWarning") # suppress warnings in glm(): "In eval(family$initialize) : non-integer #successes in a binomial glm!"
          })

          withCallingHandlers({
            fit0.bin <- glm(y ~ ., data = data0, family = "binomial", weights = weight0)
            f0.predictcv[index == k] <- predict(fit0.bin, newdata = datatot_valid, type = "response")
          },
          warning = function(w) { # don't change the = to <- in withCallingHandlers
            if (grepl("non-integer", conditionMessage(w)))
              invokeRestart("muffleWarning")
          })
        } # end of initial.predictor.method

        if (ps.method == "glm") {
          pscv[index == k] <- glm.simplereg.ps(x.ps = x_ps_train, trt = trt_train, xnew = x_ps_valid, minPS = minPS, maxPS = maxPS)
        } else {
          pscv[index == k] <- glm.ps(x.ps = x_ps_train, trt = trt_train, xnew = x_ps_valid, minPS = minPS, maxPS = maxPS)
        }
      } #end of Kfold loops

      f1.predictcv <- pmin(f1.predictcv, 0.99) * tau0
      f0.predictcv <- pmin(f0.predictcv, 0.99) * tau0

      if ("twoReg" %in% score.method) {
        ## bb-th cross fitting two regression estimator
        beta1.final <- onearmsurv.dr(ynew = ynew, dnew = dnew, trt = trt, x.cate = x.cate, tau0 = tau0, weightsurv = weightc, f.predictor = f1.predictcv, ps = pscv)
        beta0.final <- onearmsurv.dr(ynew = ynew, dnew = dnew, trt = 1 - trt, x.cate = x.cate, tau0 = tau0, weightsurv = weightc, f.predictor = f0.predictcv, ps = 1 - pscv)
        delta.twoReg.mat[bb, ] <- as.vector(beta1.final - beta0.final)
      } #end of if ("twoReg" %in% score.method)

      if ("contrastReg" %in% score.method) {
        ## bb-th cross fitting contrast regression estimator
        fit_contrast <- twoarmsurv.dr(ynew = ynew, dnew = dnew, trt = trt, x.cate = x.cate, tau0 = tau0, weightsurv = weightc, ps = pscv, f1.predictor = f1.predictcv, f0.predictor = f0.predictcv, error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune)
        delta.contrastReg.mat[bb, ] <- fit_contrast$coef
        converge[bb] <- fit_contrast$converge
      } #end of if ("contrastReg" %in% score.method)
    } #end of B loop
  } #end of if c("twoReg", "contrastReg") %in% score.method

  if ("randomForest" %in% score.method) {
    ## Random Forest method based on entire data
    data1rf <- datatotrf[trt == 1, ]
    data0rf <- datatotrf[trt == 0, ]

    fit1.rf <- rfsrc(Surv(y, d) ~ ., data = data1rf, ntree = n.trees.rf)
    fit0.rf <- rfsrc(Surv(y, d) ~ ., data = data0rf, ntree = n.trees.rf)
    result$result.randomForest <- list(fit1.rf = fit1.rf, fit0.rf = fit0.rf)
  }

  if ("boosting" %in% score.method) {
    ## Boosting method based on the entire data
    data1 <- datatot[trt == 1, ]
    weight1 <- weightc[trt == 1]
    data0 <- datatot[trt == 0, ]
    weight0 <- weightc[trt == 0]

    fit1.boosting <- gbm(y ~ ., data = data1, distribution = "gaussian", interaction.depth = tree.depth,
                         n.trees = n.trees.boosting, cv.folds = 5, weights = weight1)
    best1.iter <- max(gbm.perf(object = fit1.boosting, plot.it = plot.gbmperf, method = "cv"), 10)
    f.predict.ini <- predict(object = fit1.boosting, newdata = data1, n.trees = best1.iter)
    f.predict.ini <- jitter(f.predict.ini, factor = 0.1)
    withCallingHandlers({
      fit1.gam <- mgcv::gam(as.formula(paste0("y ~ s(f.predict.ini) + ", paste0(setdiff(colnames(data1), "y"), collapse = "+"))),
                            family = "binomial", data = data1, weight = weight1)
    },
    warning = function(w) {
      if (grepl("non-integer", conditionMessage(w)))
        invokeRestart("muffleWarning") # suppress warnings in gam(): "In eval(family$initialize) : non-integer #successes in a binomial glm!"
    })

    fit0.boosting <- gbm(y ~ ., data = data0, distribution = "gaussian", interaction.depth = tree.depth,
                         n.trees = n.trees.boosting, cv.folds = 5, weights = weight0)
    best0.iter <- max(gbm.perf(object = fit0.boosting, plot.it = plot.gbmperf, method = "cv"), 10)
    f.predict.ini <- predict(object = fit0.boosting, newdata = data0, n.trees = best0.iter)
    f.predict.ini <- jitter(f.predict.ini, factor = 0.1)
    withCallingHandlers({
      fit0.gam <- mgcv::gam(as.formula(paste0("y ~ s(f.predict.ini) + ", paste0(setdiff(colnames(data0), "y"), collapse = "+"))),
                            family = "binomial", data = data0, weight = weight0)
    },
    warning = function(w) {
      if (grepl("non-integer #successes", conditionMessage(w)))
        invokeRestart("muffleWarning") # suppress warnings in gam(): "In eval(family$initialize) : non-integer #successes in a binomial glm!"
    })

    result$result.boosting <- list(fit0.boosting = fit0.boosting, fit0.gam = fit0.gam, best0.iter = best0.iter,
                                   fit1.boosting = fit1.boosting, fit1.gam = fit1.gam, best1.iter = best1.iter)
    best.iter <- max(best.iter, best1.iter, best0.iter)
  }

  if ("poisson" %in% score.method) {
    ## Naive Poisson regression method  (score 2)
    withCallingHandlers({
      beta1.ini <- glm((tau0 - ynew) ~ x.cate, weights = weightc, family = "poisson", subset = (trt == 1))$coef
      beta0.ini <- glm((tau0 - ynew) ~ x.cate, weights = weightc, family = "poisson", subset = (trt == 0))$coef
    },
    warning = function(w) {
      if (grepl("non-integer", conditionMessage(w)))
        invokeRestart("muffleWarning") # suppress warnings in glm(): "In dpois(y, mu, log = TRUE) : non-integer x = 0.557886."
    })
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
    if(converge.contrastReg == TRUE) {
      delta.contrastReg <- colMeans(delta.contrastReg.mat[converge == TRUE, , drop = FALSE])
    } else {
      delta.contrastReg <- colMeans(delta.contrastReg.mat)
    }
    names(delta.contrastReg) <- c("(Intercept)", colnames(x.cate))
    result$result.contrastReg <- list(delta.contrastReg = delta.contrastReg, converge.contrastReg = converge.contrastReg)
  }

  result$best.iter <- best.iter

  return(result)
}



#' Calculate the log CATE score given the baseline covariates and follow-up time for specified scoring method methods for survival outcomes
#'
#' Based on intxsurv results of the CATE coefficients estimated with random forest, boosting, naive Poisson, two regression, contrast regression
#'
#' @param fit List of objects generated from intxsurv: outputs of random forest, boosting, naive Poisson, two regression, contrast regression
#' @param x.cate Matrix of \code{p.cate} baseline covariates specified in the outcome model; dimension \code{n} by \code{p.cate}.
#' @param tau0 The truncation time for defining restricted mean time lost.
#' @param score.method A vector of one or multiple methods to estimate the CATE score.
#' Allowed values are: \code{'randomForest'}, \code{'boosting'}, \code{'poisson'}, \code{'twoReg'},
#' \code{'contrastReg'}. Default specifies all 5 methods.
#'
#' @return score.randomForest: Estimated log CATE score for all \code{n} observations with the random forest method; vector of size \code{n}
#'         score.boosting: Estimated log CATE score for all \code{n} observations with the boosting method; vector of size \code{n}
#'         score.poisson: Estimated log CATE score for all \code{n} observations with the naive Poisson method; vector of size \code{n}
#'         score.twoReg: Estimated log CATE score for all \code{n} observations with the two regression method; vector of size \code{n}
#'         score.contrastReg: Estimated log CATE score for all \code{n} observations with the contrast regression method; vector of size \code{n}
#'         score = NA if the corresponding method is not called
#'

scoresurv <- function(fit, x.cate, tau0,
                      score.method = c("randomForest", "boosting", "poisson", "twoReg", "contrastReg")) {

  result <- vector("list", length(score.method))
  names(result) <- paste0("score.", score.method)

  x.aug <- cbind(1, x.cate)
  datanew <- data.frame(x.cate)
  colnames(datanew) <- colnames(x.cate)

  if ("randomForest" %in% score.method) {
    fit0.rf <- fit$result.randomForest$fit0.rf
    fit1.rf <- fit$result.randomForest$fit1.rf

    predict0 <- predict(object = fit0.rf, newdata = datanew)$survival
    time0 <- fit0.rf$time
    m0 <- length(time0)
    timegap0 <- time0 - c(0, time0[-m0])
    rmst0 <- colSums(timegap0 * rbind(1, t(predict0[,-m0])))

    predict1 <- predict(object = fit1.rf, newdata = datanew)$survival
    time1 <- fit1.rf$time
    m1 <- length(time1)
    timegap1 <- time1 - c(0, time1[-m1])
    rmst1 <- colSums(timegap1 * rbind(1, t(predict1[,-m1])))

    result$score.randomForest <- log(tau0 - rmst1) - log(tau0 - rmst0)
  }

  if ("boosting" %in% score.method) {
    fit0.boosting <- fit$result.boosting$fit0.boosting
    fit0.gam <- fit$result.boosting$fit0.gam
    best0.iter <- fit$result.boosting$best0.iter
    fit1.boosting <- fit$result.boosting$fit1.boosting
    fit1.gam <- fit$result.boosting$fit1.gam
    best1.iter <- fit$result.boosting$best1.iter

    predict0 <- predict(object = fit0.boosting, newdata = datanew, n.trees = best0.iter)
    datanew0 <- data.frame(f.predict.ini = predict0, datanew)
    colnames(datanew0) <- c("f.predict.ini", colnames(datanew))
    rmst0 <- predict(object = fit0.gam, newdata = datanew0, type = "response")

    predict1 <- predict(object = fit1.boosting, newdata = datanew, n.trees = best1.iter)
    datanew1 <- data.frame(f.predict.ini = predict1, datanew)
    colnames(datanew1) <- c("f.predict.ini", colnames(datanew))
    rmst1 <- predict(object = fit1.gam, newdata = datanew1, type = "response")

    result$score.boosting <- as.numeric(log(1 - rmst1) - log(1 - rmst0))
  }

  if ("poisson" %in% score.method) {
    delta.poisson <- fit$result.poisson
    result$score.poisson <- as.numeric(x.aug %*% delta.poisson)
  }

  if ("twoReg" %in% score.method) {
    delta.twoReg <- fit$result.twoReg
    result$score.twoReg <- as.numeric(x.aug %*% delta.twoReg)
  }

  if ("contrastReg" %in% score.method) {
    delta.contrastReg <- fit$result.contrastReg$delta.contrastReg
    result$score.contrastReg <- as.numeric(x.aug %*% delta.contrastReg)
  }

  return(result)
}
