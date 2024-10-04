#' Estimation of the conditional average treatment effect (CATE) score for survival data
#'
#' Provides singly robust and doubly robust estimation of CATE score for survival data with
#' up to 5 scoring methods among the following: Random forest, boosting, poisson regression,
#' two regressions, and contrast regression.
#'
#' @param cate.model A standard \code{Surv} formula describing the outcome model to be fitted.
#' The outcome must appear on the left-hand side.
#' @param ps.model A formula describing the propensity score (PS) model to be fitted. The treatment must
#' appear on the left-hand side. The treatment must be a numeric vector coded as 0/1.
#' If data are from a randomized controlled trial, specify \code{ps.model = ~1} as an intercept-only model.
#' @param data A data frame containing the variables in the outcome, propensity score, and inverse
#' probability of censoring models (if specified); a data frame with \code{n} rows (1 row per observation).
#' @param score.method A vector of one or multiple methods to estimate the CATE score.
#' Allowed values are: \code{'randomForest'}, \code{'boosting'}, \code{'poisson'}, \code{'twoReg'}, and
#' \code{'contrastReg'}.
#' @param ipcw.model A formula describing the inverse probability of censoring weighting (IPCW)
#' model to be fitted. The left-hand side must be empty. Default is \code{ipcw.model = NULL},
#' which corresponds to specifying the IPCW model with the same covariates as the outcome model
#' \code{cate.model} plus the treatment.
#' @param followup.time A column name in \code{data} specifying the maximum follow-up time,
#' interpreted as the potential censoring time. Default is \code{followup.time = NULL},
#' which corresponds to unknown potential censoring time.
#' @param tau0 The truncation time for defining restricted mean time lost. Default is \code{NULL},
#' which corresponds to setting the truncation time as the maximum survival time in the data.
#' @param surv.min Lower truncation limit for the probability of being censored.
#' It must be a positive value and should be chosen close to 0. Default is \code{0.025}.
#' @param ipcw.method A character value for the censoring model. Allowed values are:
#' \code{'breslow'} (Cox regression with Breslow estimator of the baseline survivor function),
#' \code{'aft (exponential)'}, \code{'aft (weibull)'}, \code{'aft (lognormal)'}
#' or \code{'aft (loglogistic)'} (accelerated failure time model with different distributions for
#' y variable). Default is \code{'breslow'}.
#' @param higher.y A logical value indicating whether higher (\code{TRUE}) or
#' lower (\code{FALSE}) values of the outcome are more desirable. Default is \code{TRUE}.
#' @param prop.cutoff A vector of numerical values (in `(0, 1]`) specifying percentiles of the
#' estimated log CATE scores to define nested subgroups. Each element represents the cutoff to
#' separate observations in nested subgroups (below vs above cutoff).
#' The length of \code{prop.cutoff} is the number of nested subgroups.
#' An equally-spaced sequence of proportions ending with 1 is recommended.
#' Default is \code{seq(0.5, 1, length = 6)}.
#' @param ps.method A character value for the method to estimate the propensity score.
#' Allowed values include one of:
#' \code{'glm'} for logistic regression with main effects only (default), or
#' \code{'lasso'} for a logistic regression with main effects and LASSO penalization on
#' two-way interactions (added to the model if interactions are not specified in \code{ps.model}).
#' Relevant only when \code{ps.model} has more than one variable.
#' @param minPS A numerical value (in `[0, 1]`) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in `(0, 1]`) above which estimated propensity scores should be
#' truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.
#' @param initial.predictor.method A character vector for the method used to get initial
#' outcome predictions conditional on the covariates specified in \code{cate.model}. Only applies
#' when \code{score.method} includes \code{'twoReg'} or \code{'contrastReg'}. Allowed values include
#' one of \code{'randomForest'}, \code{'boosting'} and \code{'logistic'} (fastest).
#' Default is \code{'randomForest'}.
#' @param tree.depth A positive integer specifying the depth of individual trees in boosting
#' (usually 2-3). Used only if \code{score.method = 'boosting'} or
#' if \code{initial.predictor.method = 'boosting'} with \code{score.method = 'twoReg'} or
#' \code{'contrastReg'}. Default is 2.
#' @param n.trees.rf A positive integer specifying the maximum number of trees in random forest.
#' Used if \code{score.method = 'ranfomForest'} or if \code{initial.predictor.method = 'randomForest'}
#' with \code{score.method = 'twoReg'} or \code{'contrastReg'}. Only applies for survival outcomes.
#' Default is \code{1000}.
#' @param n.trees.boosting A positive integer specifying the maximum number of trees in boosting
#' (usually 100-1000). Used if \code{score.method = 'boosting'} or
#' if \code{initial.predictor.method = 'boosting'} with \code{score.method = 'twoReg'} or
#' \code{'contrastReg'}. Default is \code{200}.
#' @param B A positive integer specifying the number of time cross-fitting is repeated in
#' \code{score.method = 'twoReg'} and \code{'contrastReg'}. Default is \code{3}.
#' @param Kfold A positive integer specifying the number of folds used in cross-fitting
#' to partition the data in \code{score.method = 'twoReg'} and \code{'contrastReg'}.
#' Default is \code{5}.
#' @param error.maxNR A numerical value > 0 indicating the minimum value of the mean absolute
#' error in Newton Raphson algorithm. Used only if \code{score.method = 'contrastReg'}.
#' Default is \code{0.001}.
#' @param max.iterNR A positive integer indicating the maximum number of iterations in the
#' Newton Raphson algorithm. Used only if \code{score.method = 'contrastReg'}.
#' Default is \code{150}.
#' @param tune A vector of 2 numerical values > 0 specifying tuning parameters for the
#' Newton Raphson algorithm. \code{tune[1]} is the step size, \code{tune[2]} specifies a
#' quantity to be added to diagonal of the slope matrix to prevent singularity.
#' Used only if \code{score.method = 'contrastReg'}. Default is \code{c(0.5, 2)}.
#' @param seed An optional integer specifying an initial randomization seed for reproducibility.
#' Default is \code{NULL}, corresponding to no seed.
#' @param verbose An integer value indicating what kind of intermediate progress messages should
#' be printed. \code{0} means no outputs. \code{1} means only progress and run time.
#' \code{2} means progress, run time, and all errors and warnings. Default is \code{0}.
#' @param plot.gbmperf A logical value indicating whether to plot the performance measures in
#' boosting. Used only if \code{score.method = 'boosting'} or if \code{score.method = 'twoReg'}
#' or \code{'contrastReg'} and \code{initial.predictor.method = 'boosting'}. Default is \code{TRUE}.
#' @param ... Additional arguments for \code{gbm()}
#'
#' @return Returns an object of the class \code{catefit} containing the following components:
#' \itemize{
#'  \item{\code{ate.randomForest}: } A vector of numerical values of length \code{prop.cutoff}
#'  containing the estimated ATE by the RMTL ratio in nested subgroups (defined by \code{prop.cutoff})
#'  constructed based on the estimated CATE scores with random forest method.
#'  Only provided if \code{score.method} includes \code{'randomForest'}.
#'  \item{\code{ate.boosting}: }Same as \code{ate.randomForest}, but with the nested subgroups based
#'  the estimated CATE scores with boosting. Only provided if \code{score.method}
#'  includes \code{'boosting'}.
#'  \item{\code{ate.poisson}: }Same as \code{ate.randomForest}, but with the nested subgroups based
#'  the estimated CATE scores with poisson regression.
#'  Only provided if \code{score.method} includes \code{'poisson'}.
#'  \item{\code{ate.twoReg}: }Same as \code{ate.randomForest}, but with the nested subgroups based
#'  the estimated CATE scores with two regressions.
#'  Only provided if \code{score.method} includes \code{'twoReg'}.
#'  \item{\code{ate.contrastReg}: }Same as \code{ate.randomForest}, but with the nested subgroups based
#'  the estimated CATE scores with contrast regression.
#'  Only provided if \code{score.method} includes \code{'contrastReg'}.
#'  \item{\code{hr.randomForest}: }A vector of numerical values of length \code{prop.cutoff}
#'  containing the adjusted hazard ratio in nested subgroups (defined by \code{prop.cutoff})
#'  constructed based on the estimated CATE scores with random forest method.
#'  Only provided if \code{score.method} includes \code{'randomForest'}.
#'  \item{\code{hr.boosting}: }Same as \code{hr.randomForest}, but with the nested subgroups based
#'  the estimated CATE scores with boosting. Only provided if \code{score.method}
#'  includes \code{'boosting'}.
#'  \item{\code{hr.poisson}: }Same as \code{hr.randomForest}, but with the nested subgroups based
#'  the estimated CATE scores with poisson regression.
#'  Only provided if \code{score.method} includes \code{'poisson'}.
#'  \item{\code{hr.twoReg}: }Same as \code{hr.randomForest}, but with the nested subgroups based
#'  the estimated CATE scores with two regressions.
#'  Only provided if \code{score.method} includes \code{'twoReg'}.
#'  \item{\code{hr.contrastReg}: }Same as \code{hr.randomForest}, but with the nested subgroups based
#'  the estimated CATE scores with contrast regression.
#'  Only provided if \code{score.method} includes \code{'contrastReg'}.
#'  \item{\code{score.randomForest}: }A vector of numerical values of length n
#'  (number of observations in \code{data}) containing the estimated log-CATE scores
#'  according to random forest. Only provided if \code{score.method}
#'  includes \code{'randomForest'}.
#'  \item{\code{score.boosting}: }Same as \code{score.randomForest}, but with estimated log-CATE score
#'  according to boosting. Only provided if \code{score.method} includes
#'  \code{'boosting'}.
#'  \item{\code{score.poisson}: }Same as \code{score.randomForest}, but with estimated log-CATE score
#'  according to the Poisson regression. Only provided if \code{score.method}
#'  includes \code{'poisson'}.
#'  \item{\code{score.twoReg}: }Same as \code{score.randomForest}, but with estimated log-CATE score
#'  according to two regressions. Only provided if \code{score.method} includes
#'  \code{'twoReg'}.
#'  \item{\code{score.contrastReg}: }Same as \code{score.randomForest}, but with estimated log-CATE score
#'  according to contrast regression. Only provided if \code{score.method} includes
#'  \code{'contrastReg'}.
#'  \item{\code{fit}: }Additional details on model fitting if \code{score.method}
#'  includes 'randomForest', 'boosting' or 'contrastReg':
#'  \itemize{
#'    \item{\code{result.randomForest}: }Details on the random forest model fitted to observations
#'    with treatment = 0 \code{($fit0.rf)} and to observations with treatment = 1 \code{($fit1.rf)}.
#'    Only provided if \code{score.method} includes \code{'randomForest'}.
#'    \item{\code{result.boosting}: }Details on the boosting model fitted to observations
#'    with treatment = 0, \code{($fit0.boosting)} and \code{($fit0.gam)} and to observations with treatment = 1,
#'    \code{($fit1.boosting)} and \code{($fit1.gam)}.
#'    Only provided if \code{score.method} includes \code{'boosting'}.
#'    \item{\code{result.contrastReg$converge.contrastReg}: }Whether the contrast regression algorithm converged
#'    or not. Only provided if \code{score.method} includes \code{'contrastReg'}.
#'  }
#'  \item{\code{coefficients}: }A data frame with the coefficients of the estimated log-CATE
#'  score by \code{score.method}. The data frame has number of rows equal to the number of
#'  covariates in \code{cate.model} and number of columns equal to \code{length(score.method)}.
#'  If \code{score.method} includes \code{'contrastReg'}, the data frame has an additional
#'  column containing the standard errors of the coefficients estimated with contrast regression.
#'  \code{'randomForest'} and \code{'boosting'} do not have coefficient results because
#'  tree-based methods typically do not express the log-CATE as a linear combination of coefficients
#'  and covariates.
#'  \item{\code{errors/warnings}: }A nested list of errors and warnings that were wrapped during the
#'  calculation of ATE. Errors and warnings are organized by \code{score.method}.
#' }
#'
#' @details The CATE score represents an individual-level treatment effect for survival data,
#' estimated with random forest, boosting, Poisson regression, and the doubly
#' robust estimator (two regressions, Yadlowsky, 2020) applied separately by treatment group
#' or with the other doubly robust estimators (contrast regression, Yadlowsky, 2020) applied
#' to the entire data set.
#'
#' \code{\link{catefitsurv}()} provides the coefficients of the CATE score for each scoring method requested
#' through \code{score.method}. Currently, contrast regression is the only method which allows
#' for inference of the CATE coefficients by providing standard errors of the coefficients.
#' The coefficients can be used to learn the effect size of each variable and predict the
#' CATE score for a new observation.
#'
#' \code{\link{catefitsurv}()} also provides the predicted CATE score of each observation in the data set,
#' for each scoring method. The predictions allow ranking the observations from potentially
#' high responders to the treatment to potentially low or standard responders.
#'
#' The estimated ATE among nested subgroups of high responders are also provided by scoring method.
#' Note that the ATEs in \code{\link{catefitsurv}()} are derived based on the CATE score which is estimated
#' using the same data sample. Therefore, overfitting may be an issue. \code{\link{catecvsurv}()} is more
#' suitable to inspect the estimated ATEs across scoring methods as it implements internal cross
#' validation to reduce optimism.
#'
#' @references Yadlowsky, S., Pellegrini, F., Lionetto, F., Braune, S., & Tian, L. (2020).
#' \emph{Estimation and validation of ratio-based conditional average treatment effects using
#' observational data. Journal of the American Statistical Association, 1-18.} DOI: 10.1080/01621459.2020.1772080.
#'
#' @seealso \code{\link{catecvsurv}()}
#'
#' @examples
#' \donttest{
#' library(survival)
#'
#' tau0 <- with(survivalExample, min(quantile(y[trt == "drug1"], 0.95),
#'                                quantile(y[trt == "drug0"], 0.95)))
#'
#' fit <- catefitsurv(data = survivalExample,
#'                    score.method = "randomForest",
#'                    cate.model = Surv(y, d) ~ age + female + previous_cost +
#'                                              previous_number_relapses,
#'                    ps.model = trt ~ age + previous_treatment,
#'                    ipcw.model = ~ age + previous_cost + previous_treatment,
#'                    tau0 = tau0,
#'                    seed = 999)
#'
#' coef(fit)
#' }
#'
#' @export

catefitsurv <- function(data,
                        score.method,
                        cate.model,
                        ps.model,
                        ps.method = "glm",
                        initial.predictor.method = "randomForest",
                        ipcw.model = NULL,
                        ipcw.method = "breslow",
                        minPS = 0.01,
                        maxPS = 0.99,
                        followup.time = NULL,
                        tau0 = NULL,
                        higher.y = TRUE,
                        prop.cutoff = seq(0.5, 1, length = 6),
                        surv.min = 0.025,
                        tree.depth = 2,
                        n.trees.rf = 1000,
                        n.trees.boosting = 200,
                        B = 3,
                        Kfold = 5,
                        plot.gbmperf = TRUE,
                        error.maxNR = 1e-3,
                        max.iterNR = 100,
                        tune = c(0.5, 2),
                        seed = NULL,
                        verbose = 0,
                        ...) {

  # Set seed for reproducibility
  set.seed(seed)

  if (verbose >= 1) t.start <- Sys.time()

  # CHECK ARGUMENTS
  arg.checks(
    fun = "catefit", response = "survival", data = data,
    followup.time = followup.time, tau0 = tau0, surv.min = surv.min,
    higher.y = higher.y, score.method = score.method, abc = FALSE,
    prop.cutoff = prop.cutoff,
    ps.method = ps.method, minPS = minPS, maxPS = maxPS,
    ipcw.method = ipcw.method,
    train.prop = 0.5, cv.n = 1,
    error.max = 1, max.iter = 2,
    initial.predictor.method = initial.predictor.method,
    tree.depth = tree.depth, n.trees.rf = n.trees.rf,
    n.trees.boosting = n.trees.boosting,
    B = B, Kfold = Kfold, plot.gbmperf = plot.gbmperf,
    error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune
  )

  # PRE-PROCESSING
  out <- data.preproc.surv(fun = "catefit", cate.model = cate.model,
                           ps.model = ps.model, ipcw.model = ipcw.model,
                           tau0 = tau0, data = data, prop.cutoff = prop.cutoff,
                           ps.method = ps.method, response = "survival")
  y <- out$y
  d <- out$d
  trt <- out$trt
  x.ps <- out$x.ps
  x.cate <- out$x.cate
  x.ipcw <- out$x.ipcw
  prop <- out$prop
  prop.no1 <- out$prop.no1
  yf <- NULL

  if (!is.null(followup.time)) {
    yf <- data[[followup.time]]
  }

  # Check if tau0 is large enough, i.e., if tau0 is larger than the 50% quantile of the observed or censoring time.
  if (tau0 < median(y)) warning("It is recommended to increase tau0 close to the largest observed or censored time.")

  # FUNCTION STARTS HERE
  result <- vector("list", 3 * length(score.method) + 2)
  names(result) <- c(paste0("ate.", score.method), paste0("hr.", score.method), paste0("score.", score.method), "fit", "coefficients")

  # Add names to errors/warnings
  result[["errors/warnings"]] <- vector("list", length(score.method))
  names(result[["errors/warnings"]]) <- score.method
  result[["errors/warnings"]] <- lapply(result[["errors/warnings"]],
                                        function(x) {
                                          x <- vector("list", 1)
                                          names(x) <- c("est.high")
                                          return(x)
                                        })

  # Fit the interaction model
  fit <- intxsurv(y = y, d = d, trt = trt, x.cate = x.cate, x.ps = x.ps, x.ipcw = x.ipcw, yf = yf, tau0 = tau0, surv.min = surv.min,
                  score.method = score.method,
                  ps.method = ps.method, minPS = minPS, maxPS = maxPS, ipcw.method = ipcw.method,
                  initial.predictor.method = initial.predictor.method,
                  tree.depth = tree.depth, n.trees.rf = n.trees.rf, n.trees.boosting = n.trees.boosting,
                  B = B, Kfold = Kfold, plot.gbmperf = plot.gbmperf,
                  error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune, ...)

  if (initial.predictor.method == "boosting" & fit$best.iter == n.trees.boosting) {
    warning(paste("The best boosting iteration was iteration number", n.trees.boosting, " out of ", n.trees.boosting, ". Consider increasing the maximum number of trees and turning on boosting performance plot (plot.gbmperf = TRUE).", sep = ""))
  }

  # Construct the score in the whole dataset (no training or validation for catefitsurv)
  fit.score <- scoresurv(fit = fit, x.cate = x.cate, tau0 = tau0, score.method = score.method)

  # Estimate the treatment effect in the whole dataset
  errors <- warnings <- c()
  est.prop1 <- NULL

  if (prop[length(prop)] == 1) {
    est.prop1 <- survCatch(drsurv(y = y, d = d, x.cate = x.cate, x.ps = x.ps, x.ipcw = x.ipcw, trt = trt, yf = yf,
                                  tau0 = tau0, surv.min = surv.min, ps.method = ps.method, minPS = minPS, maxPS = maxPS,
                                  ipcw.method = ipcw.method))
    est.prop1$ate.rmtl.high <- exp(est.prop1$log.rmtl.ratio)
    est.prop1$hr.high <- exp(est.prop1$log.hazard.ratio)
  }

  for (name in names(fit.score)) {
    score <- fit.score[[name]]
    name.method <- str_extract(name, "(?<=^score\\.).*$")
    result[[name]] <- score
    est.onlyhigh <- estsurv.bilevel.subgroups(y = y, d = d, x.cate = x.cate, x.ps = x.ps, x.ipcw = x.ipcw, trt = trt, yf = yf,
                                              tau0 = tau0, surv.min = surv.min,
                                              score = score, higher.y = higher.y,
                                              prop = prop.no1, onlyhigh = TRUE,
                                              ps.method = ps.method, minPS = minPS, maxPS = maxPS, ipcw.method = ipcw.method)
    result[[str_replace(name, "score", "ate")]] <- c(est.onlyhigh$ate.rmtl.high, est.prop1$ate.rmtl.high)
    result[[str_replace(name, "score", "hr")]] <- c(est.onlyhigh$hr.high, est.prop1$hr.high)

    names(result[[str_replace(name, "score", "ate")]]) <-
      names(result[[str_replace(name, "score", "hr")]]) <- paste0("prop", round(prop, 2))

    est.onlyhigh$err.high$`prop 1.0` <- est.prop1$errors
    err.onlyhigh <- est.onlyhigh$err.high[!sapply(est.onlyhigh$err.high, is.null)]
    if (length(err.onlyhigh) == 0) err.onlyhigh <- NULL

    est.onlyhigh$warn.high$`prop 1.0` <- est.prop1$warnings
    warn.onlyhigh <- est.onlyhigh$warn.high[!sapply(est.onlyhigh$warn.high, is.null)]
    if (length(warn.onlyhigh) == 0) warn.onlyhigh <- NULL

    result[['errors/warnings']][[name.method]]$est.high$errors <- err.onlyhigh
    result[['errors/warnings']][[name.method]]$est.high$warnings <- warn.onlyhigh

    if (!is.null(err.onlyhigh)) errors <- c(errors, name.method)
    if (!is.null(warn.onlyhigh)) warnings <- c(warnings, name.method)
  }# end of for (name in names(fit.score)) {}

  if (length(errors) != 0) {
    cat(paste0('    Warning: Error(s) occurred when estimating the ATEs in the nested subgroup using "', paste0(errors, collapse = '", "'), '";\n    return NAs in the corresponding subgroup.'),'\n')
    warning(paste0('Warning: Error(s) occurred when estimating the ATEs in the nested subgroup using "', paste0(errors, collapse = '", "'), "\". NAs are returned for RMTL ratio and HR in the corresponding subgroup; see 'errors/warnings'."))
  }

  if (length(warnings) != 0) {
    cat(paste0('    Warning: Warning(s) occurred when estimating the ATEs in the nested subgroup using "', paste0(warnings, collapse = '", "'), '"'),'\n')
    warning(paste0('Warning: Warning(s) occurred when estimating the ATEs in the nested subgroup using "', paste0(warnings, collapse = '", "'), "\"; see 'errors/warnings'."))
  }

  if (sum(score.method %in% c("poisson", "twoReg", "contrastReg")) > 0) {
    cf <- data.frame(matrix(NA, nrow = ncol(x.cate) + 1, ncol = 3))
    colnames(cf) <- c("poisson", "twoReg", "contrastReg")
    rownames(cf) <- c("(Intercept)", colnames(x.cate))

    if ("poisson" %in% score.method) cf$poisson <- fit$result.poisson

    if ("twoReg" %in% score.method) cf$twoReg <- fit$result.twoReg

    if ("contrastReg" %in% score.method) {
      cf$contrastReg <- fit$result.contrastReg$delta.contrastReg
    }

    result$coefficients <- cf[, colSums(is.na(cf)) != nrow(cf), drop = FALSE]
  }

  if ("randomForest" %in% score.method) result$fit$result.randomForest <-
    fit$result.randomForest
  if ("boosting" %in% score.method) result$fit$result.boosting <-
    fit$result.boosting
  if ("contrastReg" %in% score.method) result$fit$result.contrastReg$converge.contrastReg <-
    fit$result.contrastReg$converge.contrastReg

  if (verbose >= 1) {
    t.end <- Sys.time()
    t.diff <- round(difftime(t.end, t.start),2)
    cat('Total runtime :',as.numeric(t.diff), attributes(t.diff)$units, '\n')
  }

  result$response <- "survival"
  result$score.method <- score.method

  class(result) <- "catefit"
  return(result)
}

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
#' @param minPS A numerical value (in `[0, 1]`) below which estimated propensity scores should be
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
#' @param ... Additional arguments for \code{gbm()}
#'
#' @return Depending on what score.method is, the outputs is a combination of the following:
#'           result.randomForest: Results of random forest fit, for trt = 0 and trt = 1 separately
#'           result.boosting: Results of boosting fit, for trt = 0 and trt = 1 separately
#'           result.poisson: Naive Poisson estimator (beta1 - beta0); vector of length \code{p.cate} + 1
#'           result.twoReg: Two regression estimator (beta1 - beta0); vector of length \code{p.cate} + 1
#'           result.contrastReg: A list of the contrast regression results with 2 elements:
#'               $delta.contrastReg: Contrast regression DR estimator; vector of length \code{p.cate} + 1
#'               $converge.contrastReg: Indicator that the Newton Raphson algorithm converged for \code{delta_0}; boolean

intxsurv <- function(y, d, trt, x.cate, x.ps, x.ipcw, yf = NULL, tau0,
                     surv.min = 0.025,
                     score.method = c("randomForest", "boosting", "poisson",
                                      "twoReg", "contrastReg"),
                     ps.method = "glm", minPS = 0.01, maxPS = 0.99,
                     ipcw.method = "breslow",
                     initial.predictor.method = "randomForest",
                     tree.depth = 3, n.trees.rf = 1000, n.trees.boosting = 150,
                     B = 3, Kfold = 5, plot.gbmperf = TRUE,
                     error.maxNR = 1e-3, max.iterNR = 100,
                     tune = c(0.5, 2), ...) {

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

  # Cross fitting

  index1 <- rep(1:Kfold, floor(N1 / Kfold))
  if (N1 > Kfold * floor(N1 / Kfold)) index1 <- c(index1, 1:(N1 - Kfold * floor(N1 / Kfold)))

  index0 <- rep(1:Kfold, floor(N0 / Kfold))
  if (N0 > Kfold * floor(N0 / Kfold)) index0 <- c(index0, Kfold + 1 - 1:(N0 - Kfold * floor(N0 / Kfold)))

  best.iter <- 0

  if (any(c("twoReg", "contrastReg") %in% score.method)) {

    delta.twoReg.mat <- delta.contrastReg.mat <- matrix(NA, B, p.aug)
    converge <- rep(NA, B)

    for (bb in seq(B)) {
      index1cv <- sample(x = index1, size = N1, replace = FALSE)
      index0cv <- sample(x = index0, size = N0, replace = FALSE)
      index <- rep(NA, N)
      index[trt == 1] <- index1cv
      index[trt == 0] <- index0cv

      f1.predictcv <- f0.predictcv <- pscv <- rep(NA, N)
      for (k in seq(Kfold)) {
        if (initial.predictor.method == "randomForest") {

          datatot_train <- datatotrf[index != k, ]
          x_ps_train <- x.ps[index != k, , drop = FALSE]
          trt_train <- trt[index != k]

          datatot_valid <- datatotrf[index == k, ]
          x_ps_valid <- x.ps[index == k, , drop = FALSE]

          data1rf <- datatot_train[trt_train == 1,]
          fit1 <- rfsrc(Surv(y, d) ~ ., data = data1rf, ntree = n.trees.rf)
          surv1.prd <- predict(object = fit1, newdata = datatot_valid)$survival
          time1 <- fit1$time
          m1 <- length(time1)
          timegap1 <- time1 - c(0, time1[-m1])
          f1.predictcv[index == k] <- colSums(timegap1 * rbind(1, t(surv1.prd[, -m1]))) / tau0

          data0rf <- datatot_train[trt_train == 0,]
          fit0 <- rfsrc(Surv(y, d) ~ ., data = data0rf, ntree = n.trees.rf)
          surv0.prd <- predict(object = fit0, newdata = datatot_valid)$survival
          time0 <- fit0$time
          m0 <- length(time0)
          timegap0 <- time0 - c(0, time0[-m0])
          f0.predictcv[index == k] <- colSums(timegap0 * rbind(1, t(surv0.prd[, -m0]))) / tau0

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
                               cv.folds = 5, weights = weight1, ...)
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
                               cv.folds = 5, weights = weight0, ...)
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
                         n.trees = n.trees.boosting, cv.folds = 5, weights = weight1, ...)
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
                         n.trees = n.trees.boosting, cv.folds = 5, weights = weight0, ...)
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
    if (converge.contrastReg == TRUE) {
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
