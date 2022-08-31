# ------------------------------------------------------------------
#
# Project: Precision Medicine MS (precmed) - Comprehensive R package
#
# Purpose: Output functions for Continuous outcomes
#
# Platform: Windows
# R Version: 4.1.0
#



#' Cross-validation of the conditional average treatment effect (CATE) score for continuous outcomes
#'
#' Provides doubly robust estimation of the average treatment effect (ATE) in nested and
#' mutually exclusive subgroups of patients defined by an estimated conditional average
#' treatment effect (CATE) score via cross-validation (CV). The CATE score can be estimated
#' with up to 6 methods among the following: Linear regression, boosting, two regressions,
#' contrast regression, random forest and generalized additive model (see \code{score.method}).
#'
#' @param cate.model A formula describing the outcome model to be fitted.
#' The outcome must appear on the left-hand side.
#' @param init.model A formula describing the initial predictor model. The outcome must appear on the left-hand side.
#' It must be specified when \code{score.method = contrastReg} or \code{twoReg}.
#' @param ps.model A formula describing the propensity score model to be fitted.
#' The treatment must appear on the left-hand side. The treatment must be a numeric vector
#' coded as 0/1. If data are from a RCT, specify \code{ps.model} as an intercept-only model.
#' @param data A data frame containing the variables in the outcome and propensity score models;
#' a data frame with \code{n} rows (1 row per observation).
#' @param score.method A vector of one or multiple methods to estimate the CATE score.
#' Allowed values are: \code{'boosting'}, \code{'gaussian'}, \code{'twoReg'}, \code{'contrastReg'},
#' \code{'randomForest'}, \code{'gam'}.
#' @param higher.y A logical value indicating whether higher (\code{TRUE}) or lower (\code{FALSE})
#' values of the outcome are more desirable. Default is \code{TRUE}.
#' @param abc A logical value indicating whether the area between curves (ABC) should be calculated
#' at each cross-validation iterations, for each \code{score.method}. Default is \code{TRUE}.
#' @param prop.cutoff A vector of numerical values (in (0, 1]) specifying percentiles of the
#' estimated CATE scores to define nested subgroups. Each element represents the cutoff to
#' separate observations in nested subgroups (below vs above cutoff).
#' The length of \code{prop.cutoff} is the number of nested subgroups.
#' An equally-spaced sequence of proportions ending with 1 is recommended.
#' Default is \code{seq(0.5, 1, length = 6)}.
#' @param prop.multi A vector of numerical values (in [0, 1]) specifying percentiles of the
#' estimated CATE scores to define mutually exclusive subgroups.
#' It should start with 0, end with 1, and be of \code{length(prop.multi) > 2}.
#' Each element represents the cutoff to separate the observations into
#' \code{length(prop.multi) - 1} mutually exclusive subgroups.
#' Default is \code{c(0, 1/3, 2/3, 1)}.
#' @param ps.method A character value for the method to estimate the propensity score.
#' Allowed values include one of:
#' \code{'glm'} for logistic regression with main effects only (default), or
#' \code{'lasso'} for a logistic regression with main effects and LASSO penalization on
#' two-way interactions (added to the model if interactions are not specified in \code{ps.model}).
#' Relevant only when \code{ps.model} has more than one variable.
#' @param minPS A numerical value (in [0, 1]) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in (0, 1]) above which estimated propensity scores should be
#' truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.
#' @param train.prop A numerical value (in (0, 1)) indicating the proportion of total data used
#' for training. Default is \code{3/4}.
#' @param cv.n A positive integer value indicating the number of cross-validation iterations.
#' Default is \code{10}.
#' @param error.max A numerical value > 0 indicating the tolerance (maximum value of error)
#' for the largest standardized absolute difference in the covariate distributions or in the
#' doubly robust estimated rate ratios between the training and validation sets. This is used
#' to define a balanced training-validation splitting. Default is \code{0.1}.
#' @param max.iter A positive integer value indicating the maximum number of iterations when
#' searching for a balanced training-validation split. Default is \code{5,000}.
#' @param initial.predictor.method A character vector for the method used to get initial
#' outcome predictions conditional on the covariates in \code{cate.model}
#' in \code{score.method = 'twoReg'} and \code{'contrastReg'}. Allowed values include
#' one of \code{'poisson'} (fastest), \code{'boosting'} and \code{'gam'}.
#' Default is \code{'boosting'}.
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
#' @param n.trees.rf A positive integer specifying the maximum number of trees in random forest.
#' Used if \code{score.method = 'ranfomForest'} or if \code{initial.predictor.method = 'randomForest'}
#' with \code{score.method = 'twoReg'} or \code{'contrastReg'}. Only applies for survival outcomes.
#' Default is \code{1000}.
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
#' @param error.maxNR A numerical value > 0 indicating the minimum value of the mean absolute
#' error in Newton Raphson algorithm. Used only if \code{score.method = 'contrastReg'}.
#' Default is \code{0.001}.
#' @param tune A vector of 2 numerical values > 0 specifying tuning parameters for the
#' Newton Raphson algorithm. \code{tune[1]} is the step size, \code{tune[2]} specifies a
#' quantity to be added to diagonal of the slope matrix to prevent singularity.
#' Used only if \code{score.method = 'contrastReg'}. Default is \code{c(0.5, 2)}.
#' @param seed An optional integer specifying an initial randomization seed for reproducibility.
#' Default is \code{NULL}, corresponding to no seed.
#' @param verbose An integer value indicating what kind of intermediate progress messages should
#' be printed. \code{0} means no outputs. \code{1} means only progress bar and run time.
#' \code{2} means progress bar, run time, and all errors and warnings. Default is \code{2}.
#' @param ... Additional arguments for \code{gbm()}
#'
#' @return Returns a list containing the following components saved as a \code{"precmed"} object:
#' \itemize{
#'  \item{\code{ate.gaussian}: }{A list of results output if \code{score.method} includes
#'  \code{'gaussian'}:}
#'  \itemize{
#'     \item{\code{ate.est.train.high.cv}: }{A matrix of numerical values with
#'     \code{length(prop.cutoff)} rows and \code{cv.n} columns.
#'     The ith column/jth row cell contains the estimated ATE in the nested subgroup of high responders
#'     defined by CATE score above (if \code{higher.y = TRUE}) or below (if \code{higher.y = FALSE}) the
#'     \code{prop.cutoff[j]}x100\% percentile of the estimated CATE score in the training set in the ith
#'     cross-validation iteration.}
#'     \item{\code{ate.est.train.low.cv}: }{A matrix of numerical values with
#'     \code{length(prop.cutoff) - 1} rows and \code{cv.n} columns.
#'     The ith column/jth row cell contains the estimated ATE in the nested subgroup of low responders
#'     defined by CATE score below (if \code{higher.y = TRUE}) or above (if \code{higher.y = FALSE}) the
#'     \code{prop.cutoff[j]}x100\% percentile of the estimated CATE score in the training set in the ith
#'     cross-validation iteration.}
#'     \item{\code{ate.est.valid.high.cv}: }{Same as \code{ate.est.train.high.cv},
#'     but in the validation set.}
#'     \item{\code{ate.est.valid.low.cv}: }{Same as \code{ate.est.train.low.cv},
#'     but in the validation set.}
#'     \item{\code{ate.est.train.group.cv}: }{A matrix of numerical values with
#'     \code{length(prop.multi) - 1} rows and \code{cv.n} columns.
#'     The ith column contains the estimated ATE in \code{length(prop.multi) - 1}
#'     mutually exclusive subgroups defined by \code{prop.multi} in the training set in ith
#'     cross-validation iteration.}
#'     \item{\code{ate.est.valid.group.cv}: }{Same as \code{ate.est.train.group.cv}, but in the
#'     validation set.}
#'     \item{\code{abc.valid}: }{A vector of numerical values of length \code{cv.n},
#'     The ith element returns the ABC of the validation curve in the ith cross-validation
#'     iteration. Only returned if \code{abc = TRUE}.}
#'  }
#'  \item{\code{ate.boosting}: }{A list of results similar to \code{ate.gaussian} output
#'  if \code{score.method} includes \code{'boosting'}.}
#'  \item{\code{ate.twoReg}: }{A list of results similar to \code{ate.gaussian} output
#'  if \code{score.method} includes \code{'twoReg'}.}
#'  \item{\code{ate.contrastReg}: }{A list of results similar to \code{ate.gaussian} output
#'  if \code{score.method} includes \code{'contrastReg'}.}
#'  \item{\code{ate.randomForest}: }{A list of results similar to \code{ate.gaussian} output
#'  if \code{score.method} includes \code{'randomForest'}.}
#'  \item{\code{ate.gam}: }{A list of results similar to \code{ate.gaussian} output
#'  if \code{score.method} includes \code{'gam'}.}
#'  \item{\code{props}: }{A list of 3 elements:}
#'  \itemize{
#'     \item{\code{prop.onlyhigh}: }{The original argument \code{prop.cutoff},
#'     reformatted as necessary.}
#'     \item{\code{prop.bi}: }{The original argument \code{prop.cutoff},
#'     similar to \code{prop.onlyhigh} but reformatted to exclude 1.}
#'     \item{\code{prop.multi}: }{The original argument \code{prop.multi},
#'     reformatted as necessary.}
#'  }
#'  \item{\code{overall.ate.train}: }{A vector of numerical values of length \code{cv.n}.
#'  The ith element contains the ATE in the training set of the ith cross-validation
#'  iteration, estimated with the doubly robust estimator.}
#'  \item{\code{overall.ate.valid}: }{A vector of numerical values of length \code{cv.n}.
#'  The ith element contains the ATE in the validation set of the ith cross-validation
#'  iteration, estimated with the doubly robust estimator.}
#'  \item{\code{higher.y}: }{The original \code{higher.y} argument.}
#'  \item{\code{abc}: }{The original \code{abc} argument.}
#'  \item{\code{cv.n}: }{The original \code{cv.n} argument.}
#'  \item{\code{response}: }{The type of response. Always 'continuous' for this function.}
#'  \item{\code{formulas}:}{A list of 3 elements: (1) \code{cate.model} argument,
#'  (2) \code{ps.model} argument and (3) original labels of the left-hand side variable in
#'  \code{ps.model} (treatment) if it was not 0/1.}
#' }
#'
#' @details The CATE score represents an individual-level treatment effect for continuous data,
#' estimated with boosting, linear regression, random forest, generalized additive model and the doubly
#' robust estimator (two regressions, Yadlowsky, 2020) applied separately by treatment group
#' or with the other doubly robust estimators (contrast regression, Yadlowsky, 2020) applied
#' to the entire data set.
#'
#' Internal CV is applied to reduce optimism in choosing the CATE estimation method that
#' captures the most treatment effect heterogeneity. The CV is applied by repeating the
#' following steps \code{cv.n} times:
#'
#' \enumerate{
#'  \item Split the data into a training and validation set according to \code{train.prop}.
#'  The training and validation sets must be balanced with respect to covariate distributions
#'  and doubly robust rate ratio estimates (see \code{error.max}).
#'  \item Estimate the CATE score in the training set with the specified scoring method.
#'  \item Predict the CATE score in the validation set using the scoring model fitted from
#'  the training set.
#'  \item Build nested subgroups of treatment responders in the training and validation sets,
#'  separately, and estimate the ATE within each nested subgroup. For each element i of
#'  \code{prop.cutoff} (e.g., \code{prop.cutoff[i]} = 0.6), take the following steps:
#'  \enumerate{
#'    \item Identify high responders as observations with the 60\%
#'    (i.e., \code{prop.cutoff[i]}x100\%) highest (if \code{higher.y = TRUE}) or
#'    lowest (if \code{higher.y = FALSE}) estimated CATE scores.
#'    \item Estimate the ATE in the subgroup of high responders using a doubly robust estimator.
#'    \item Conversely, identify low responders as observations with the 40\%
#'    (i.e., 1 - \code{prop.cutoff[i]}x100\%) lowest (if \code{higher.y} = TRUE) or
#'    highest (if \code{higher.y} = FALSE) estimated CATE scores.
#'    \item Estimate the ATE in the subgroup of low responders using a doubly robust estimator.
#'  }
#'  \item Build mutually exclusive subgroups of treatment responders in the training and
#'  validation sets, separately, and estimate the ATE within each subgroup. Mutually exclusive
#'  subgroups are built by splitting the estimated CATE scores according to \code{prop.multi}.
#'  \item If \code{abc} = TRUE, calculate the area between the ATE and the series of ATEs in
#'  nested subgroups of high responders in the validation set.
#' }
#'
#' @references Yadlowsky, S., Pellegrini, F., Lionetto, F., Braune, S., & Tian, L. (2020).
#' \emph{Estimation and validation of ratio-based conditional average treatment effects using
#' observational data. Journal of the American Statistical Association, 1-18.}
#' \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1772080}
#'
#' @seealso \code{\link{plot.precmed}()}, \code{\link{boxplot.precmed}()}, \code{\link{abc}()} methods for \code{"precmed"} objects,
#' and \code{\link{pmcount}()} function.
#'
#' @examples
#' \dontrun{
#' cv <- cvmean(cate.model = y ~ age  + previous_treatment + previous_cost + previous_status_measure,
#'              ps.model = trt ~ age + previous_treatment,
#'              data = meanExample,
#'              higher.y = FALSE,
#'              score.method = "",
#'              cv.n = 5,
#'              plot.gbmperf = FALSE,
#'              seed = 999)
#'
#' # Try:
#' # plot(x = cv, ylab = "Mean difference of drug1 vs drug0 in each subgroup")
#' # boxplot(x = cv, ylab = "mean difference of drug1 vs drug0 in each subgroup")
#' # abc(x = cv)
#' }
#' @export
#'
#' @importFrom graphics hist lines
#' @importFrom stats as.formula coef glm median model.frame model.matrix model.offset model.response na.omit optim pchisq predict qnorm quantile sd var
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom gbm gbm gbm.perf
#' @importFrom MESS auc
#' @importFrom stringr str_replace str_extract str_detect
#' @importFrom dplyr %>%
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom MASS glm.nb ginv
#' @importFrom randomForestSRC rfsrc predict.rfsrc

cvmean <- function(cate.model, init.model = NULL, ps.model, data, score.method,
                   higher.y = TRUE,
                   abc = TRUE,
                   prop.cutoff = seq(0.5, 1, length = 6),
                   prop.multi = c(0, 1/3, 2/3, 1),
                   ps.method = "glm", minPS = 0.01, maxPS = 0.99,
                   train.prop = 3/4, cv.n = 10,
                   error.max = 0.1, max.iter = 5000,
                   initial.predictor.method = "boosting", xvar.smooth.score = NULL, xvar.smooth.init = NULL,
                   tree.depth = 2, n.trees.rf = 1000, n.trees.boosting = 200, B = 3, Kfold = 6, plot.gbmperf = TRUE,
                   error.maxNR = 1e-3, tune = c(0.5, 2),
                   seed = NULL, verbose = 1, ...) {

  # TODO: now score.method has no default (mandatory argument)

  # Set seed for reproducibility
  set.seed(seed)

  if (verbose >= 1) t.start <- Sys.time()

  #### CHECK ARGUMENTS ####
  ##TODO: Insert n.trees.rf
  arg.checks(
    fun = "cv", response = "continuous", data = data, higher.y = higher.y, score.method = score.method, abc = abc,
    prop.cutoff = prop.cutoff, prop.multi = prop.multi,
    ps.method = ps.method, minPS = minPS, maxPS = maxPS,
    train.prop = train.prop, cv.n = cv.n,
    error.max = error.max, max.iter = max.iter,
    initial.predictor.method = initial.predictor.method,
    tree.depth = tree.depth, n.trees.boosting = n.trees.boosting, B = B, Kfold = Kfold, plot.gbmperf = plot.gbmperf,
    error.maxNR = error.maxNR, tune = tune
  )

  #### PRE-PROCESSING ####
  out <- data.preproc.mean(fun = "cv", cate.model = cate.model, init.model = init.model, ps.model = ps.model,
                           data = data, prop.cutoff = prop.cutoff, prop.multi = prop.multi, score.method = score.method,
                           ps.method = ps.method, initial.predictor.method = initial.predictor.method)
  y <- out$y
  trt <- out$trt
  x.ps <- out$x.ps
  x.cate <- out$x.cate
  prop.onlyhigh <- out$prop.onlyhigh
  prop.bi <- out$prop.bi
  prop.multi <- out$prop.multi
  initial.predictor.method <- out$initial.predictor.method

  if (any(score.method %in% c("contrastReg", "twoReg"))) x.init <- out$x.init
  if (abc == TRUE) prop.abc <- seq(prop.onlyhigh[1], prop.onlyhigh[length(prop.onlyhigh)], length = 100)

  #### FUNCTION STARTS HERE ####
  n.subgroup.onlyhigh <- length(prop.onlyhigh)
  n.subgroup.bi <- length(prop.bi)
  n.subgroup.multi <- length(prop.multi) - 1

  result <- vector("list", length(score.method) + 1)
  names(result) <- c(paste0("ate.", score.method), "props")
  for (name in paste0("ate.", score.method)) {
    result[[name]] <- vector("list", 6)
    names(result[[name]]) <- c("ate.est.train.high.cv", "ate.est.train.low.cv",
                               "ate.est.valid.high.cv", "ate.est.valid.low.cv",
                               "ate.est.train.group.cv", "ate.est.valid.group.cv")
    result[[name]]$ate.est.train.high.cv <-
      result[[name]]$ate.est.valid.high.cv <- matrix(NA, n.subgroup.onlyhigh, cv.n)
    result[[name]]$ate.est.train.low.cv <-
      result[[name]]$ate.est.valid.low.cv <- matrix(NA, n.subgroup.bi, cv.n)
    result[[name]]$ate.est.train.group.cv <-
      result[[name]]$ate.est.valid.group.cv <- matrix(NA, n.subgroup.multi, cv.n)

    # Add names to columns
    result[[name]] <- lapply(result[[name]], function(x) {
      colnames(x) <- paste0("cv", seq(1, cv.n))
      return(x)
    })

    # # Add names to rows
    rownames(result[[name]]$ate.est.train.high.cv) <-
      rownames(result[[name]]$ate.est.valid.high.cv) <- paste0("prop", round(prop.onlyhigh, 2))
    rownames(result[[name]]$ate.est.train.low.cv) <-
      rownames(result[[name]]$ate.est.valid.low.cv) <- paste0("prop", round(1 - prop.bi, 2))
    rownames(result[[name]]$ate.est.train.group.cv) <-
      rownames(result[[name]]$ate.est.valid.group.cv) <- paste0("prop", round(prop.multi[-1], 2))
  }
  result$props <- vector("list", 3)
  names(result$props) <- c("prop.onlyhigh", "prop.bi", "prop.multi")
  result$overall.ate.train <- result$overall.ate.valid <- rep(NA, cv.n)

  if ("contrastReg" %in% score.method) converge.contrastReg.cv <- rep(NA, cv.n)

  # Set progress bar
  if (verbose >= 1) pb <- txtProgressBar(min = 0, max = cv.n, style = 3)

  # Begin CV iteration
  for (cv.i in 1:cv.n) {

    ##### Split the data ------------------------------------------------------------------------------
    if (verbose >= 1) {
      cat("\ncv =", cv.i, "\n")
      cat("  splitting the data..\n")
    }

    datacv <- balancemean.split(y = y, trt = trt, x.cate = x.cate, x.ps = x.ps,
                                minPS = minPS, maxPS = maxPS, train.prop = train.prop,
                                error.max = error.max, max.iter = max.iter)

    y.train <- datacv$y.train
    trt.train <- datacv$trt.train
    x.cate.train <- datacv$x.cate.train
    x.ps.train <- datacv$x.ps.train

    y.valid <- datacv$y.valid
    trt.valid <- datacv$trt.valid
    x.cate.valid <- datacv$x.cate.valid
    x.ps.valid <- datacv$x.ps.valid

    if (any(score.method %in% c("contrastReg", "twoReg"))) {
      x.init.train <- x.init[-datacv$bestid.valid, , drop = FALSE]
      x.init.valid <- x.init[datacv$bestid.valid, , drop = FALSE]
    } else {
      x.init.train <- NA
      x.init.valid <- NA
    }

    #    if (any(score.method %in% c("contrastReg", "twoReg"))) {
    #      x.init.train <- x.init[-datacv$bestid.valid, , drop = FALSE]
    #      x.init.valid <- x.init[datacv$bestid.valid, , drop = FALSE]
    #    }

    # Calculate the ATE in the training and validation set (for abc, plot, boxplot)
    ## Changed: remove exp, drcount to drmean
    result$overall.ate.train[cv.i] <- drmean(y = y.train,
                                             x.cate = x.cate.train, x.ps = x.ps.train,
                                             trt = trt.train,
                                             ps.method = ps.method, minPS = minPS, maxPS = maxPS,
                                             interactions = FALSE)$mean.diff
    result$overall.ate.valid[cv.i] <- drmean(y = y.valid,
                                             x.cate = x.cate.valid, x.ps = x.ps.valid,
                                             trt = trt.valid,
                                             ps.method = ps.method, minPS = minPS, maxPS = maxPS,
                                             interactions = FALSE)$mean.diff


    ####### Fit the interaction model --------------------------------------------------------------
    if (verbose >= 1) cat("  training..\n")
    fit <- intxmean(y = y.train, trt = trt.train,
                    x.cate = x.cate.train, x.init = x.init.train, x.ps = x.ps.train,
                    score.method = score.method,
                    initial.predictor.method = initial.predictor.method,
                    xvar.smooth.init = xvar.smooth.init, xvar.smooth.score = xvar.smooth.score,
                    ps.method = ps.method, minPS = minPS, maxPS = maxPS,
                    tree.depth = tree.depth,  n.trees.rf = n.trees.rf, n.trees.boosting = n.trees.boosting,
                    Kfold = Kfold, B = B, plot.gbmperf = plot.gbmperf, ...)

    #    if ("contrastReg" %in% score.method) converge.contrastReg.cv[cv.i] <-
    #      fit$result.contrastReg$converge.contrastReg



    if (verbose == 2) {
      if (fit$best.iter == n.trees.boosting) {
        warning(paste("The best boosting iteration was iteration number", n.trees.boosting, " out of ", n.trees.boosting, ". Consider increasing the maximum number of trees and turning on boosting performance plot (plot.gbmperf = TRUE).", sep = ""))
      }

      # Check NA in coefficients of the score
      ##Changed poisson to gaussian
      if ("gaussian" %in% score.method & sum(is.na(fit$result.gaussian)) > 0) {
        warning("One or more coefficients in the score (Gaussian) are NA.
              Consider inspecting the distribution of the covariates in cate.model.")
      }
      if ("twoReg" %in% score.method & sum(is.na(fit$result.twoReg)) > 0) {
        warning("One or more coefficients in the score (two regressions) are NA.
              Consider inspecting the distribution of the covariates in cate.model.")
      }
      if ("contrastReg" %in% score.method & sum(is.na(fit$result.contrastReg)) > 0) {
        warning("One or more coefficients in the score (contrast regression) are NA.
              Consider inspecting the distribution of the covariates in cate.model.")
      }
    }


    if (verbose == 2) {
      if (length(names(fit$err.fit)) != 0) {
        if (verbose == 2) cat(paste0('    Error(s) occurred when fitting "', paste0(names(fit$err.fit), collapse = '", "'), '";\n    return NAs in the corresponding parameters'),'\n')
        warning(paste0('Error(s) occurred when fitting "', paste0(names(fit$err.fit), collapse = '", "'), '" in cross-validation iteration ', cv.i, ". NAs are returned for the corresponding parameters; see 'errors/warnings'."))
        result[['errors/warnings']][[names(fit$err.fit)]][[paste0("cv", cv.i)]]$errors <- fit$err.fit
      }

      if (length(names(fit$warn.fit)) != 0) {
        if (verbose == 2) cat(paste0('    Warning(s) occurred when fitting "', paste0(names(fit$warn.fit), collapse = '", "'), '"'),'\n')
        warning(paste0('Warning(s) occurred when fitting "', paste0(names(fit$warn.fit), collapse = '", "'), '" in cross-validation iteration ', cv.i, "; see 'errors/warnings'."))
        result[['errors/warnings']][[names(fit$warn.fit)]][[paste0("cv", cv.i)]]$warnings <- fit$warn.fit
      }
    }



    # Remove the method that produced the error, stop the analysis if all methods produced the error

    if (length(names(fit$err.fit)) > 0) {
      score.method.updated <- score.method[-which(score.method %in% names(fit$err.fit))]
    } else {score.method.updated <- score.method}
    #    score.method.updated <- score.method[-which(score.method %in% names(fit$err.fit))]
    if (length(score.method.updated) == 0) {stop("All methods produced error in fitting.")}

    ####### Construct the score in training set ------------------------------------------
    fit.score.train <- scoremean(fit = fit, x.cate = x.cate.train,
                                 score.method = score.method.updated)

    ####### Estimate the treatment effect in training set --------------------------------
    # if prop.onlyhigh ends with 1, estimate ATE once across all methods in the training
    if (prop.onlyhigh[n.subgroup.onlyhigh] == 1) {
      est.prop1 <- drmean(y = y.train,
                          x.cate = x.cate.train, x.ps = x.ps.train, trt = trt.train,
                          ps.method = ps.method, minPS = minPS, maxPS = maxPS, interactions = FALSE)$mean.diff
    } else {
      est.prop1 <- NULL
    }

    for (name in names(fit.score.train)) {
      score.train <- fit.score.train[[name]]
      ate.train.cv.i <-
        estmean.bilevel.subgroups(y = y.train,
                                  x.cate = x.cate.train, x.ps = x.ps.train,
                                  trt = trt.train,
                                  score = score.train, higher.y = higher.y,
                                  prop = prop.bi, onlyhigh = FALSE, ps.method = ps.method,
                                  minPS = minPS, maxPS = maxPS)
      result[[stringr::str_replace(name, "score", "ate")]]$ate.est.train.high.cv[, cv.i] <-
        c(ate.train.cv.i$ate.est.high, est.prop1)
      result[[stringr::str_replace(name, "score", "ate")]]$ate.est.train.low.cv[, cv.i] <-
        ate.train.cv.i$ate.est.low
      result[[stringr::str_replace(name, "score", "ate")]]$ate.est.train.group.cv[, cv.i] <-
        estmean.multilevel.subgroup(y = y.train,
                                    x.cate = x.cate.train, x.ps = x.ps.train,
                                    trt = trt.train,
                                    score = score.train, higher.y = higher.y,
                                    prop = prop.multi, ps.method = ps.method,
                                    minPS = minPS, maxPS = maxPS)
    }

    ####### Construct the score in validation set ----------------------------------------------
    if (verbose >= 1) cat("  validating..\n")

    fit.score.valid <- scoremean(fit = fit,
                                 x.cate = x.cate.valid,
                                 score.method = score.method.updated)

    ####### Estimate the treatment effect in validation set --------------------------------------
    # if prop.onlyhigh ends with 1, estimate ATE once across all methods in the validation
    if (prop.onlyhigh[n.subgroup.onlyhigh] == 1) {
      est.prop1 <- drmean(y = y.valid,
                          x.cate = x.cate.valid, x.ps = x.ps.valid, trt = trt.valid,
                          ps.method = ps.method, minPS = minPS, maxPS = maxPS, interactions = FALSE)$mean.diff
    }

    for (name in names(fit.score.valid)) {
      score.valid <- fit.score.valid[[name]]
      ate.valid.cv.i <-
        estmean.bilevel.subgroups(y = y.valid,
                                  x.cate = x.cate.valid, x.ps = x.ps.valid,
                                  trt = trt.valid,
                                  score = score.valid, higher.y = higher.y,
                                  prop = prop.bi, onlyhigh = FALSE,
                                  ps.method = ps.method, minPS = minPS, maxPS = maxPS)
      result[[stringr::str_replace(name, "score", "ate")]]$ate.est.valid.high.cv[, cv.i] <-
        c(ate.valid.cv.i$ate.est.high, est.prop1)
      result[[stringr::str_replace(name, "score", "ate")]]$ate.est.valid.low.cv[, cv.i] <-
        ate.valid.cv.i$ate.est.low
      result[[stringr::str_replace(name, "score", "ate")]]$ate.est.valid.group.cv[, cv.i] <-
        estmean.multilevel.subgroup(y = y.valid,
                                    x.cate = x.cate.valid, x.ps = x.ps.valid,
                                    trt = trt.valid,
                                    score = score.valid, higher.y = higher.y,
                                    prop = prop.multi,
                                    ps.method = ps.method, minPS = minPS, maxPS = maxPS)

      if (abc == TRUE) {
        temp.abc <- estmean.bilevel.subgroups(y = y.valid,
                                              x.cate = x.cate.valid, x.ps = x.ps.valid,
                                              trt = trt.valid,
                                              score = score.valid, higher.y = higher.y,
                                              prop = prop.abc, onlyhigh = TRUE,
                                              ps.method = ps.method, minPS = minPS, maxPS = maxPS)

        if (higher.y == TRUE) {
          # if higher scores are better, ABC requires subtracting the ATE from validation curve
          #ref <- temp.abc - log(result$overall.ate.valid[cv.i])
          ref <- temp.abc - (result$overall.ate.valid[cv.i])
        } else {
          #  if lower scores are better, ABC requires subtracting the validation curve from ATE
          #ref <- log(result$overall.ate.valid[cv.i]) - temp.abc
          ref <- (result$overall.ate.valid[cv.i]) - temp.abc
        }
        result[[stringr::str_replace(name, "score", "ate")]]$abc.valid[cv.i] <-
          MESS::auc(x = prop.abc, y = ref, from = prop.abc[1], to = prop.abc[length(prop.abc)], type = "spline")
      }

    }

    #    if (verbose >= 1 & ("contrastReg" %in% score.method)) {
    #      cat("  convergence", converge.contrastReg.cv[cv.i], "\n")
    #      cat("  ", date(), "\n")
    #      setTxtProgressBar(pb, cv.i)
    #    }

    # TODO: is it right?
    #    if (verbose >= 1 & ("contrastReg" %in% score.method.updated)) {
    #      cat("  convergence", converge.contrastReg.cv[cv.i], "\n")
    #      cat("  ", date(), "\n")
    #      setTxtProgressBar(pb, cv.i)
    #    }

  }
  if (verbose == 2) {
    if (any(is.na(unlist(result[stringr::str_replace(names(fit.score.train), "score", "ate")])))) {
      warning("Missing log rate ratio detected in the subgroups due to negative doubly robust estimator of
            the outcome for one or both treatment group(s).")
    }


  }

  if (verbose >= 1) {
    close(pb)
    t.end <- Sys.time()
    t.diff <- round(difftime(t.end, t.start),2)
    cat('Total runtime :',as.numeric(t.diff), attributes(t.diff)$units, '\n')
  }

  result$props$prop.onlyhigh <- prop.onlyhigh
  result$props$prop.bi <- prop.bi
  result$props$prop.multi <- prop.multi
  result$fgam <- fit$fgam.init
  result$higher.y <- higher.y
  result$abc <- abc
  result$cv.n <- cv.n
  result$response <- "continuous"
  result$formulas <- list(cate.model = cate.model, ps.model = ps.model, trt_labels = out$cat.trt)

  class(result) <- "precmed"

  return(result)
}



#' Estimation of the conditional average treatment effect (CATE) score for continuous data
#'
#' Provides singly robust and doubly robust estimation of CATE score with up to 6 scoring methods
#' among the following: Linear regression, boosting, two regressions, contrast regression, random forest and
#' generalized additive model.
#'
#' @param cate.model A formula describing the outcome model to be fitted.
#' The outcome must appear on the left-hand side.
#' @param init.model A formula describing the initial predictor model. The outcome must appear on the left-hand side.
#' It must be specified when \code{score.method = contrastReg} or \code{twoReg}.
#' @param ps.model A formula describing the propensity score model to be fitted.
#' The treatment must appear on the left-hand side.
#' The treatment must be a numeric vector coded as 0/1.
#' If data are from a RCT, specify \code{ps.model} as an intercept-only model.
#' @param data A data frame containing the variables in the outcome and propensity score models;
#' a data frame with \code{n} rows (1 row per observation).
#' @param score.method A vector of one or multiple methods to estimate the CATE score.
#' Allowed values are: \code{'boosting'}, \code{'gaussian'}, \code{'twoReg'},
#' \code{'contrastReg'}, \code{'randomForest'}, \code{'gam'}.
#' @param higher.y A logical value indicating whether higher (\code{TRUE}) or
#' lower (\code{FALSE}) values of the outcome are more desirable. Default is \code{TRUE}.
#' @param prop.cutoff A vector of numerical values (in (0, 1]) specifying percentiles of
#' the estimated log CATE scores to define nested subgroups. Each element represents the
#' cutoff to separate observations in nested subgroups (below vs above cutoff).
#' The length of \code{prop.cutoff} is the number of nested subgroups.
#' An equally-spaced sequence of proportions ending with 1 is recommended.
#' Default is \code{seq(0.5, 1, length = 6)}.
#' @param ps.method A character value for the method to estimate the propensity score.
#' Allowed values include one of: \code{'glm'} for logistic regression with main effects only
#' (default), or \code{'lasso'} for a logistic regression with main effects and LASSO penalization
#' on two-way interactions (added to the model if not specified in \code{ps.model}).
#' @param minPS A numerical value (in [0, 1]) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in (0, 1]) above which estimated propensity scores should be
#' truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.
#' @param initial.predictor.method A character vector for the method used to get initial outcome
#' predictions conditional on the covariates in \code{init.model} in
#' \code{score.method = 'twoReg'} and \code{'contrastReg'}. Allowed values include one of
#' \code{'gaussian'} (fastest), \code{'boosting'} and \code{'gam'}.
#' Default is \code{'boosting'}.
#' @param xvar.smooth.init A vector of characters indicating the name of the variables used as the
#' smooth terms if \code{initial.predictor.method = 'gam'}. The variables must be selected from
#' the variables listed in \code{init.model}. Default is \code{NULL}, which uses all variables
#' in \code{init.model}.
#' @param xvar.smooth.score A vector of characters indicating the name of the variables used as the
#' smooth terms if \code{score.method = 'gam'}. The variables must be selected from
#' the variables listed in \code{cate.model}. Default is \code{NULL}, which uses all variables
#' in \code{cate.model}.
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
#' @param plot.gbmperf A logical value indicating whether to plot the performance measures
#' in boosting. Used only if \code{score.method = 'boosting'} or if
#' \code{score.method = 'twoReg'} or \code{'contrastReg'} and
#' \code{initial.predictor.method = 'boosting'}. Default is \code{TRUE}.
#' @param error.maxNR A numerical value > 0 indicating the minimum value of the mean absolute
#' error in Newton Raphson algorithm. Used only if \code{score.method = 'contrastReg'}.
#' Default is \code{0.001}.
#' @param tune A vector of 2 numerical values > 0 specifying tuning parameters for the
#' Newton Raphson algorithm. \code{tune[1]} is the step size, \code{tune[2]} specifies a
#' quantity to be added to diagonal of the slope matrix to prevent singularity.
#' Used only if \code{score.method = 'contrastReg'}. Default is \code{c(0.5, 2)}.
#' @param seed An optional integer specifying an initial randomization seed for reproducibility.
#' Default is \code{NULL}, corresponding to no seed.
#' @param ... Additional arguments for \code{gbm()}
#'
#' @return Returns a list containing the following components:
#' \itemize{
#'  \item{\code{ate.gaussian}: }{A vector of numerical values of length \code{prop.cutoff}
#'  containing the estimated ATE in nested subgroups (defined by \code{prop.cutoff})
#'  constructed based on the estimated CATE scores with Poisson regression.
#'  Only provided if \code{score.method} includes \code{'gaussian'}.}
#'  \item{\code{ate.boosting}: }{Same as \code{ate.gaussian}, but with the nested subgroups based
#'  the estimated CATE scores with boosting. Only provided if \code{score.method}
#'  includes \code{'boosting'}.}
#'  \item{\code{ate.twoReg}: }{Same as \code{ate.gaussian}, but with the nested subgroups based
#'  the estimated CATE scores with two regressions.
#'  Only provided if \code{score.method} includes \code{'twoReg'}.}
#'  \item{\code{ate.contrastReg}: }{Same as \code{ate.gaussian}, but with the nested subgroups based
#'  the estimated CATE scores with contrast regression.
#'  Only provided if \code{score.method} includes \code{'contrastReg'}.}
#'  \item{\code{ate.randomForest}: }{Same as \code{ate.gaussian}, but with the nested subgroups based
#'  the estimated CATE scores with random forest.
#'  Only provided if \code{score.method} includes \code{'gam'}.}
#'  \item{\code{ate.gam}: }{Same as \code{ate.gaussian}, but with the nested subgroups based
#'  the estimated CATE scores with generalized additive model.
#'  Only provided if \code{score.method} includes \code{'gam'}.}
#'  \item{\code{score.gaussian}: }{A vector of numerical values of length n
#'  (number of observations in \code{data}) containing the estimated CATE scores
#'  according to the linear regression. Only provided if \code{score.method}
#'  includes \code{'gaussian'}.}
#'  \item{\code{score.boosting}: }{Same as \code{score.gaussian}, but with estimated CATE score
#'  according to boosting. Only provided if \code{score.method} includes
#'  \code{'boosting'}.}
#'  \item{\code{score.twoReg}: }{Same as \code{score.gaussian}, but with estimated CATE score
#'  according to two regressions. Only provided if \code{score.method} includes
#'  \code{'twoReg'}.}
#'  \item{\code{score.contrastReg}: }{Same as \code{score.gaussian}, but with estimated CATE score
#'  according to contrast regression. Only provided if \code{score.method} includes
#'  \code{'contrastReg'}.}
#'  \item{\code{score.randomForest}: }{Same as \code{score.gaussian}, but with estimated CATE score
#'  according to random forest. Only provided if \code{score.method}
#'  includes \code{'randomForest'}.}
#'  \item{\code{score.gam}: }{Same as \code{score.gaussian}, but with estimated CATE score
#'  according to generalized additive model. Only provided if \code{score.method}
#'  includes \code{'gam'}.}
#'  \item{\code{fit}: }{Additional details on model fitting if \code{score.method}
#'  includes 'boosting' or 'contrastReg':}
#'  \itemize{
#'    \item{\code{result.boosting}: }{Details on the boosting model fitted to observations
#'    with treatment = 0 \code{($fit0.boosting)} and to observations with treatment = 1 \code{($fit1.boosting)}.
#'    Only provided if \code{score.method} includes \code{'boosting'}.}
#'    \item{\code{result.randomForest}: }{Details on the boosting model fitted to observations
#'    with treatment = 0 \code{($fit0.randomForest)} and to observations with treatment = 1 \code{($fit1.randomForest)}.
#'    Only provided if \code{score.method} includes \code{'randomForest'}.}
#'    \item{\code{result.gam}: }{Details on the boosting model fitted to observations
#'    with treatment = 0 \code{($fit0.gam)} and to observations with treatment = 1 \code{($fit1.gam)}.
#'    Only provided if \code{score.method} includes \code{'gam'}.}
#'    \item{\code{result.contrastReg$sigma.contrastReg}: }{Variance-covariance matrix of
#'    the estimated CATE coefficients in contrast regression.
#'    Only provided if \code{score.method} includes \code{'contrastReg'}.}
#'  }
#'  \item{\code{coefficients}: }{A data frame with the coefficients of the estimated CATE
#'  score by \code{score.method}. The data frame has number of rows equal to the number of
#'  covariates in \code{cate.model} and number of columns equal to \code{length(score.method)}.
#'  If \code{score.method} includes \code{'contrastReg'}, the data frame has an additional
#'  column containing the standard errors of the coefficients estimated with contrast regression.
#'  \code{'boosting'}, \code{'randomForest'}, \code{'gam'} do not have coefficient results because these methods do not
#'  express the CATE as a linear combination of coefficients and covariates.}
#' }
#'
#' @details The CATE score represents an individual-level treatment effect, estimated with
#' either linear regression, boosting, random forest and generalized additive model applied separately by
#' treatment group or with two doubly robust estimators, two regressions and contrast regression
#' (Yadlowsky, 2020) applied to the entire dataset.
#'
#' \code{\link{pmmean}()} provides the coefficients of the CATE score for each scoring method requested
#' through \code{score.method}. Currently, contrast regression is the only method which allows
#' for inference of the CATE coefficients by providing standard errors of the coefficients.
#' The coefficients can be used to learn the effect size of each variable and predict the
#' CATE score for a new observation.
#'
#' \code{\link{pmmean}()} also provides the predicted CATE score of each observation in the data set,
#' for each scoring method. The predictions allow ranking the observations from potentially
#' high responders to the treatment to potentially low or standard responders.
#'
#' The estimated ATE among nested subgroups of high responders are also provided by scoring method.
#' Note that the ATEs in \code{\link{pmmean}()} are derived based on the CATE score which is estimated
#' using the same data sample. Therefore, overfitting may be an issue. \code{\link{pmmean}()} is more
#' suitable to inspect the estimated ATEs across scoring methods as it implements internal cross
#' validation to reduce optimism.
#'
#' @references Yadlowsky, S., Pellegrini, F., Lionetto, F., Braune, S., & Tian, L. (2020).
#' \emph{Estimation and validation of ratio-based conditional average treatment effects using
#' observational data. Journal of the American Statistical Association, 1-18.}
#' \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1772080}
#'
#' @seealso \code{\link{cvmean}()} function
#'
#' @examples
#'\dontrun{
#' pm <- pmmean(cate.model = y ~ age  + previous_treatment + previous_cost + previous_status_measure,
#'              init.model = y ~ age  + previous_treatment + previous_cost + previous_status_measure,
#'              ps.model = trt ~ age,
#'              data = meanExample,
#'              higher.y = FALSE,
#'              score.method = "gaussian",
#'              seed = 999)
#'}
#'@export

pmmean <- function(cate.model, init.model, ps.model, data, score.method,
                   higher.y = TRUE,
                   prop.cutoff = seq(0.5, 1, length = 6),
                   ps.method = "glm", minPS = 0.01, maxPS = 0.99,
                   initial.predictor.method = "boosting",
                   xvar.smooth.score = NULL, xvar.smooth.init = NULL,
                   tree.depth = 2, n.trees.rf = 1000, n.trees.boosting = 200, B = 3, Kfold = 6, plot.gbmperf = FALSE,
                   error.maxNR = 1e-3, tune = c(0.5, 2),
                   seed = NULL, ...) {

  # TODO: score.method is now a mandatory argument

  # Set seed once for reproducibility
  set.seed(seed)

  t.start <- Sys.time()

  #### CHECK ARGUMENTS ####
  arg.checks(
    fun = "pm", response = "continuous", data = data, higher.y = higher.y, score.method = score.method, prop.cutoff = prop.cutoff,
    ps.method = ps.method, minPS = minPS, maxPS = maxPS,
    initial.predictor.method = initial.predictor.method,
    tree.depth = tree.depth, n.trees.boosting = n.trees.boosting, B = B, Kfold = Kfold, plot.gbmperf = plot.gbmperf,
    error.maxNR = error.maxNR, tune = tune
  )

  #### PRE-PROCESSING ####
  out <- data.preproc.mean(fun = "pm", cate.model = cate.model, init.model = init.model, ps.model = ps.model,
                           score.method = score.method, data = data, prop.cutoff = prop.cutoff, ps.method = ps.method)
  y <- out$y
  trt <- out$trt
  x.ps <- out$x.ps
  x.cate <- out$x.cate
  prop <- out$prop

  if (any(score.method %in% c("contrastReg", "twoReg"))) x.init <- out$x.init

  #### FUNCTION STARTS HERE ####
  result <- vector("list", length(score.method) * 2 + 2)
  names(result) <- c(paste0("ate.", score.method),
                     paste0("score.", score.method), "fit", "coefficients")

  fit <- intxmean(y = y, trt = trt, x.cate = x.cate, x.init = x.init, x.ps = x.ps,
                  score.method = score.method,
                  ps.method = ps.method, minPS = minPS, maxPS = maxPS,
                  initial.predictor.method = initial.predictor.method,
                  xvar.smooth.init = xvar.smooth.init, xvar.smooth.score = xvar.smooth.score,
                  tree.depth = tree.depth, n.trees.rf = n.trees.rf, n.trees.boosting = n.trees.boosting,
                  Kfold = Kfold, B = B, plot.gbmperf = plot.gbmperf,
                  error.maxNR = error.maxNR, tune = tune, ...)

  if (fit$best.iter == n.trees.boosting) {
    warning(paste("The best boosting iteration was iteration number", n.trees.boosting, " out of ", n.trees.boosting, ". Consider increasing the maximum number of trees and turning on boosting performance plot (plot.gbmperf = TRUE).", sep = ""))
  }

  # Check NA in coefficients of the score
  if ("gaussian" %in% score.method & sum(is.na(fit$result.gaussian)) > 0) {
    warning("One or more coefficients in the score (Gaussian) are NA.
            Consider inspecting the distribution of the covariates in cate.model.")
  }
  if ("twoReg" %in% score.method & sum(is.na(fit$result.twoReg)) > 0) {
    warning("One or more coefficients in the score (two regressions) are NA.
            Consider inspecting the distribution of the covariates in cate.model.")
  }
  if ("contrastReg" %in% score.method & sum(is.na(fit$result.contrastReg)) > 0) {
    warning("One or more coefficients in the score (contrast regression) are NA.
            Consider inspecting the distribution of the covariates in cate.model.")
  }

  if (length(names(fit$err.fit)) > 0) {
    score.method.updated <- score.method[-which(score.method %in% names(fit$err.fit))]
  } else {score.method.updated <- score.method}
  #    score.method.updated <- score.method[-which(score.method %in% names(fit$err.fit))]
  if (length(score.method.updated) == 0) {stop("All methods produced error in fitting.")}

  fit.score <- scoremean(fit = fit,
                         x.cate = x.cate,
                         score.method = score.method.updated)

  for (name in names(fit.score)) {
    score <- fit.score[[name]]
    result[[name]] <- score
    ate <- estmean.bilevel.subgroups(y = y,
                                     x.cate = x.cate, x.ps = x.ps,
                                     trt = trt,
                                     score = score, higher.y = higher.y,
                                     prop = prop, onlyhigh = TRUE,
                                     ps.method = ps.method, minPS = minPS, maxPS = maxPS)
    result[[str_replace(name, "score", "ate")]] <- ate
    names(result[[str_replace(name, "score", "ate")]]) <- paste0("prop", round(prop, 2))
  }

  if (sum(score.method %in% c("gaussian", "twoReg", "contrastReg")) > 0) {
    cf <- data.frame(matrix(NA, nrow = ncol(x.cate) + 1, ncol = 4))
    colnames(cf) <- c("gaussian", "twoReg", "contrastReg", "SE_contrastReg")
    rownames(cf) <- c("(Intercept)", colnames(x.cate))

    if ("gaussian" %in% score.method) cf$gaussian <- fit$result.gaussian

    if ("twoReg" %in% score.method) cf$twoReg <- fit$result.twoReg

    if ("contrastReg" %in% score.method) {
      cf$contrastReg <- fit$result.contrastReg$delta.contrastReg
      cf$SE_contrastReg <- sqrt(diag(fit$result.contrastReg$sigma.contrastReg))
    }

    result$coefficients <- cf[, colSums(is.na(cf)) != nrow(cf), drop = FALSE]
  }


  if (any(is.na(unlist(result[str_replace(names(fit.score), "score", "ate")])))) {
    warning("Missing log rate ratio detected due to negative doubly robust estimator of y|x
            for one or both treatment group(s).")
  }

  if ("boosting" %in% score.method) result$fit$result.boosting <-
    fit$result.boosting
  if ("randomForest" %in% score.method) result$fit$result.randomForest <-
    fit$result.randomForest
  if ("gam" %in% score.method) result$fit$result.gam <-
    fit$result.gam
  if ("contrastReg" %in% score.method) result$fit$result.contrastReg$sigma.contrastReg <-
    fit$result.contrastReg$sigma.contrastReg

  t.end <- Sys.time()
  t.diff <- round(difftime(t.end, t.start),2)
  cat('Total runtime :',as.numeric(t.diff), attributes(t.diff)$units, '\n')

  class(result) <- "precmed"
  return(result)
}


#' Doubly robust estimator of and inference for the average treatment effect for continuous data
#'
#' Doubly robust estimator of the average treatment effect between two treatments, which is the rate ratio
#' of treatment 1 over treatment 0 for count outcomes. Bootstrap is used for inference.
#'
#' @param cate.model A formula describing the outcome model to be fitted.
#' The outcome must appear on the left-hand side.
#' @param ps.model A formula describing the propensity score model to be fitted.
#' The treatment must appear on the left-hand side. The treatment must be a numeric vector
#' coded as 0/1. If data are from a RCT, specify \code{ps.model} as an intercept-only model.
#' @param data A data frame containing the variables in the outcome and propensity score models;
#' a data frame with \code{n} rows (1 row per observation).
#' @param ps.method A character value for the method to estimate the propensity score. Allowed values include one of:
#' \code{'glm'} for logistic regression with main effects only (default), or
#' \code{'lasso'} for a logistic regression with main effects and LASSO penalization on
#' two-way interactions (added to the model if interactions are not specified in \code{ps.model}).
#' Relevant only when \code{ps.model} has more than one variable.
#' @param minPS A numerical value (in [0, 1]) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in (0, 1]) above which estimated propensity scores should be
#' truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.
#' @param interactions A logical value indicating whether the outcome model should be fitted separately by treatment arm
#' with the variables in \code{cate.model}, which is equivalent to assuming interaction terms between \code{trt} and
#' all the variables in \code{cate.model} in a model fitted with both treatment arms. If \code{TRUE}, the outcome
#' model will be fitted separately by treatment arms only if at least 10 patients received each treatment option.
#' Default is \code{TRUE}.
#' @param n.boot A numeric value indicating the number of bootstrap samples used. Default is \code{500}.
#' @param verbose An integer value indicating whether intermediate progress messages and histograms should
#' be printed. \code{1} indicates messages are printed and \code{0} otherwise. Default is \code{1}.
#' @param plot.boot A logical value indicating whether histograms of the bootstrapped log(rate ratio) should
#' be produced at every \code{n.boot/10}-th iteration and whether the final histogram should be outputted.
#' Default is \code{FALSE}.
#' @param seed An optional integer specifying an initial randomization seed for reproducibility.
#' Default is \code{NULL}, corresponding to no seed.
#'
#' @return Return a list of 8 elements:
#' \itemize{
#'   \item{\code{log.rate.ratio}: } A numeric value of the estimated log rate ratio.
#'   \item{\code{se.boot.log.rate.ratio}: } A numeric value of the bootstrap standard error of log rate ratio.
#'   \item{\code{rate.ratio}: } A numeric value of the estimated rate ratio.
#'   \item{\code{rate.ratio0}: } A numeric value of the estimated rate in the group trt=0.
#'   \item{\code{rate.ratio1}: } A numeric value of the estimated rate in the group trt=1.
#'   \item{\code{rate.ratio.CIl}: } A numeric value of the lower limit 95% bootstrap confidence interval
#'     for estimated rate ratio.
#'   \item{\code{rate.ratio.CIu}: } A numeric value of the upper limit 95% bootstrap confidence interval
#'     for estimated rate ratio.
#'   \item{\code{pvalue}: } A numeric value of the p-value derived from the bootstrapped values
#'     based on a Chi-squared distribution.
#'   \item{\code{warning}: } A warning message produced if the treatment variable was not coded as 0/1. The key
#'   to map the original coding of the variable to a 0/1 key is displayed in the warning to facilitate the
#'   interpretation of the remaining of the output.
#'   \item{\code{plot}: } If \code{plot.boot} is \code{TRUE}, a histogram displaying the distribution of the bootstrapped log rate ratios.
#'   The red vertical reference line in the histogram represents the estimated log rate ratio.
#' }
#'
#' @details This helper function estimates the average treatment effect (ATE) between two
#'  treatment groups in a given dataset specified by \code{y, trt, x.cate, x.ps, time}. The ATE is
#'  estimated with a doubly robust estimator that accounts for imbalances in covariate distributions
#'  between the two treatment groups with inverse probability treatment weighting.
#'  For count outcomes, the estimated ATE is the estimated
#'  rate ratio between treatment 1 versus treatment 0. Both original and log-transformed ATEs are
#'  returned, as well as the rate in either treatment group.
#'  If \code{inference = TRUE}, the variability of the estimated rate ratio is also calculated
#'  using bootstrap. Additional variability outputs include standard error of the log rate ratio,
#'  95% confidence interval of the rate ratio, p-value, and a histogram of the log rate ratio.
#'
#' @examples
#'\dontrun{
#' output <- drmean.inference(cate.model = y ~ age + female + previous_treatment +
#'                                previous_cost + previous_number_relapses,
#'                             ps.model = trt ~ age + previous_treatment,
#'                             data = countExample,
#'                             plot.boot = TRUE,
#'                             seed = 999)
#' print(output)
#' output$plot
#'}
#'
#' @importFrom ggplot2 ggplot geom_histogram geom_vline
#' @importFrom dplyr mutate
#' @export

drmean.inference <- function(cate.model, ps.model, data,
                             ps.method = "glm", minPS = 0.01, maxPS = 0.99,
                             interactions = TRUE, n.boot = 500, verbose = 1, plot.boot = FALSE, seed = NULL) {

  # Set seed once for reproducibility
  set.seed(seed)

  #### CHECK ARGUMENTS ####
  # TODO: Phoebe I think we removed the arg.check on interactions but here it is?
  arg.checks(fun = "drinf", response = "continuous", data = data,
             ps.method = ps.method,
             minPS = minPS,
             maxPS = maxPS,
             interactions = interactions,
             n.boot = n.boot,
             plot.boot = plot.boot)

  #### PRE-PROCESSING ####
  preproc <- data.preproc.mean(fun = "drinf", cate.model = cate.model,
                               init.model = NULL, #TODO: Needs to be specified
                               ps.model = ps.model,
                               data = data, ps.method = ps.method)
  y <- preproc$y
  trt <- preproc$trt
  x.ps <- preproc$x.ps
  x.cate <- preproc$x.cate

  # Point estimate
  est <- drmean(y = y, x.cate = x.cate, x.ps = x.ps, trt = trt, ps.method = ps.method, minPS = minPS, maxPS = maxPS, interactions = interactions)
  meandiff <- est$mean.diff

  # Apply bootstrap
  n <- length(y)
  if (is.na(meandiff) == TRUE) {
    stop("Impossible to use bootstrap, Mean Difference is NA.")
  } else {
    trt.boot <- rep(0, n.boot)
    for (i in 1:n.boot) {
      idsub.boot <- sample(n, size = n, replace = TRUE)
      trt.boot[i] <- drmean(y = y[idsub.boot], x.cate = x.cate[idsub.boot, , drop = FALSE], x.ps = x.ps[idsub.boot, , drop = FALSE], trt = trt[idsub.boot], ps.method = ps.method, minPS = minPS, maxPS = maxPS, interactions = interactions)$mean.diff

      if (i %% (n.boot %/% 10) == 0) {
        if (verbose == 1) cat("Bootstrap iteration", i, "\n")
        if (plot.boot == TRUE) {
          p <- trt.boot %>%
            as.data.frame() %>%
            ggplot(aes(x = .data$.)) +
            geom_histogram(bins = 50) +
            theme_classic() +
            geom_vline(xintercept = meandiff, linetype = 1, color = "red") +
            labs(x = "Bootstrap values of log rate ratio", y = "Frequency", title = paste0(i, " bootstrap iterations"))
          print(p)
        } # end of if (plot.boot == TRUE)
      } # end of if (i %% (n.boot %/% 10) == 0)
    } # end of for (i in 1:n.boot)

    out <- c()
    out$mean.diff <- meandiff
    out$se.mean.diff <- se.est <- sd(trt.boot, na.rm = TRUE)
    #out$rate.ratio <- exp(logrr)
    out$mean.diff0 <- est$mean.diff0
    out$mean.diff1 <- est$mean.diff1
    out$mean.diff.CIl <- meandiff - qnorm(0.975) * se.est
    out$mean.diff.CIu <- meandiff + qnorm(0.975) * se.est
    out$pval <- 1 - pchisq(meandiff^2 / se.est^2, 1) #TODO: modify this..can we use chisq?
    out$warning <- preproc$warning

    if (plot.boot == TRUE) {
      plot <- trt.boot %>% as.data.frame() %>%
        ggplot(aes(x = .data$.)) +
        geom_histogram(bins = 50) +
        theme_classic() +
        geom_vline(xintercept = meandiff, linetype = 1, color = "red") +
        labs(x = "Bootstrap values of mean difference", y = "Frequency", title = paste0(n.boot, " bootstrap iterations"))
      out$plot <- plot
    }
  }

  return(out)
}
