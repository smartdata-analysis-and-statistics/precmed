# ------------------------------------------------------------------
#
# Project: Precision Medicine MS (precmed) - Comprehensive R package
#
# Purpose: Output functions for Count outcomes
#
# Platform: Windows
# R Version: 4.1.0
#



#' Cross-validation of the conditional average treatment effect (CATE) score for count outcomes
#'
#' Provides doubly robust estimation of the average treatment effect (ATE) in nested and
#' mutually exclusive subgroups of patients defined by an estimated conditional average
#' treatment effect (CATE) score via cross-validation (CV). The CATE score can be estimated
#' with up to 5 methods among the following: Poisson regression, boosting, two regressions,
#' contrast regression, and negative binomial (see \code{score.method}).
#'
#' @param cate.model A formula describing the outcome model to be fitted.
#' The outcome must appear on the left-hand side.
#' @param ps.model A formula describing the propensity score model to be fitted.
#' The treatment must appear on the left-hand side. The treatment must be a numeric vector
#' coded as 0/1. If data are from a randomized trial, specify \code{ps.model} as an intercept-only model.
#' @param data A data frame containing the variables in the outcome and propensity score model; a data frame with \code{n} rows
#' (1 row per observation).
#' @param score.method A vector of one or multiple methods to estimate the CATE score.
#' Allowed values are: \code{'boosting'}, \code{'poisson'}, \code{'twoReg'}, \code{'contrastReg'},
#' and \code{'negBin'}.
#' @param higher.y A logical value indicating whether higher (\code{TRUE}) or
#' lower (\code{FALSE}) values of the outcome are more desirable. Default is \code{TRUE}.
#' @param abc A logical value indicating whether the area between curves (ABC) should be calculated
#' at each cross-validation iterations, for each \code{score.method}. Default is \code{TRUE}.
#' @param prop.cutoff A vector of numerical values (in (0, 1]) specifying percentiles of the
#' estimated log CATE scores to define nested subgroups. Each element represents the cutoff to
#' separate observations in nested subgroups (below vs above cutoff).
#' The length of \code{prop.cutoff} is the number of nested subgroups.
#' An equally-spaced sequence of proportions ending with 1 is recommended.
#' Default is \code{seq(0.5, 1, length = 6)}.
#' @param prop.multi A vector of numerical values (in [0, 1]) specifying percentiles of the
#' estimated log CATE scores to define mutually exclusive subgroups.
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
#' outcome predictions conditional on the covariates in \code{cate.model}. Only applies
#' when \code{score.method} includes \code{'twoReg'} or \code{'contrastReg'}. Allowed values
#' include one of \code{'poisson'} (fastest), \code{'boosting'} and \code{'gam'}.
#' Default is \code{'boosting'}.
#' @param xvar.smooth A vector of characters indicating the name of the variables used as
#' the smooth terms if \code{initial.predictor.method = 'gam'}. The variables must be selected
#' from the variables listed in \code{cate.model}.
#' Default is \code{NULL}, which uses all variables in \code{cate.model}.
#' @param tree.depth A positive integer specifying the depth of individual trees in boosting
#' (usually 2-3). Used only if \code{score.method = 'boosting'} or
#' if \code{initial.predictor.method = 'boosting'} with \code{score.method = 'twoReg'} or
#' \code{'contrastReg'}. Default is 2.
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
#' @param plot.gbmperf A logical value indicating whether to plot the performance measures in
#' boosting. Used only if \code{score.method = 'boosting'} or if \code{score.method = 'twoReg'}
#' or \code{'contrastReg'} and \code{initial.predictor.method = 'boosting'}. Default is \code{TRUE}.
#' @param verbose An integer value indicating what kind of intermediate progress messages should
#' be printed. \code{0} means no outputs. \code{1} means only progress bar and run time.
#' \code{2} means progress bar, run time, and all errors and warnings. Default is \code{0}.
#' @param ... Additional arguments for \code{gbm()}
#'
#' @return Returns a list containing the following components saved as a \code{"precmed"} object:
#' \itemize{
#'  \item{\code{ate.poisson}: }{A list of results output if \code{score.method} includes
#'  \code{'poisson'}:}
#'  \itemize{
#'     \item{\code{ate.est.train.high.cv}: }{A matrix of numerical values with
#'     \code{length(prop.cutoff)} rows and \code{cv.n} columns.
#'     The ith row/jth column cell contains the estimated ATE in the nested subgroup of high responders
#'     defined by CATE score above (if \code{higher.y = TRUE}) or below (if \code{higher.y = FALSE}) the
#'     \code{prop.cutoff[i]}x100\% percentile of the estimated CATE score in the training set in the jth
#'     cross-validation iteration.}
#'     \item{\code{ate.est.train.low.cv}: }{A matrix of numerical values with
#'     \code{length(prop.cutoff) - 1} rows and \code{cv.n} columns.
#'     The ith row/jth column cell contains the estimated ATE in the nested subgroup of low responders
#'     defined by CATE score below (if \code{higher.y = TRUE}) or above (if \code{higher.y = FALSE}) the
#'     \code{prop.cutoff[i]}x100\% percentile of the estimated CATE score in the training set in the jth
#'     cross-validation iteration.}
#'     \item{\code{ate.est.valid.high.cv}: }{Same as \code{ate.est.train.high.cv},
#'     but in the validation set.}
#'     \item{\code{ate.est.valid.low.cv}: }{Same as \code{ate.est.train.low.cv},
#'     but in the validation set.}
#'     \item{\code{ate.est.train.group.cv}: }{A matrix of numerical values with
#'     \code{length(prop.multi) - 1} rows and \code{cv.n} columns.
#'     The jth column contains the estimated ATE in \code{length(prop.multi) - 1}
#'     mutually exclusive subgroups defined by \code{prop.multi} in the training set in jth
#'     cross-validation iteration.}
#'     \item{\code{ate.est.valid.group.cv}: }{Same as \code{ate.est.train.group.cv}, but in the
#'     validation set.}
#'     \item{\code{abc.valid}: }{A vector of numerical values of length \code{cv.n}.
#'     The ith element returns the ABC of the validation curve in the ith cross-validation
#'     iteration. Only returned if \code{abc = TRUE}.}
#'  }
#'  \item{\code{ate.boosting}: }{A list of results similar to \code{ate.poisson} output
#'  if \code{score.method} includes \code{'boosting'}.}
#'  \item{\code{ate.twoReg}: }{A list of results similar to \code{ate.poisson} output
#'  if \code{score.method} includes \code{'twoReg'}.}
#'  \item{\code{ate.contrastReg}: }{A list of results similar to \code{ate.poisson} output
#'  if \code{score.method} includes \code{'contrastReg'}.
#'  This method has an additional element in the list of results:}
#'  \itemize{
#'     \item{\code{converge.contrastReg.cv}: }{A vector of logical value of length \code{cv.n}.
#'     The ith element indicates whether the algorithm converged in the ith cross-validation
#'     iteration.}
#'  }
#'  \item{\code{ate.negBin}: }{A list of results similar to \code{ate.poisson} output
#'  if \code{score.method} includes \code{'negBin'}.}
#'  \item{\code{props}: }{A list of 3 elements:}
#'  \itemize{
#'     \item{\code{prop.onlyhigh}: }{The original argument \code{prop.cutoff},
#'     reformatted as necessary.}
#'     \item{\code{prop.bi}: }{The original argument \code{prop.cutoff},
#'     similar to \code{prop.onlyhigh} but reformatted to exclude 1.}
#'     \item{\code{prop.multi}: }{The original argument \code{prop.multi},
#'     reformatted as necessary to include 0 and 1.}
#'  }
#'  \item{\code{overall.ate.valid}: }{A vector of numerical values of length \code{cv.n}.
#'  The ith element contains the ATE in the validation set of the ith cross-validation
#'  iteration, estimated with the doubly robust estimator.}
#'  \item{\code{overall.ate.train}: }{A vector of numerical values of length \code{cv.n}.
#'  The ith element contains the ATE in the training set of the ith cross-validation
#'  iteration, estimated with the doubly robust estimator.}
#'  \item{\code{fgam}: }{The formula used in GAM if \code{initial.predictor.method = 'gam'}.}
#'  \item{\code{higher.y}: }{The original \code{higher.y} argument.}
#'  \item{\code{abc}: }{The original \code{abc} argument.}
#'  \item{\code{cv.n}: }{The original \code{cv.n} argument.}
#'  \item{\code{response}: }{The type of response. Always 'count' for this function.}
#'  \item{\code{formulas}:}{A list of 3 elements: (1) \code{cate.model} argument,
#'  (2) \code{ps.model} argument and (3) original labels of the left-hand side variable in
#'  \code{ps.model} (treatment) if it was not 0/1.}
#' }
#'
#' @details The CATE score represents an individual-level treatment effect expressed as a
#' rate ratio for count outcomes. It can be estimated with boosting, Poisson regression,
#' negative binomial regression, and the doubly robust estimator two regressions (Yadlowsky,
#' 2020) applied separately by treatment group or with the other doubly robust estimator
#' contrast regression (Yadlowsky, 2020) applied to the entire data set.
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
#'  \item If \code{abc} = TRUE, calculate the area between the ATE and the series of ATEs in
#'  nested subgroups of high responders in the validation set.
#'  \item Build mutually exclusive subgroups of treatment responders in the training and
#'  validation sets, separately, and estimate the ATE within each subgroup. Mutually exclusive
#'  subgroups are built by splitting the estimated CATE scores according to \code{prop.multi}.
#' }
#'
#' @references Yadlowsky, S., Pellegrini, F., Lionetto, F., Braune, S., & Tian, L. (2020).
#' \emph{Estimation and validation of ratio-based conditional average treatment effects using
#' observational data. Journal of the American Statistical Association, 1-18.}
#' \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1772080}
#'
#' @seealso \code{\link{plot.precmed}()}, \code{\link{boxplot.precmed}()}, \code{\link{abc}()} methods for \code{"precmed"} objects,
#' and \code{\link{catefitcount}()} function.
#'
#' @examples
#' catecv <- catecvcount(cate.model = y ~ age + female + previous_treatment +
#'                                previous_cost + previous_number_relapses + offset(log(years)),
#'               ps.model = trt ~ age + previous_treatment,
#'               data = countExample,
#'               higher.y = FALSE,
#'               score.method = "poisson",
#'               cv.n = 5,
#'               plot.gbmperf = FALSE,
#'               seed = 999)
#'
#' # Try:
#' # plot(x = catecv, ylab = "Rate ratio of drug1 vs drug0 in each subgroup")
#' # boxplot(x = catecv, ylab = "Rate ratio of drug1 vs drug0 in each subgroup")
#' # abc(x = catecv)
#'
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
#' @importFrom MASS glm.nb

catecvcount <- function(data,
                        score.method, # Mandatory argument that needs to be specified
                        cate.model,
                        ps.model,
                        ps.method = "glm",
                        initial.predictor.method = "boosting",
                        minPS = 0.01,
                        maxPS = 0.99,
                        verbose = 0,
                        higher.y = TRUE,
                        prop.cutoff = seq(0.5, 1, length = 6),
                        prop.multi = c(0, 1/3, 2/3, 1),
                        abc = TRUE,
                        train.prop = 3/4,
                        cv.n = 10,
                        error.max = 0.1,
                        max.iter = 5000,
                        xvar.smooth = NULL,
                        tree.depth = 2,
                        n.trees.boosting = 200,
                        B = 3,
                        Kfold = 5,
                        error.maxNR = 1e-3,
                        max.iterNR = 150,
                        tune = c(0.5, 2),
                        seed = NULL,
                        plot.gbmperf = TRUE,
                        ...) {

  # Set seed for reproducibility
  set.seed(seed)

  if (verbose >= 1) t.start <- Sys.time()

  #### CHECK ARGUMENTS ####
  arg.checks(
    fun = "cv", response = "count", data = data, higher.y = higher.y, score.method = score.method, abc = abc,
    prop.cutoff = prop.cutoff, prop.multi = prop.multi,
    ps.method = ps.method, minPS = minPS, maxPS = maxPS,
    train.prop = train.prop, cv.n = cv.n,
    error.max = error.max, max.iter = max.iter,
    initial.predictor.method = initial.predictor.method,
    tree.depth = tree.depth, n.trees.boosting = n.trees.boosting, B = B, Kfold = Kfold, plot.gbmperf = plot.gbmperf,
    error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune
  )

  #### PRE-PROCESSING ####
  out <- data.preproc(fun = "cv", cate.model = cate.model, ps.model = ps.model,
                      data = data, prop.cutoff = prop.cutoff, prop.multi = prop.multi,
                      ps.method = ps.method, initial.predictor.method = initial.predictor.method)
  y <- out$y
  trt <- out$trt
  x.ps <- out$x.ps
  x.cate <- out$x.cate
  time <- out$time
  prop.onlyhigh <- out$prop.onlyhigh
  prop.bi <- out$prop.bi
  prop.multi <- out$prop.multi
  initial.predictor.method <- out$initial.predictor.method

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

    datacv <- balance.split(y = y, trt = trt, x.cate = x.cate, x.ps = x.ps, time = time,
                            minPS = minPS, maxPS = maxPS, train.prop = train.prop,
                            error.max = error.max, max.iter = max.iter)

    y.train <- datacv$y.train
    trt.train <- datacv$trt.train
    x.cate.train <- datacv$x.cate.train
    x.ps.train <- datacv$x.ps.train
    time.train <- datacv$time.train

    y.valid <- datacv$y.valid
    trt.valid <- datacv$trt.valid
    x.cate.valid <- datacv$x.cate.valid
    x.ps.valid <- datacv$x.ps.valid
    time.valid <- datacv$time.valid

    # Calculate the ATE in the training and validation set (for abc, plot, boxplot)
    result$overall.ate.train[cv.i] <- exp(drcount(y = y.train,
                                                  x.cate = x.cate.train, x.ps = x.ps.train,
                                                  trt = trt.train, time = time.train,
                                                  ps.method = ps.method, minPS = minPS, maxPS = maxPS,
                                                  interactions = FALSE)$log.rate.ratio)
    result$overall.ate.valid[cv.i] <- exp(drcount(y = y.valid,
                                                  x.cate = x.cate.valid, x.ps = x.ps.valid,
                                                  trt = trt.valid, time = time.valid,
                                                  ps.method = ps.method, minPS = minPS, maxPS = maxPS,
                                                  interactions = FALSE)$log.rate.ratio)


    ####### Fit the interaction model --------------------------------------------------------------
    if (verbose >= 1) cat("  training..\n")
    fit <- intxcount(y = y.train, trt = trt.train,
                     x.cate = x.cate.train, x.ps = x.ps.train,
                     time = time.train,
                     score.method = score.method,
                     initial.predictor.method = initial.predictor.method,
                     xvar.smooth = xvar.smooth,
                     ps.method = ps.method, minPS = minPS, maxPS = maxPS,
                     tree.depth = tree.depth, n.trees.boosting = n.trees.boosting,
                     Kfold = Kfold, B = B, plot.gbmperf = plot.gbmperf,
                     error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune, ...)

    if ("contrastReg" %in% score.method) converge.contrastReg.cv[cv.i] <-
      fit$result.contrastReg$converge.contrastReg

    if (verbose == 2){
      if (fit$best.iter == n.trees.boosting) {
        warning(paste("The best boosting iteration was iteration number", n.trees.boosting, " out of ", n.trees.boosting, ". Consider increasing the maximum number of trees and turning on boosting performance plot (plot.gbmperf = TRUE).", sep = ""))
      }

      # Check NA in coefficients of the score
      if ("poisson" %in% score.method & sum(is.na(fit$result.poisson)) > 0) {
        warning("One or more coefficients in the score (Poisson) are NA.
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
      if ("negBin" %in% score.method & sum(is.na(fit$result.negBin)) > 0) {
        warning("One or more coefficients in the score (negative binomial) are NA.
              Consider inspecting the distribution of the covariates in cate.model.")
      }
    }


    ####### Construct the score in training set ------------------------------------------
    fit.score.train <- scorecount(fit = fit, x.cate = x.cate.train,
                                  time = time.train, score.method = score.method)

    ####### Estimate the treatment effect in training set --------------------------------
    # if prop.onlyhigh ends with 1, estimate ATE once across all methods in the training
    if (prop.onlyhigh[n.subgroup.onlyhigh] == 1) {
      est.prop1 <- exp(drcount(y = y.train,
                               x.cate = x.cate.train, x.ps = x.ps.train, trt = trt.train,
                               time = time.train,
                               ps.method = ps.method, minPS = minPS, maxPS = maxPS, interactions = FALSE)$log.rate.ratio)
    } else {
      est.prop1 <- NULL
    }

    for (name in names(fit.score.train)) {
      score.train <- fit.score.train[[name]]
      ate.train.cv.i <-
        estcount.bilevel.subgroups(y = y.train,
                                   x.cate = x.cate.train, x.ps = x.ps.train,
                                   time = time.train, trt = trt.train,
                                   score = score.train, higher.y = higher.y,
                                   prop = prop.bi, onlyhigh = FALSE, ps.method = ps.method,
                                   minPS = minPS, maxPS = maxPS)
      result[[str_replace(name, "score", "ate")]]$ate.est.train.high.cv[, cv.i] <-
        c(ate.train.cv.i$ate.est.high, est.prop1)
      result[[str_replace(name, "score", "ate")]]$ate.est.train.low.cv[, cv.i] <-
        ate.train.cv.i$ate.est.low
      result[[str_replace(name, "score", "ate")]]$ate.est.train.group.cv[, cv.i] <-
        estcount.multilevel.subgroup(y = y.train,
                                     x.cate = x.cate.train, x.ps = x.ps.train,
                                     time = time.train, trt = trt.train,
                                     score = score.train, higher.y = higher.y,
                                     prop = prop.multi, ps.method = ps.method,
                                     minPS = minPS, maxPS = maxPS)
    }

    ####### Construct the score in validation set ----------------------------------------------
    if (verbose >= 1) cat("  validating..\n")

    fit.score.valid <- scorecount(fit = fit,
                                  x.cate = x.cate.valid, time = time.valid,
                                  score.method = score.method)

    ####### Estimate the treatment effect in validation set --------------------------------------
    # if prop.onlyhigh ends with 1, estimate ATE once across all methods in the validation
    if (prop.onlyhigh[n.subgroup.onlyhigh] == 1) {
      est.prop1 <- exp(drcount(y = y.valid,
                               x.cate = x.cate.valid, x.ps = x.ps.valid, trt = trt.valid,
                               time = time.valid,
                               ps.method = ps.method, minPS = minPS, maxPS = maxPS, interactions = FALSE)$log.rate.ratio)
    }

    for (name in names(fit.score.valid)) {
      score.valid <- fit.score.valid[[name]]
      ate.valid.cv.i <-
        estcount.bilevel.subgroups(y = y.valid,
                                   x.cate = x.cate.valid, x.ps = x.ps.valid,
                                   time = time.valid, trt = trt.valid,
                                   score = score.valid, higher.y = higher.y,
                                   prop = prop.bi, onlyhigh = FALSE,
                                   ps.method = ps.method, minPS = minPS, maxPS = maxPS)
      result[[str_replace(name, "score", "ate")]]$ate.est.valid.high.cv[, cv.i] <-
        c(ate.valid.cv.i$ate.est.high, est.prop1)
      result[[str_replace(name, "score", "ate")]]$ate.est.valid.low.cv[, cv.i] <-
        ate.valid.cv.i$ate.est.low
      result[[str_replace(name, "score", "ate")]]$ate.est.valid.group.cv[, cv.i] <-
        estcount.multilevel.subgroup(y = y.valid,
                                     x.cate = x.cate.valid, x.ps = x.ps.valid,
                                     time = time.valid, trt = trt.valid,
                                     score = score.valid, higher.y = higher.y,
                                     prop = prop.multi,
                                     ps.method = ps.method, minPS = minPS, maxPS = maxPS)

      if (abc == TRUE) {
        temp.abc <- log(estcount.bilevel.subgroups(y = y.valid,
                                                   x.cate = x.cate.valid, x.ps = x.ps.valid,
                                                   time = time.valid, trt = trt.valid,
                                                   score = score.valid, higher.y = higher.y,
                                                   prop = prop.abc, onlyhigh = TRUE,
                                                   ps.method = ps.method, minPS = minPS, maxPS = maxPS))

        if (higher.y == TRUE) {
          # if higher scores are better, ABC requires subtracting the ATE from validation curve
          ref <- temp.abc - log(result$overall.ate.valid[cv.i])
        } else {
          #  if lower scores are better, ABC requires subtracting the validation curve from ATE
          ref <- log(result$overall.ate.valid[cv.i]) - temp.abc
        }
        result[[str_replace(name, "score", "ate")]]$abc.valid[cv.i] <-
          auc(x = prop.abc, y = ref, from = prop.abc[1], to = prop.abc[length(prop.abc)], type = "spline")
      }

    }

    if (verbose >= 1 & ("contrastReg" %in% score.method)) {
      cat("  convergence", converge.contrastReg.cv[cv.i], "\n")
      cat("  ", date(), "\n")
      setTxtProgressBar(pb, cv.i)
    }

  }

  if (verbose == 2) {
    if (any(is.na(unlist(result[str_replace(names(fit.score.train), "score", "ate")])))) {
      warning("Missing log rate ratio detected in the subgroups due to negative doubly robust estimator of
            the outcome for one or both treatment group(s).")
    }

    if ("contrastReg" %in% score.method) {
      result$ate.contrastReg <- append(result$ate.contrastReg,
                                       list(converge.contrastReg.cv = converge.contrastReg.cv))
      if (any(is.na(converge.contrastReg.cv))) {
        warning("Contrast regression algorithm did not converge in at least one cross-validation iteration.")
      }
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
  result$fgam <- fit$fgam
  result$higher.y <- higher.y
  result$abc <- abc
  result$cv.n <- cv.n
  result$response <- "count"
  result$formulas <- list(cate.model = cate.model, ps.model = ps.model, trt_labels = out$cat.trt)

  class(result) <- "precmed"

  return(result)
}


#' Estimation of the conditional average treatment effect (CATE) score for count data
#'
#' Provides singly robust and doubly robust estimation of CATE score with up to 5 scoring methods
#' among the following: Poisson regression, boosting, two regressions, contrast regression, and
#' negative binomial.
#'
#' @param cate.model A formula describing the outcome model to be fitted.
#' The outcome must appear on the left-hand side.
#' @param ps.model A formula describing the propensity score model to be fitted.
#' The treatment must appear on the left-hand side. The treatment must be a numeric vector
#' coded as 0/1. If data are from a randomized trial, specify \code{ps.model} as an intercept-only model.
#' @param data A data frame containing the variables in the outcome and propensity score model; a data frame with \code{n} rows
#' (1 row per observation).
#' @param score.method A vector of one or multiple methods to estimate the CATE score.
#' Allowed values are: \code{'boosting'}, \code{'poisson'}, \code{'twoReg'}, \code{'contrastReg'},
#' and \code{'negBin'}.
#' @param higher.y A logical value indicating whether higher (\code{TRUE}) or
#' lower (\code{FALSE}) values of the outcome are more desirable. Default is \code{TRUE}.
#' @param prop.cutoff A vector of numerical values (in (0, 1]) specifying percentiles of the
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
#' @param minPS A numerical value (in [0, 1]) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in (0, 1]) above which estimated propensity scores should be
#' truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.
#' @param initial.predictor.method A character vector for the method used to get initial
#' outcome predictions conditional on the covariates in \code{cate.model}. Only applies
#' when \code{score.method} includes \code{'twoReg'} or \code{'contrastReg'}. Allowed values
#' include one of \code{'poisson'} (fastest), \code{'boosting'} and \code{'gam'}.
#' Default is \code{'boosting'}.
#' @param xvar.smooth A vector of characters indicating the name of the variables used as
#' the smooth terms if \code{initial.predictor.method = 'gam'}. The variables must be selected
#' from the variables listed in \code{cate.model}.
#' Default is \code{NULL}, which uses all variables in \code{cate.model}.
#' @param tree.depth A positive integer specifying the depth of individual trees in boosting
#' (usually 2-3). Used only if \code{score.method = 'boosting'} or
#' if \code{initial.predictor.method = 'boosting'} with \code{score.method = 'twoReg'} or
#' \code{'contrastReg'}. Default is 2.
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
#' @return Returns a list containing the following components:
#' \itemize{
#'  \item{\code{ate.poisson}: }{A vector of numerical values of length \code{prop.cutoff}
#'  containing the estimated ATE in nested subgroups (defined by \code{prop.cutoff})
#'  constructed based on the estimated CATE scores with poisson regression.
#'  Only provided if \code{score.method} includes \code{'poisson'}.}
#'  \item{\code{ate.boosting}: }{Same as \code{ate.poisson}, but with the nested subgroups based
#'  the estimated CATE scores with boosting. Only provided if \code{score.method}
#'  includes \code{'boosting'}.}
#'  \item{\code{ate.twoReg}: }{Same as \code{ate.poisson}, but with the nested subgroups based
#'  the estimated CATE scores with two regressions.
#'  Only provided if \code{score.method} includes \code{'twoReg'}.}
#'  \item{\code{ate.contrastReg}: }{Same as \code{ate.poisson}, but with the nested subgroups based
#'  the estimated CATE scores with contrast regression.
#'  Only provided if \code{score.method} includes \code{'contrastReg'}.}
#'  \item{\code{ate.negBin}: }{Same as \code{ate.poisson}, but with the nested subgroups based
#'  the estimated CATE scores with negative binomial regression.
#'  Only provided if \code{score.method} includes \code{'negBin'}.}
#'  \item{\code{score.poisson}: }{A vector of numerical values of length n
#'  (number of observations in \code{data}) containing the estimated log-CATE scores
#'  according to the Poisson regression. Only provided if \code{score.method}
#'  includes \code{'poisson'}.}
#'  \item{\code{score.boosting}: }{Same as \code{score.poisson}, but with estimated log-CATE score
#'  according to boosting. Only provided if \code{score.method} includes
#'  \code{'boosting'}.}
#'  \item{\code{score.twoReg}: }{Same as \code{score.poisson}, but with estimated log-CATE score
#'  according to two regressions. Only provided if \code{score.method} includes
#'  \code{'twoReg'}.}
#'  \item{\code{score.contrastReg}: }{Same as \code{score.poisson}, but with estimated log-CATE score
#'  according to contrast regression. Only provided if \code{score.method} includes
#'  \code{'contrastReg'}.}
#'  \item{\code{score.negBin}: }{Same as \code{score.poisson}, but with estimated log-CATE score
#'  according to negative binomial regression. Only provided if \code{score.method}
#'  includes \code{'negBin'}.}
#'  \item{\code{fit}: }{Additional details on model fitting if \code{score.method}
#'  includes 'boosting' or 'contrastReg':}
#'  \itemize{
#'    \item{\code{result.boosting}: }{Details on the boosting model fitted to observations
#'    with treatment = 0 \code{($fit0.boosting)} and to observations with treatment = 1 \code{($fit1.boosting)}.
#'    Only provided if \code{score.method} includes \code{'boosting'}.}
#'    \item{\code{result.contrastReg$sigma.contrastReg}: }{Variance-covariance matrix of
#'    the estimated log-CATE coefficients in contrast regression.
#'    Only provided if \code{score.method} includes \code{'contrastReg'}.}
#'  }
#'  \item{\code{coefficients}: }{A data frame with the coefficients of the estimated log-CATE
#'  score by \code{score.method}. The data frame has number of rows equal to the number of
#'  covariates in \code{cate.model} and number of columns equal to \code{length(score.method)}.
#'  If \code{score.method} includes \code{'contrastReg'}, the data frame has an additional
#'  column containing the standard errors of the coefficients estimated with contrast regression.
#'  \code{'boosting'} does not have coefficient results because tree-based methods typically do not
#'  express the log-CATE as a linear combination of coefficients and covariates.}
#' }
#'
#' @details The CATE score represents an individual-level treatment effect, estimated with
#' either Poisson regression, boosting or negative binomial regression applied separately by
#' treatment group or with two doubly robust estimators, two regressions and contrast regression
#' (Yadlowsky, 2020) applied to the entire dataset.
#'
#' \code{\link{catefitcount}()} provides the coefficients of the CATE score for each scoring method requested
#' through \code{score.method}. Currently, contrast regression is the only method which allows
#' for inference of the CATE coefficients by providing standard errors of the coefficients.
#' The coefficients can be used to learn the effect size of each variable and predict the
#' CATE score for a new observation.
#'
#' \code{\link{catefitcount}()} also provides the predicted CATE score of each observation in the data set,
#' for each scoring method. The predictions allow ranking the observations from potentially
#' high responders to the treatment to potentially low or standard responders.
#'
#' The estimated ATE among nested subgroups of high responders are also provided by scoring method.
#' Note that the ATEs in \code{\link{catefitcount}()} are derived based on the CATE score which is estimated
#' using the same data sample. Therefore, overfitting may be an issue. \code{\link{catecvcount}()} is more
#' suitable to inspect the estimated ATEs across scoring methods as it implements internal cross
#' validation to reduce optimism.
#'
#' @references Yadlowsky, S., Pellegrini, F., Lionetto, F., Braune, S., & Tian, L. (2020).
#' \emph{Estimation and validation of ratio-based conditional average treatment effects using
#' observational data. Journal of the American Statistical Association, 1-18.}
#' \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1772080}
#'
#' @seealso \code{\link{catecvcount}()}
#'
#' @examples
#' catefit <- catefitcount(cate.model = y ~ age + female + previous_treatment +
#'                                previous_cost + previous_number_relapses + offset(log(years)),
#'               ps.model = trt ~ age + previous_treatment,
#'               data = countExample,
#'               higher.y = FALSE,
#'               score.method = "poisson",
#'               seed = 999)
#' @export

catefitcount <- function(data,
                         score.method,
                         cate.model,
                         ps.model,
                         ps.method = "glm",
                         initial.predictor.method = "boosting",
                         minPS = 0.01,
                         maxPS = 0.99,
                         verbose = 0,
                         higher.y = TRUE,
                         prop.cutoff = seq(0.5, 1, length = 6),
                         xvar.smooth = NULL,
                         tree.depth = 2,
                         n.trees.boosting = 200,
                         B = 3,
                         Kfold = 6,
                         error.maxNR = 1e-3,
                         max.iterNR = 150,
                         tune = c(0.5, 2),
                         seed = NULL,
                         plot.gbmperf = FALSE,
                         ...) {


  # Set seed once for reproducibility
  set.seed(seed)

  if (verbose >= 1) t.start <- Sys.time()

  #### CHECK ARGUMENTS ####
  arg.checks(
    fun = "pm", response = "count", data = data, higher.y = higher.y, score.method = score.method, prop.cutoff = prop.cutoff,
    ps.method = ps.method, minPS = minPS, maxPS = maxPS,
    initial.predictor.method = initial.predictor.method,
    tree.depth = tree.depth, n.trees.boosting = n.trees.boosting, B = B, Kfold = Kfold, plot.gbmperf = plot.gbmperf,
    error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune
  )

  #### PRE-PROCESSING ####
  out <- data.preproc(fun = "pm", cate.model = cate.model, ps.model = ps.model,
                      data = data, prop.cutoff = prop.cutoff, ps.method = ps.method)
  y <- out$y
  trt <- out$trt
  x.ps <- out$x.ps
  x.cate <- out$x.cate
  time <- out$time
  prop <- out$prop

  #### FUNCTION STARTS HERE ####
  result <- vector("list", length(score.method) * 2 + 2)
  names(result) <- c(paste0("ate.", score.method),
                     paste0("score.", score.method), "fit", "coefficients")

  fit <- intxcount(y = y, trt = trt, x.cate = x.cate, x.ps = x.ps, time = time,
                   score.method = score.method,
                   ps.method = ps.method, minPS = minPS, maxPS = maxPS,
                   initial.predictor.method = initial.predictor.method,
                   xvar.smooth = xvar.smooth,
                   tree.depth = tree.depth, n.trees.boosting = n.trees.boosting,
                   Kfold = Kfold, B = B, plot.gbmperf = plot.gbmperf,
                   error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune, ...)

  if (fit$best.iter == n.trees.boosting) {
    warning(paste("The best boosting iteration was iteration number", n.trees.boosting, " out of ", n.trees.boosting, ". Consider increasing the maximum number of trees and turning on boosting performance plot (plot.gbmperf = TRUE).", sep = ""))
  }

  # Check NA in coefficients of the score
  if ("poisson" %in% score.method & sum(is.na(fit$result.poisson)) > 0) {
    warning("One or more coefficients in the score (Poisson) are NA.
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
  if ("negBin" %in% score.method & sum(is.na(fit$result.negBin)) > 0) {
    warning("One or more coefficients in the score (negative binomial) are NA.
            Consider inspecting the distribution of the covariates in cate.model.")
  }

  fit.score <- scorecount(fit = fit,
                          x.cate = x.cate, time = time,
                          score.method = score.method)

  for (name in names(fit.score)) {
    score <- fit.score[[name]]
    result[[name]] <- score
    ate <- estcount.bilevel.subgroups(y = y,
                                      x.cate = x.cate, x.ps = x.ps,
                                      time = time, trt = trt,
                                      score = score, higher.y = higher.y,
                                      prop = prop, onlyhigh = TRUE,
                                      ps.method = ps.method, minPS = minPS, maxPS = maxPS)
    result[[str_replace(name, "score", "ate")]] <- ate
    names(result[[str_replace(name, "score", "ate")]]) <- paste0("prop", round(prop, 2))
  }

  if (sum(score.method %in% c("poisson", "twoReg", "contrastReg", "negBin")) > 0) {
    cf <- data.frame(matrix(NA, nrow = ncol(x.cate) + 1, ncol = 5))
    colnames(cf) <- c("poisson", "twoReg", "contrastReg", "SE_contrastReg", "negBin")
    rownames(cf) <- c("(Intercept)", colnames(x.cate))

    if ("poisson" %in% score.method) cf$poisson <- fit$result.poisson

    if ("twoReg" %in% score.method) cf$twoReg <- fit$result.twoReg

    if ("contrastReg" %in% score.method) {
      cf$contrastReg <- fit$result.contrastReg$delta.contrastReg
      cf$SE_contrastReg <- sqrt(diag(fit$result.contrastReg$sigma.contrastReg))
    }

    if ("negBin" %in% score.method) cf$negBin <- fit$result.negBin

    result$coefficients <- cf[, colSums(is.na(cf)) != nrow(cf), drop = FALSE]
  }


  if (any(is.na(unlist(result[str_replace(names(fit.score), "score", "ate")])))) {
    warning("Missing log rate ratio detected due to negative doubly robust estimator of y|x
            for one or both treatment group(s).")
  }

  if ("boosting" %in% score.method) result$fit$result.boosting <-
    fit$result.boosting
  if ("contrastReg" %in% score.method) result$fit$result.contrastReg$sigma.contrastReg <-
    fit$result.contrastReg$sigma.contrastReg

  if (verbose >= 1) {
    t.end <- Sys.time()
    t.diff <- round(difftime(t.end, t.start),2)
    cat('Total runtime :',as.numeric(t.diff), attributes(t.diff)$units, '\n')
  }

  class(result) <- "precmed"
  return(result)
}


#' Doubly robust estimator of and inference for the average treatment effect for count data
#'
#' Doubly robust estimator of the average treatment effect between two treatments, which is the rate ratio
#' for count outcomes. Bootstrap is used for inference.
#'
#' @param cate.model A formula describing the outcome model to be fitted.
#' The outcome must appear on the left-hand side.
#' @param ps.model A formula describing the propensity score (PS) model to be fitted. The treatment must
#' appear on the left-hand side. The treatment must be a numeric vector coded as 0/1.
#' If data are from a randomized controlled trial, specify \code{ps.model = ~1} as an intercept-only model.
#' @param data A data frame containing the variables in the outcome, propensity score, and inverse
#' probability of censoring models (if specified); a data frame with \code{n} rows (1 row per observation).
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
#' @param interactions A logical value indicating whether the outcome model should assume interactions
#' between \code{x} and \code{trt}. If \code{TRUE}, interactions will be assumed only if at least 10 patients
#' received each treatment option. Default is \code{TRUE}.
#' @param n.boot A numeric value indicating the number of bootstrap samples used. Default is \code{500}.
#' @param seed An optional integer specifying an initial randomization seed for reproducibility.
#' Default is \code{NULL}, corresponding to no seed.
#' @param verbose An integer value indicating whether intermediate progress messages should
#' be printed. \code{1} indicates messages are printed and \code{0} otherwise. Default is \code{0}.
#'
#' @return Return an item of the class \code{atefit} with the following elements:
#' \itemize{
#'   \item{\code{log.rate.ratio}: } A vector of numeric values of the estimated log rate ratio of \code{trt=1}
#'   over \code{trt=0}, bootstrap standard error, lower and upper limits of 95% confidence interval, and the p-value.
#'   \item{\code{rate0}: } A numeric value of the estimated rate in the group \code{trt=0}.
#'   \item{\code{rate1}: } A numeric value of the estimated rate in the group \code{trt=1}.
#'   \item{\code{trt.boot}: } Estimated log rate ratios in each bootstrap sample.
#'   \item{\code{warning}: } A warning message produced if the treatment variable was not coded as 0/1. The key
#'   to map the original coding of the variable to a 0/1 key is displayed in the warning to facilitate the
#'   interpretation of the remaining of the output.
#' }
#'
#' @details This helper function estimates the average treatment effect (ATE) between two treatment groups in a given
#' dataset. The ATE is estimated with a doubly robust estimator that accounts for imbalances in covariate distributions
#' between the two treatment groups with inverse probability treatment weighting. For count outcomes, the estimated ATE
#' is the estimated rate ratio between treatment 1 versus treatment 0. The log-transformed ATEs are returned, as well
#' as the rate in either treatment group. The variability of the estimated rate ratio is also calculated using bootstrap.
#' Additional outputs include standard error of the log rate ratio, 95% confidence interval, p-value, and a histogram of
#' the bootstrap estimates.
#'
#' @examples
#' \dontrun{
#' output <- atefitcount(cate.model = y ~ age + female + previous_treatment +
#'                                previous_cost + previous_number_relapses + offset(log(years)),
#'                             ps.model = trt ~ age + previous_treatment,
#'                             data = countExample,
#'                             seed = 999, verbose = 1)
#' output
#' plot(output)
#'}
#' @export
#'
#' @importFrom dplyr mutate
#'

atefitcount <- function(data,
                        cate.model,
                        ps.model,
                        ps.method = "glm",
                        minPS = 0.01,
                        maxPS = 0.99,
                        verbose = 0,
                        interactions = TRUE,
                        n.boot = 500,
                        seed = NULL) {

  # Set seed once for reproducibility
  set.seed(seed)

  #### CHECK ARGUMENTS ####
  arg.checks(fun = "drinf",
             response = "count",
             data = data,
             ps.method = ps.method,
             minPS = minPS,
             maxPS = maxPS,
             interactions = interactions,
             n.boot = n.boot,
             plot.boot = FALSE)

  #### PRE-PROCESSING ####
  preproc <- data.preproc(fun = "drinf", cate.model = cate.model, ps.model = ps.model,
                          data = data, ps.method = ps.method)
  y <- preproc$y
  trt <- preproc$trt
  x.ps <- preproc$x.ps
  x.cate <- preproc$x.cate
  time <- preproc$time

  # Point estimate
  est <- drcount(y = y,
                 x.cate = x.cate,
                 x.ps = x.ps,
                 trt = trt,
                 time = time,
                 ps.method = ps.method,
                 minPS = minPS,
                 maxPS = maxPS,
                 interactions = interactions)
  logrr <- est$log.rate.ratio

  # Apply bootstrap
  n <- length(y)
  if (is.na(logrr)) {
    stop("Impossible to use bootstrap, log(RR)=NA.")
  } else {
    trt.boot <- rep(0, n.boot)

    pb   <- txtProgressBar(min = 1,
                           max = n.boot,
                           style = 3)

    for (i in seq(n.boot)) {
      idsub.boot <- sample(n, size = n, replace = TRUE)
      trt.boot[i] <- drcount(y = y[idsub.boot],
                             x.cate = x.cate[idsub.boot, , drop = FALSE],
                             x.ps = x.ps[idsub.boot, , drop = FALSE],
                             trt = trt[idsub.boot],
                             time = time[idsub.boot],
                             ps.method = ps.method,
                             minPS = minPS,
                             maxPS = maxPS,
                             interactions = interactions)$log.rate.ratio

      if (verbose == 1) setTxtProgressBar(pb, i)
    }
    close(pb)

    out <- c()
    out$response <- "count"
    se.est <- sd(trt.boot, na.rm = TRUE)
    out$log.rate.ratio <- data.frame(estimate = logrr,
                                     SE = se.est,
                                     CI.lower = logrr - qnorm(0.975) * se.est,
                                     CI.upper = logrr + qnorm(0.975) * se.est,
                                     pvalue = 1 - pchisq(logrr^2 / se.est^2, 1))
    rownames(out$log.rate.ratio) <- "log.rate.ratio"
    out$rate0 <- data.frame(estimate = est$rate0)
    rownames(out$rate0) <- "rate0"
    out$rate1 <- data.frame(estimate = est$rate1)
    rownames(out$rate1) <- "rate1"
    out$n.boot <- n.boot # Number of  bootstrap samples
    out$trt.boot <- trt.boot #bootstrap estimates of the treatment effect
    out$warning <- preproc$warning


  }

  class(out) <- "atefit"

  return(out)
}
