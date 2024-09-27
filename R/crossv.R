#' Cross-validation of the conditional average treatment effect (CATE) score for count, survival or continuous outcomes
#'
#' Provides (doubly robust) estimation of the average treatment effect (ATE) for count, survival or continuous
#' outcomes in nested and mutually exclusive subgroups of patients defined by an estimated conditional
#' average treatment effect (CATE) score via cross-validation (CV).
#'
#' @param response A string describing the type of outcome in the data.
#' Allowed values include "count" (see \code{\link{catecvcount}()}), "survival"
#' (see \code{\link{catecvsurv}()}) and "continuous" (see
#' \code{\link{catecvmean}()}) .
#' @param data A data frame containing the variables in the outcome, propensity score, and inverse
#' probability of censoring models (if specified); a data frame with \code{n} rows (1 row per observation).
#' @param score.method A vector of one or multiple methods to estimate the CATE score.
#' Allowed values are: \code{'boosting'}, \code{'twoReg'}, \code{'contrastReg'}, \code{'poisson'} (count and survival outcomes only),
#' \code{'randomForest'} (survival, continuous outcomes only), \code{negBin} (count outcomes only), \code{'gam'} (continuous outcomes only),
#' \code{'gaussian'} (continuous outcomes only).
#' @param cate.model A formula describing the outcome model to be fitted.
#' The outcome must appear on the left-hand side. For survival outcomes, a
#' \code{Surv} object must be used to describe the outcome.
#' @param ps.model A formula describing the propensity score (PS) model to be fitted. The treatment must
#' appear on the left-hand side. The treatment must be a numeric vector coded as 0/1.
#' If data are from a randomized controlled trial, specify \code{ps.model = ~1} as an intercept-only model.
#' @param ps.method A character value for the method to estimate the propensity score.
#' Allowed values include one of:
#' \code{'glm'} for logistic regression with main effects only (default), or
#' \code{'lasso'} for a logistic regression with main effects and LASSO penalization on
#' two-way interactions (added to the model if interactions are not specified in \code{ps.model}).
#' Relevant only when \code{ps.model} has more than one variable.
#' @param init.model A formula describing the initial predictor model. The outcome must appear on the left-hand side.
#' It must be specified when \code{score.method = contrastReg} or \code{twoReg}.
#' @param initial.predictor.method A character vector for the method used to get initial
#' outcome predictions conditional on the covariates specified in \code{cate.model}. Only applies
#' when \code{score.method} includes \code{'twoReg'} or \code{'contrastReg'}. Allowed values include one of
#' \code{'randomForest'} (survival outcomes only), \code{'boosting'}, \code{'logistic'}
#' (survival outcomes only, fast), \code{'poisson'} (count outcomes only, fast), \code{'gaussian'} (continuous outcomes only) and
#' \code{'gam'} (count and continuous outcomes only). Default is \code{NULL}, which assigns \code{'boosting'}
#' for count outcomes and \code{'randomForest'} for survival outcomes.
#' @param ipcw.model A formula describing the inverse probability of censoring weighting (IPCW)
#' model to be fitted. The left-hand side must be empty. Only applies for survival outcomes.
#' Default is \code{NULL}, which corresponds to specifying the IPCW with the same covariates
#' as the outcome model \code{cate.model}, plus the treatment.
#' @param ipcw.method A character value for the censoring model. Only applies for survival
#' outcomes. Allowed values are: \code{'breslow'} (Cox regression with Breslow estimator of t
#' he baseline survivor function), \code{'aft (exponential)'}, \code{'aft (weibull)'},
#' \code{'aft (lognormal)'} or \code{'aft (loglogistic)'} (accelerated failure time model
#' with different distributions for y variable). Default is \code{'breslow'}.
#' @param minPS A numerical value (in `[0, 1]`) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in `(0, 1]`) above which estimated propensity scores should be
#' truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.
#' @param followup.time A column name in \code{data} specifying the maximum follow-up time,
#' interpreted as the potential censoring time. Only applies for survival outcomes.
#' Default is \code{NULL}, which corresponds to unknown potential censoring time.
#' @param tau0 The truncation time for defining restricted mean time lost. Only applies for
#' survival outcomes. Default is \code{NULL}, which corresponds to setting the truncation time as the
#' maximum survival time in the data.
#' @param higher.y A logical value indicating whether higher (\code{TRUE}) or
#' lower (\code{FALSE}) values of the outcome are more desirable. Default is \code{TRUE}.
#' @param prop.cutoff A vector of numerical values (in `(0, 1]`) specifying percentiles of the
#' estimated log CATE scores to define nested subgroups. Each element represents the cutoff to
#' separate observations in nested subgroups (below vs above cutoff).
#' The length of \code{prop.cutoff} is the number of nested subgroups.
#' An equally-spaced sequence of proportions ending with 1 is recommended.
#' Default is \code{seq(0.5, 1, length = 6)}.
#' @param prop.multi A vector of numerical values (in `[0, 1]`) specifying percentiles of the
#' estimated log CATE scores to define mutually exclusive subgroups.
#' It should start with 0, end with 1, and be of \code{length(prop.multi) > 2}.
#' Each element represents the cutoff to separate the observations into
#' \code{length(prop.multi) - 1} mutually exclusive subgroups.
#' Default is \code{c(0, 1/3, 2/3, 1)}.
#' @param abc A logical value indicating whether the area between curves (ABC) should be calculated
#' at each cross-validation iterations, for each \code{score.method}. Default is \code{TRUE}.
#' @param train.prop A numerical value (in `(0, 1)`) indicating the proportion of total data used
#' for training. Default is \code{3/4}.
#' @param cv.n A positive integer value indicating the number of cross-validation iterations.
#' Default is \code{10}.
#' @param error.max A numerical value > 0 indicating the tolerance (maximum value of error)
#' for the largest standardized absolute difference in the covariate distributions or in the
#' doubly robust estimated rate ratios between the training and validation sets. This is used
#' to define a balanced training-validation splitting. Default is \code{0.1}.
#' @param max.iter A positive integer value indicating the maximum number of iterations when
#' searching for a balanced training-validation split. Default is \code{5,000}.
#' @param surv.min Lower truncation limit for the probability of being censored.
#' It must be a positive value and should be chosen close to 0. Only applies for survival
#' outcomes. Default is \code{0.025}.
#' @param xvar.smooth.score A vector of characters indicating the name of the variables used as
#' the smooth terms if \code{score.method = 'gam'}. The variables must be selected
#' from the variables listed in \code{cate.model}.
#' @param xvar.smooth.init A vector of characters indicating the name of the variables used as
#' the smooth terms if \code{initial.predictor.method = 'gam'}. The variables must be selected
#' from the variables listed in \code{init.model}.
#' Default is \code{NULL}, which uses all variables in \code{init.model}.
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
#' @param plot.gbmperf A logical value indicating whether to plot the performance measures in
#' boosting. Used only if \code{score.method = 'boosting'} or if \code{score.method = 'twoReg'}
#' or \code{'contrastReg'} and \code{initial.predictor.method = 'boosting'}. Default is \code{TRUE}.
#' @param verbose An integer value indicating what kind of intermediate progress messages should
#' be printed. \code{0} means no outputs. \code{1} means only progress bar and run time.
#' \code{2} means progress bar, run time, and all errors and warnings. Default is \code{0}.
#'
#' @return For count response, see description of outputs in \code{\link{catecvcount}()}.
#' For survival response, see description of outputs in \code{\link{catecvsurv}()}.
#' For continuous response, see description of outputs in \code{\link{catecvmean}()}.
#'
#' @details For count response, see details in \code{\link{catecvcount}()}.
#' For survival response, see details in \code{\link{catecvsurv}()}.
#' For continuous response, see details in \code{\link{catecvmean}()}.
#'
#' @references Yadlowsky, S., Pellegrini, F., Lionetto, F., Braune, S., & Tian, L. (2020).
#' \emph{Estimation and validation of ratio-based conditional average treatment effects using
#' observational data. Journal of the American Statistical Association, 1-18.}
#' \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1772080}
#'
#' @seealso \code{\link{catefit}()} function and \code{\link{boxplot}()}, \code{\link{abc}} methods for
#' \code{"precmed"} objects.
#'
#' @examples
#' \donttest{
#' cate_1 <- catecv(response = "count",
#'                  data = countExample,
#'                  score.method = "poisson",
#'                  cate.model = y ~ age + female + previous_treatment +
#'                               previous_cost + previous_number_relapses +
#'                               offset(log(years)),
#'                  ps.model = trt ~ age + previous_treatment,
#'                  higher.y = FALSE, cv.n = 5, seed = 999, verbose = 1)
#'
#' plot(cate_1, ylab = "RMTL ratio of drug1 vs drug0 in each subgroup")
#' boxplot(cate_1, ylab = "RMTL ratio of drug1 vs drug0 in each subgroup")
#' abc(cate_1)
#'
#' # Survival outcome
#' library(survival)
#' tau0 <- with(survivalExample,
#'              min(quantile(y[trt == "drug1"], 0.95), quantile(y[trt == "drug0"], 0.95)))
#'
#' cate_2 <- catecv(response = "survival",
#'                  data = survivalExample,
#'                  score.method = c("poisson", "randomForest"),
#'                  cate.model = Surv(y, d) ~ age + female + previous_cost +
#'                               previous_number_relapses,
#'                  ps.model = trt ~ age + previous_treatment,
#'                  initial.predictor.method = "randomForest",
#'                  ipcw.model = ~ age + previous_cost + previous_treatment,
#'                  tau0 = tau0,
#'                  higher.y = TRUE,
#'                  surv.min = 0.025,
#'                  cv.n = 5,
#'                  seed = 999)
#'
#' plot(cate_2, ylab = "RMTL ratio of drug1 vs drug0 in each subgroup")
#' boxplot(cate_2, ylab = "RMTL ratio of drug1 vs drug0 in each subgroup")
#' abc(cate_2)
#'
#'
#' }
#'
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom gam gam
#' @importFrom gbm gbm gbm.perf predict.gbm
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom graphics hist lines
#' @importFrom MASS glm.nb
#' @importFrom randomForestSRC rfsrc predict.rfsrc
#' @importFrom stats as.formula coef glm median model.frame model.matrix model.offset model.response na.omit optim pchisq predict qnorm quantile sd var
#' @importFrom stringr str_replace str_extract str_detect
#' @importFrom survival Surv coxph coxph.detail survreg
#' @importFrom utils setTxtProgressBar txtProgressBar

catecv <- function(response,
                   data,
                   score.method,  # Mandatory arguments (count & survival & continuous)
                   cate.model,
                   ps.model,
                   ps.method = "glm",
                   init.model = NULL, # This is used in continuous data, when contrast/two regression is used
                   initial.predictor.method = NULL,
                   ipcw.model = NULL,
                   ipcw.method = "breslow",
                   minPS = 0.01,
                   maxPS = 0.99,
                   followup.time = NULL,
                   tau0 = NULL, # Non-mandatory arguments survival only
                   higher.y = TRUE,
                   prop.cutoff = seq(0.5, 1, length = 6),
                   prop.multi = c(0, 1/3, 2/3, 1),
                   abc = TRUE, # Non-mandatory arguments survival & count & continuous (except xvar.smooth)
                   train.prop = 3/4,
                   cv.n = 10,
                   error.max = 0.1,
                   max.iter = 5000,
                   surv.min = 0.025,
                   xvar.smooth.score = NULL,
                   xvar.smooth.init = NULL,
                   tree.depth = 2,
                   n.trees.rf = 1000,
                   n.trees.boosting = 200,
                   B = 3,
                   Kfold = 5,
                   error.maxNR = 1e-3,
                   max.iterNR = 150,
                   tune = c(0.5, 2),
                   seed = NULL,
                   plot.gbmperf = TRUE,
                   verbose = 0) {

  stopifnot("`response` must be either `count`, `survival`, or `continuous`." = any(response == c("count", "survival", "continuous")))

  .args <- as.list(match.call())[-1]
  .args$response <- NULL

  switch(response,
         count = {
           cvout <- do.call(catecvcount, .args)
         },
         survival = {
           cvout <- do.call(catecvsurv, .args)
         },
         continuous = {
           cvout <- do.call(catecvmean, .args)
         }
  )
  return(cvout)
}

#' Cross-validation of the conditional average treatment effect (CATE) score for survival outcomes
#'
#' Provides doubly robust estimation of the average treatment effect (ATE) by the
#' RMTL (restricted mean time lost) ratio in nested and mutually exclusive subgroups of patients
#' defined by an estimated conditional average treatment effect (CATE) score via
#' cross-validation (CV). The CATE score can be estimated with up to 5 methods among the following:
#' Random forest, boosting, poisson regression, two regressions, and contrast regression
#' (see \code{score.method}).
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
#' @param abc A logical value indicating whether the area between curves (ABC) should be calculated
#' at each cross-validation iterations, for each \code{score.method}. Default is \code{TRUE}.
#' @param prop.cutoff A vector of numerical values (in `(0, 1]`) specifying percentiles of the
#' estimated log CATE scores to define nested subgroups. Each element represents the cutoff to
#' separate observations in nested subgroups (below vs above cutoff).
#' The length of \code{prop.cutoff} is the number of nested subgroups.
#' An equally-spaced sequence of proportions ending with 1 is recommended.
#' Default is \code{seq(0.5, 1, length = 6)}.
#' @param prop.multi A vector of numerical values (in `[0, 1]`) specifying percentiles of the
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
#' @param minPS A numerical value (in `[0, 1]`) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in `(0, 1]`) above which estimated propensity scores should be
#' truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.
#' @param train.prop A numerical value (in `(0, 1)`) indicating the proportion of total data used
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
#' @param plot.gbmperf A logical value indicating whether to plot the performance measures in
#' boosting. Used only if \code{score.method = 'boosting'} or if \code{score.method = 'twoReg'}
#' or \code{'contrastReg'} and \code{initial.predictor.method = 'boosting'}. Default is \code{TRUE}.
#' @param verbose An integer value indicating what kind of intermediate progress messages should
#' be printed. \code{0} means no outputs. \code{1} means only progress bar and run time.
#' \code{2} means progress bar, run time, and all errors and warnings. Default is \code{0}.
#'
#' @return Returns a list containing the following components saved as a \code{"precmed"} object:
#' \itemize{
#'  \item{\code{ate.randomForest}: }{A list of ATE output measured by the RMTL ratio if
#'  \code{score.method} includes \code{'randomForest'}:}
#'  \itemize{
#'     \item{\code{ate.est.train.high.cv}: }{A matrix of numerical values with
#'     \code{length(prop.cutoff)} rows and \code{cv.n} columns.
#'     The ith column/jth row cell contains the estimated ATE in the nested subgroup of high responders
#'     defined by CATE score above (if \code{higher.y = FALSE}) or below (if \code{higher.y = TRUE}) the
#'     \code{prop.cutoff[j]}x100\% percentile of the estimated CATE score in the training set in the ith
#'     cross-validation iteration.}
#'     \item{\code{ate.est.train.low.cv}: }{A matrix of numerical values with
#'     \code{length(prop.cutoff) - 1} rows and \code{cv.n} columns.
#'     TThe ith column/jth row cell contains the estimated ATE in the nested subgroup of low responders
#'     defined by CATE score below (if \code{higher.y = FALSE}) or above (if \code{higher.y = TRUE}) the
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
#'  \item{\code{ate.boosting}: }{A list of results similar to \code{ate.randomForest} output
#'  if \code{score.method} includes \code{'boosting'}.}
#'  \item{\code{ate.poisson}: }{A list of results similar to \code{ate.randomForest} output
#'  if \code{score.method} includes \code{'poisson'}.}
#'  \item{\code{ate.twoReg}: }{A list of results similar to \code{ate.randomForest} output
#'  if \code{score.method} includes \code{'twoReg'}.}
#'  \item{\code{ate.contrastReg}: }{A list of results similar to \code{ate.randomForest} output
#'  if \code{score.method} includes \code{'contrastReg'}.
#'  This method has an additional element in the list of results:}
#'  \itemize{
#'     \item{\code{converge.contrastReg.cv}: }{A vector of logical value of length \code{cv.n}.
#'     The ith element indicates whether the algorithm converged in the ith cross-validation
#'     iteration.}
#'  }
#'  \item{\code{hr.randomForest}: }{A list of adjusted hazard ratio if \code{score.method} includes
#'  \code{'randomForest'}:}
#'  \itemize{
#'     \item{\code{hr.est.train.high.cv}: }{A matrix of numerical values with
#'     \code{length(prop.cutoff)} rows and \code{cv.n} columns.
#'     The ith column/jth row cell contains the estimated HR in the nested subgroup of high responders
#'     defined by CATE score above (if \code{higher.y = FALSE}) or below (if \code{higher.y = TRUE}) the
#'     \code{prop.cutoff[j]}x100\% percentile of the estimated CATE score in the training set in the ith
#'     cross-validation iteration.}
#'     \item{\code{hr.est.train.low.cv}: }{A matrix of numerical values with
#'     \code{length(prop.cutoff) - 1} rows and \code{cv.n} columns.
#'     TThe ith column/jth row cell contains the estimated HR in the nested subgroup of low responders
#'     defined by CATE score below (if \code{higher.y = FALSE}) or above (if \code{higher.y = TRUE}) the
#'     \code{prop.cutoff[j]}x100\% percentile of the estimated CATE score in the training set in the ith
#'     cross-validation iteration.}
#'     \item{\code{hr.est.valid.high.cv}: }{Same as \code{hr.est.train.high.cv},
#'     but in the validation set.}
#'     \item{\code{hr.est.valid.low.cv}: }{Same as \code{hr.est.train.low.cv},
#'     but in the validation set.}
#'     \item{\code{hr.est.train.group.cv}: }{A matrix of numerical values with
#'     \code{length(prop.multi) - 1} rows and \code{cv.n} columns.
#'     The ith column contains the estimated HR in \code{length(prop.multi) - 1}
#'     mutually exclusive subgroups defined by \code{prop.multi} in the training set in ith
#'     cross-validation iteration.}
#'     \item{\code{hr.est.valid.group.cv}: }{Same as \code{hr.est.train.group.cv}, but in the
#'     validation set.}
#'  }
#'  \item{\code{hr.boosting}: }{A list of results similar to \code{hr.randomForest} output
#'  if \code{score.method} includes \code{'boosting'}.}
#'  \item{\code{hr.poisson}: }{A list of results similar to \code{hr.randomForest} output
#'  if \code{score.method} includes \code{'poisson'}.}
#'  \item{\code{hr.twoReg}: }{A list of results similar to \code{hr.randomForest} output
#'  if \code{score.method} includes \code{'twoReg'}.}
#'  \item{\code{hr.contrastReg}: }{A list of results similar to \code{hr.randomForest} output
#'  if \code{score.method} includes \code{'contrastReg'}.
#'  \item{\code{props}: }{A list of 3 elements:}
#'  \itemize{
#'     \item{\code{prop.onlyhigh}: }{The original argument \code{prop.cutoff},
#'     reformatted as necessary.}
#'     \item{\code{prop.bi}: }{The original argument \code{prop.cutoff},
#'     similar to \code{prop.onlyhigh} but reformatted to exclude 1.}
#'     \item{\code{prop.multi}: }{The original argument \code{prop.multi},
#'     reformatted as necessary to include 0 and 1.}
#'  }
#'  \item{\code{overall.ate.train}: }{A vector of numerical values of length \code{cv.n}.
#'  The ith element contains the ATE (RMTL ratio) in the training set of the ith cross-validation
#'  iteration, estimated with the doubly robust estimator.}
#'  \item{\code{overall.hr.train}: }{A vector of numerical values of length \code{cv.n}.
#'  The ith element contains the ATE (HR) in the training set of the ith cross-validation
#'  iteration.}
#'  \item{\code{overall.ate.valid}: }{A vector of numerical values of length \code{cv.n}.
#'  The ith element contains the ATE (RMTL ratio) in the validation set of the ith cross-validation
#'  iteration, estimated with the doubly robust estimator.}
#' \item{\code{overall.hr.valid}: }{A vector of numerical values of length \code{cv.n}.
#' The ith element contains the ATE (HR) in the validation set of the ith cross-validation
#' iteration.}
#'  \item{\code{errors/warnings}: }{A nested list of errors and warnings that were wrapped during the
#'  calculation of ATE. Errors and warnings are organized by \code{score.method} and
#'  position in the CV flow.}
#'  \item{\code{higher.y}: }{The original \code{higher.y} argument.}
#'  \item{\code{abc}: }{The original \code{abc} argument.}
#'  \item{\code{cv.n}: }{The original \code{cv.n} argument.}
#'  \item{\code{response}: }{The type of response. Always 'survival' for this function.}
#'  \item{\code{formulas}:}{A list of 3 elements: (1) \code{cate.model} argument,
#'  (2) \code{ps.model} argument and (3) original labels of the left-hand side variable in
#'  \code{ps.model} (treatment) if it was not 0/1.}
#' }
#' }
#' @details The CATE score represents an individual-level treatment effect expressed as the
#' restricted mean survival time (RMTL) ratio) for survival outcomes. It can be estimated with boosting,
#' Poisson regression, random forest, and the doubly robust estimator two regressions (Yadlowsky, 2020)
#' applied separately by treatment group or with the other doubly robust estimator contrast regression
#' (Yadlowsky, 2020) applied to the entire data set.
#'
#' Internal CV is applied to reduce optimism in choosing the CATE estimation method that
#' captures the most treatment effect heterogeneity. The CV is applied by repeating the
#' following steps \code{cv.n} times:
#'
#' \enumerate{
#'  \item Split the data into a training and validation set according to \code{train.prop}.
#'  The training and validation sets must be balanced with respect to covariate distributions
#'  and doubly robust RMTL ratio estimates (see \code{error.max}).
#'  \item Estimate the CATE score in the training set with the specified scoring method.
#'  \item Predict the CATE score in the validation set using the scoring model fitted from
#'  the training set.
#'  \item Build nested subgroups of treatment responders in the training and validation sets,
#'  separately, and estimate the ATE within each nested subgroup. For each element i of
#'  \code{prop.cutoff} (e.g., \code{prop.cutoff[i]} = 0.6), take the following steps:
#'  \enumerate{
#'    \item Identify high responders as observations with the 60\%
#'    (i.e., \code{prop.cutoff[i]}x100\%) highest (if \code{higher.y = FALSE}) or
#'    lowest (if \code{higher.y = TRUE}) estimated CATE scores.
#'    \item Estimate the ATE in the subgroup of high responders using a doubly robust estimator.
#'    \item Conversely, identify low responders as observations with the 40\%
#'    (i.e., 1 - \code{prop.cutoff[i]}x100\%) lowest (if \code{higher.y} = FALSE) or
#'    highest (if \code{higher.y} = TRUE) estimated CATE scores.
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
#' @seealso \code{\link{catefitsurv}()} function and \code{\link{boxplot}()}, \code{\link{abc}} methods for
#' \code{"precmed"} objects.
#'
#' @examples
#' \donttest{
#' library(survival)
#'
#' tau0 <- with(survivalExample,
#'              min(quantile(y[trt == "drug1"], 0.95), quantile(y[trt == "drug0"], 0.95)))
#'
#' catecv <- catecvsurv(data = survivalExample,
#'                      score.method = "poisson",
#'                      cate.model = Surv(y, d) ~ age + female + previous_cost +
#'                                                previous_number_relapses,
#'                      ps.model = trt ~ age + previous_treatment,
#'                      initial.predictor.method = "logistic",
#'                      ipcw.model = ~ age + previous_cost + previous_treatment,
#'                      tau0 = tau0,
#'                      higher.y = TRUE,
#'                      cv.n = 5, seed = 999, verbose = 1)
#'
#' # Try:
#' plot(catecv, ylab = "RMTL ratio of drug1 vs drug0 in each subgroup")
#' boxplot(catecv, ylab = "RMTL ratio of drug1 vs drug0 in each subgroup")
#' abc(catecv)
#'}
#' @export
#'
#' @importFrom stats as.formula coef glm median model.frame model.matrix model.response na.omit optim pchisq predict qnorm quantile sd var
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom survival Surv coxph coxph.detail survreg
#' @importFrom randomForestSRC rfsrc predict.rfsrc
#' @importFrom gbm gbm gbm.perf predict.gbm
#' @importFrom gam gam
#' @importFrom stringr str_replace
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom dplyr %>%
#'
catecvsurv <- function(data,
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
                       prop.multi = c(0, 1/3, 2/3, 1),
                       abc = TRUE,
                       train.prop = 3/4,
                       cv.n = 10,
                       error.max = 0.1,
                       max.iter = 5000,
                       surv.min = 0.025,
                       tree.depth = 2,
                       n.trees.rf = 1000,
                       n.trees.boosting = 200,
                       B = 3,
                       Kfold = 5,
                       error.maxNR = 1e-3,
                       max.iterNR = 150,
                       tune = c(0.5, 2),
                       seed = NULL,
                       plot.gbmperf = TRUE,
                       verbose = 0) {

  # Set seed for reproducibility
  set.seed(seed)

  if (verbose >= 1) t.start <- Sys.time()

  # Check arguments
  arg.checks(
    fun = "crossv", response = "survival", data = data, followup.time = followup.time, tau0 = tau0, surv.min = surv.min,
    higher.y = higher.y, score.method = score.method, abc = abc,
    prop.cutoff = prop.cutoff, prop.multi = prop.multi,
    ps.method = ps.method, minPS = minPS, maxPS = maxPS, ipcw.method = ipcw.method,
    train.prop = train.prop, cv.n = cv.n,
    error.max = error.max, max.iter = max.iter,
    initial.predictor.method = initial.predictor.method,
    tree.depth = tree.depth, n.trees.rf = n.trees.rf, n.trees.boosting = n.trees.boosting,
    B = B, Kfold = Kfold, plot.gbmperf = plot.gbmperf,
    error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune
  )


  # Pre-processing
  out <- data.preproc.surv(fun = "crossv", cate.model = cate.model, ps.model = ps.model, ipcw.model = ipcw.model, tau0 = tau0,
                           data = data, prop.cutoff = prop.cutoff, prop.multi = prop.multi,
                           ps.method = ps.method, initial.predictor.method = initial.predictor.method, response = "survival")
  y <- out$y
  d <- out$d
  trt <- out$trt
  x.ps <- out$x.ps
  x.cate <- out$x.cate
  x.ipcw <- out$x.ipcw
  prop.onlyhigh <- out$prop.onlyhigh
  prop.bi <- out$prop.bi
  prop.multi <- out$prop.multi
  tau0 <- out$tau0
  initial.predictor.method <- out$initial.predictor.method

  if (is.null(followup.time)) {
    yf <- NULL
  } else {
    yf <- data[[followup.time]]
  }

  # Check if tau0 is large enough, i.e., if tau0 is larger than the 50% quantile of the observed or censoring time.
  if (tau0 < median(y)) warning("It is recommended to increase tau0 close to the largest observed or censored time.")

  if (abc == TRUE) prop.abc <- seq(prop.onlyhigh[1], prop.onlyhigh[length(prop.onlyhigh)], length = 100)

  #### FUNCTION STARTS HERE ####
  n.subgroup.onlyhigh <- length(prop.onlyhigh)
  n.subgroup.bi <- length(prop.bi)
  n.subgroup.multi <- length(prop.multi) - 1

  result <- vector("list", 2 * length(score.method) + 1)
  names(result) <- c(paste0("ate.", score.method), paste0("hr.", score.method), "props")
  for (name in c(paste0("ate.", score.method), paste0("hr.", score.method))) {
    name1 <- gsub("\\..*", ".", name)
    result[[name]] <- vector("list", 6)
    names(result[[name]]) <- paste0(name1, c("est.train.high.cv", "est.train.low.cv",
                                             "est.valid.high.cv", "est.valid.low.cv",
                                             "est.train.group.cv", "est.valid.group.cv"))
    result[[name]][[paste0(name1, "est.train.high.cv")]] <-
      result[[name]][[paste0(name1, "est.valid.high.cv")]] <- matrix(NA, n.subgroup.onlyhigh, cv.n)
    result[[name]][[paste0(name1, "est.train.low.cv")]] <-
      result[[name]][[paste0(name1, "est.valid.low.cv")]] <- matrix(NA, n.subgroup.bi, cv.n)
    result[[name]][[paste0(name1, "est.train.group.cv")]] <-
      result[[name]][[paste0(name1, "est.valid.group.cv")]] <- matrix(NA, n.subgroup.multi, cv.n)

    # Add names to columns
    result[[name]] <- lapply(result[[name]], function(x) {
      colnames(x) <- paste0("cv", seq(1, cv.n))
      return(x)
    })

    # Add names to rows
    rownames(result[[name]][[paste0(name1, "est.train.high.cv")]]) <-
      rownames(result[[name]][[paste0(name1, "est.valid.high.cv")]]) <- paste0("prop", round(prop.onlyhigh, 2))
    rownames(result[[name]][[paste0(name1, "est.train.low.cv")]]) <-
      rownames(result[[name]][[paste0(name1, "est.valid.low.cv")]]) <- paste0("prop", round(1 - prop.bi, 2))
    rownames(result[[name]][[paste0(name1, "est.train.group.cv")]]) <-
      rownames(result[[name]][[paste0(name1, "est.valid.group.cv")]]) <- paste0("prop", round(prop.multi[-1], 2))
  }

  result$props <- vector("list", 3)
  names(result$props) <- c("prop.onlyhigh", "prop.bi", "prop.multi")
  result$overall.ate.train <- result$overall.hr.train <-
    result$overall.ate.valid <- result$overall.hr.valid <- rep(NA, cv.n)

  # Add names to errors/warnings
  result[["errors/warnings"]] <- vector("list", length(score.method))
  names(result[["errors/warnings"]]) <- score.method
  result[["errors/warnings"]] <- lapply(result[["errors/warnings"]],
                                        function(x) {
                                          x <- vector("list", 6)
                                          names(x) <- c("est.train.high.cv", "est.train.low.cv",
                                                        "est.valid.high.cv", "est.valid.low.cv",
                                                        "est.train.group.cv", "est.valid.group.cv")
                                          return(x)
                                        })

  if ("contrastReg" %in% score.method) converge.contrastReg.cv <- rep(NA, cv.n)

  # Save coefficients
  if (sum(score.method %in% c("poisson", "twoReg", "contrastReg")) > 0) {
    cf <- vector("list", length = 3)

    cf <- lapply(1:length(cf), function(x) {
      cf[[x]] <- matrix(NA, nrow = ncol(x.cate) + 1, ncol = cv.n)
      colnames(cf[[x]]) <- paste0("cv", seq(1, cv.n))
      rownames(cf[[x]]) <- c("(Intercept)", colnames(x.cate))
      return(cf[[x]])
    })
    names(cf) <- c("poisson", "twoReg", "contrastReg")
  }

  # Set progress bar
  if (verbose >= 1) pb <- txtProgressBar(min = 0, max = cv.n, style = 3)

  # Begin CV iteration
  for(cv.i in 1:cv.n) {

    ##### Split the data ------------------------------------------------------------------------------
    if (verbose >= 1){
      cat("\ncv =", cv.i, "\n")
      cat("  splitting the data..\n")
    }

    datacv <- balancesurv.split(y = y, d = d, trt = trt, x.cate = x.cate, x.ps = x.ps, x.ipcw = x.ipcw, yf = yf,
                                train.prop = train.prop, error.max = error.max, max.iter = max.iter)

    y.train <- datacv$y.train
    d.train <- datacv$d.train
    trt.train <- datacv$trt.train
    x.cate.train <- datacv$x.cate.train
    x.ps.train <- datacv$x.ps.train
    x.ipcw.train <- datacv$x.ipcw.train
    yf.train <- datacv$yf.train

    y.valid <- datacv$y.valid
    d.valid <- datacv$d.valid
    trt.valid <- datacv$trt.valid
    x.cate.valid <- datacv$x.cate.valid
    x.ps.valid <- datacv$x.ps.valid
    x.ipcw.valid <- datacv$x.ipcw.valid
    yf.valid <- datacv$yf.valid

    # Calculate the ATE in the training and validation sets (for abc, plot, boxplot)
    overall.train <- drsurv(y = y.train, d = d.train, x.cate = x.cate.train, x.ps = x.ps.train,
                            x.ipcw = x.ipcw.train, trt = trt.train, yf = yf.train, tau0 = tau0, surv.min = surv.min,
                            ps.method = ps.method, minPS = minPS, maxPS = maxPS, ipcw.method = ipcw.method)
    result$overall.ate.train[cv.i] <- exp(overall.train$log.rmtl.ratio)
    result$overall.hr.train[cv.i] <- exp(overall.train$log.hazard.ratio)

    overall.valid <- drsurv(y = y.valid, d = d.valid, x.cate = x.cate.valid, x.ps = x.ps.valid,
                            x.ipcw = x.ipcw.valid, trt = trt.valid, yf = yf.valid, tau0 = tau0, surv.min = surv.min,
                            ps.method = ps.method, minPS = minPS, maxPS = maxPS, ipcw.method = ipcw.method)
    result$overall.ate.valid[cv.i] <- exp(overall.valid$log.rmtl.ratio)
    result$overall.hr.valid[cv.i] <- exp(overall.valid$log.hazard.ratio)

    ####### Fit the interaction model --------------------------------------------------------------
    if (verbose >= 1) cat("  training..\n")
    fit <- intxsurv(y = y.train, d = d.train, trt = trt.train, x.cate = x.cate.train, x.ps = x.ps.train, x.ipcw = x.ipcw.train,
                    yf = yf.train, tau0 = tau0, surv.min = surv.min,
                    score.method = score.method,
                    ps.method = ps.method, minPS = minPS, maxPS = maxPS, ipcw.method = ipcw.method,
                    initial.predictor.method = initial.predictor.method,
                    tree.depth = tree.depth, n.trees.rf = n.trees.rf, n.trees.boosting = n.trees.boosting,
                    B = B, Kfold = Kfold, plot.gbmperf = plot.gbmperf,
                    error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune)


    if ("contrastReg" %in% score.method) {
      cf$contrastReg[, cv.i] <- fit$result.contrastReg$delta.contrastReg
      converge.contrastReg.cv[cv.i] <- fit$result.contrastReg$converge.contrastReg
    }

    if ("twoReg" %in% score.method) {
      cf$twoReg[, cv.i] <- fit$result.twoReg
    }

    if ("poisson" %in% score.method) {
      cf$poisson[, cv.i] <- fit$result.poisson
    }

    if ((initial.predictor.method == "boosting") & (fit$best.iter == n.trees.boosting) & (verbose == 2)) {
      warning(paste("The best boosting iteration was iteration number", n.trees.boosting, " out of ", n.trees.boosting, ". Consider increasing the maximum number of trees and turning on boosting performance plot (plot.gbmperf = TRUE).", sep = ""))
    }


    ####### Construct the score in training set ------------------------------------------
    fit.score.train <- scoresurv(fit = fit, x.cate = x.cate.train, tau0 = tau0, score.method = score.method)

    ####### Estimate the treatment effect in training set --------------------------------
    errors.train <- warnings.train <- c()
    if (prop.onlyhigh[n.subgroup.onlyhigh] == 1) {
      est.prop1 <- survCatch(drsurv(y = y.train, d = d.train, x.cate = x.cate.train, x.ps = x.ps.train, x.ipcw = x.ipcw.train,
                                    trt = trt.train, yf = yf.train,
                                    tau0 = tau0, surv.min = surv.min, ps.method = ps.method, minPS = minPS, maxPS = maxPS, ipcw.method = ipcw.method))
      est.prop1$ate.rmtl.high <- exp(est.prop1$log.rmtl.ratio)
      est.prop1$hr.high <- exp(est.prop1$log.hazard.ratio)
    } else {
      est.prop1 <- NULL
    }

    for (name in names(fit.score.train)) {
      score.train <- fit.score.train[[name]]
      name.method <- str_extract(name, "(?<=^score\\.).*$")

      est.bi <- estsurv.bilevel.subgroups(y = y.train, d = d.train, x.cate = x.cate.train, x.ps = x.ps.train, x.ipcw = x.ipcw.train,
                                          trt = trt.train, yf = yf.train,
                                          tau0 = tau0, surv.min = surv.min,
                                          score = score.train, higher.y = higher.y,
                                          prop = prop.bi, onlyhigh = FALSE,
                                          ps.method = ps.method, minPS = minPS, maxPS = maxPS, ipcw.method = ipcw.method)
      result[[str_replace(name, "score", "ate")]]$ate.est.train.high.cv[, cv.i] <- c(est.bi$ate.rmtl.high, est.prop1$ate.rmtl.high)
      result[[str_replace(name, "score", "ate")]]$ate.est.train.low.cv[, cv.i] <- est.bi$ate.rmtl.low
      result[[str_replace(name, "score", "hr")]]$hr.est.train.high.cv[, cv.i] <- c(est.bi$hr.high, est.prop1$hr.high)
      result[[str_replace(name, "score", "hr")]]$hr.est.train.low.cv[, cv.i] <- est.bi$hr.low

      est.bi$err.high$`prop 1.0` <- est.prop1$errors
      err.bi.high <- est.bi$err.high[!sapply(est.bi$err.high, is.null)]
      err.bi.low <- est.bi$err.low[!sapply(est.bi$err.low, is.null)]
      if (length(err.bi.high) == 0) err.bi.high <- NULL
      if (length(err.bi.low) == 0) err.bi.low <- NULL

      est.bi$warn.high$`prop 1.0` <- est.prop1$warnings
      warn.bi.high <- est.bi$warn.high[!sapply(est.bi$warn.high, is.null)]
      warn.bi.low <- est.bi$warn.low[!sapply(est.bi$warn.low, is.null)]
      if (length(warn.bi.high) == 0) warn.bi.high <- NULL
      if (length(warn.bi.low) == 0) warn.bi.low <- NULL
      result[['errors/warnings']][[name.method]]$est.train.high.cv[[paste0("cv", cv.i)]]$errors <- err.bi.high
      result[['errors/warnings']][[name.method]]$est.train.high.cv[[paste0("cv", cv.i)]]$warnings <- warn.bi.high
      result[['errors/warnings']][[name.method]]$est.train.low.cv[[paste0("cv", cv.i)]]$errors <- err.bi.low
      result[['errors/warnings']][[name.method]]$est.train.low.cv[[paste0("cv", cv.i)]]$warnings <- warn.bi.low

      est.multi <- estsurv.multilevel.subgroups(y = y.train, d = d.train, x.cate = x.cate.train, x.ps = x.ps.train, x.ipcw = x.ipcw.train,
                                                trt = trt.train, yf = yf.train,
                                                tau0 = tau0, surv.min = surv.min,
                                                score = score.train, higher.y = higher.y,
                                                prop = prop.multi,
                                                ps.method = ps.method, minPS = minPS, maxPS = maxPS, ipcw.method = ipcw.method)
      result[[str_replace(name, "score", "ate")]]$ate.est.train.group.cv[, cv.i] <- est.multi$ate.rmtl
      result[[str_replace(name, "score", "hr")]]$hr.est.train.group.cv[, cv.i] <- est.multi$hr

      err.multi <- est.multi$err[!sapply(est.multi$err, is.null)]
      warn.multi <- est.multi$warn[!sapply(est.multi$warn, is.null)]
      if (length(err.multi) == 0) err.multi <- NULL
      if (length(warn.multi) == 0) warn.multi <- NULL
      result[['errors/warnings']][[name.method]]$est.train.group.cv[[paste0("cv", cv.i)]]$errors <- err.multi
      result[['errors/warnings']][[name.method]]$est.train.group.cv[[paste0("cv", cv.i)]]$warnings <- warn.multi

      if (any(!is.null(c(err.bi.high, err.bi.low, err.multi)))) errors.train <- c(errors.train, name.method)
      if (any(!is.null(c(warn.bi.high, warn.bi.low, warn.multi)))) warnings.train <- c(warnings.train, name.method)
    } # end of for (name in names(fit.score.train)) {}

    if (verbose == 2){
      if (length(errors.train) != 0) {
        if (verbose == 2) cat(paste0('    Warning: Error(s) occurred when estimating the ATEs in the nested subgroup using "', paste0(errors.train, collapse = '", "'), '";\n    return NAs in the corresponding subgroup.'),'\n')
        warning(paste0('Error(s) occurred when estimating the ATEs in the nested subgroup in the training set using "', paste0(errors.train, collapse = '", "'), '" in cross-validation iteration ', cv.i, ". NAs are returned for RMTL ratio and HR in the corresponding subgroup; see 'errors/warnings'."))
      }

      if (length(warnings.train) != 0) {
        if (verbose == 2) cat(paste0('     Warning: Warning(s) occurred when estimating the ATEs in the nested subgroup using "', paste0(warnings.train, collapse = '", "'), '"'),'\n')
        warning(paste0('Warning(s) occurred when estimating the ATEs in the nested subgroup in the training set using "', paste0(warnings.train, collapse = '", "'), '" in cross-validation iteration ', cv.i, "; see 'errors/warnings'."))
      }
    }

    ####### Construct the score in validation set ------------------------------------------
    if (verbose >= 1) cat("  validating..\n")

    fit.score.valid <- scoresurv(fit = fit, x.cate = x.cate.valid, tau0 = tau0, score.method = score.method)

    ####### Estimate the treatment effect in validation set --------------------------------
    errors.valid <- warnings.valid  <- c()

    if (prop.onlyhigh[n.subgroup.onlyhigh] == 1) {
      est.prop1 <- survCatch(drsurv(y = y.valid, d = d.valid, x.cate = x.cate.valid, x.ps = x.ps.valid, x.ipcw = x.ipcw.valid,
                                    trt = trt.valid, yf = yf.valid,
                                    tau0 = tau0, surv.min = surv.min, ps.method = ps.method, minPS = minPS, maxPS = maxPS, ipcw.method = ipcw.method))
      est.prop1$ate.rmtl.high <- exp(est.prop1$log.rmtl.ratio)
      est.prop1$hr.high <- exp(est.prop1$log.hazard.ratio)
    } else {
      est.prop1 <- NULL
    }

    for (name in names(fit.score.valid)) {
      score.valid <- fit.score.valid[[name]]
      name.method <- str_extract(name, "(?<=^score\\.).*$")

      est.bi <- estsurv.bilevel.subgroups(y = y.valid, d = d.valid, x.cate = x.cate.valid, x.ps = x.ps.valid, x.ipcw = x.ipcw.valid,
                                          trt = trt.valid, yf = yf.valid,
                                          tau0 = tau0, surv.min = surv.min,
                                          score = score.valid, higher.y = higher.y,
                                          prop = prop.bi, onlyhigh = FALSE,
                                          ps.method = ps.method, minPS = minPS, maxPS = maxPS, ipcw.method = ipcw.method)
      result[[str_replace(name, "score", "ate")]]$ate.est.valid.high.cv[, cv.i] <- c(est.bi$ate.rmtl.high, est.prop1$ate.rmtl.high)
      result[[str_replace(name, "score", "ate")]]$ate.est.valid.low.cv[, cv.i] <- est.bi$ate.rmtl.low
      result[[str_replace(name, "score", "hr")]]$hr.est.valid.high.cv[, cv.i] <- c(est.bi$hr.high, est.prop1$hr.high)
      result[[str_replace(name, "score", "hr")]]$hr.est.valid.low.cv[, cv.i] <- est.bi$hr.low


      est.bi$err.high$`prop 1.0` <- est.prop1$errors
      err.bi.high <- est.bi$err.high[!sapply(est.bi$err.high, is.null)]
      err.bi.low <- est.bi$err.low[!sapply(est.bi$err.low, is.null)]
      if (length(err.bi.high) == 0) err.bi.high <- NULL
      if (length(err.bi.low) == 0) err.bi.low <- NULL

      est.bi$warn.high$`prop 1.0` <- est.prop1$warnings
      warn.bi.high <- est.bi$warn.high[!sapply(est.bi$warn.high, is.null)]
      warn.bi.low <- est.bi$warn.low[!sapply(est.bi$warn.low, is.null)]
      if (length(warn.bi.high) == 0) warn.bi.high <- NULL
      if (length(warn.bi.low) == 0) warn.bi.low <- NULL
      result[['errors/warnings']][[name.method]]$est.valid.high.cv[[paste0("cv", cv.i)]]$errors <- err.bi.high
      result[['errors/warnings']][[name.method]]$est.valid.high.cv[[paste0("cv", cv.i)]]$warnings <- warn.bi.high
      result[['errors/warnings']][[name.method]]$est.valid.low.cv[[paste0("cv", cv.i)]]$errors <- err.bi.low
      result[['errors/warnings']][[name.method]]$est.valid.low.cv[[paste0("cv", cv.i)]]$warnings <- warn.bi.low


      est.multi <- estsurv.multilevel.subgroups(y = y.valid, d = d.valid, x.cate = x.cate.valid, x.ps = x.ps.valid, x.ipcw = x.ipcw.valid,
                                                trt = trt.valid, yf = yf.valid,
                                                tau0 = tau0, surv.min = surv.min,
                                                score = score.valid, higher.y = higher.y,
                                                prop = prop.multi,
                                                ps.method = ps.method, minPS = minPS, maxPS = maxPS, ipcw.method = ipcw.method)
      result[[str_replace(name, "score", "ate")]]$ate.est.valid.group.cv[, cv.i] <- est.multi$ate.rmtl
      result[[str_replace(name, "score", "hr")]]$hr.est.valid.group.cv[, cv.i] <- est.multi$hr

      err.multi <- est.multi$err[!sapply(est.multi$err, is.null)]
      warn.multi <- est.multi$warn[!sapply(est.multi$warn, is.null)]
      if (length(err.multi) == 0) err.multi <- NULL
      if (length(warn.multi) == 0) warn.multi <- NULL
      result[['errors/warnings']][[name.method]]$est.valid.group.cv[[paste0("cv", cv.i)]]$errors <- err.multi
      result[['errors/warnings']][[name.method]]$est.valid.group.cv[[paste0("cv", cv.i)]]$warnings <- warn.multi

      if (any(!is.null(c(err.bi.high, err.bi.low, err.multi)))) errors.valid <- c(errors.valid, name.method)
      if (any(!is.null(c(warn.bi.high, warn.bi.low, warn.multi)))) warnings.valid <- c(warnings.valid, name.method)

      if (abc == TRUE) {
        temp.abc <- log(estsurv.bilevel.subgroups(y = y.valid, d = d.valid, x.cate = x.cate.valid, x.ps = x.ps.valid, x.ipcw = x.ipcw.valid,
                                                  trt = trt.valid, yf = yf.valid, tau0 = tau0, surv.min = surv.min,
                                                  score = score.valid, higher.y = higher.y, prop = prop.abc, onlyhigh = TRUE,
                                                  ps.method = ps.method, minPS = minPS, maxPS = maxPS, ipcw.method = ipcw.method)$ate.rmtl.high)

        if (higher.y == FALSE) {
          # if higher.y=FALSE, higher scores are better and ABC requires subtracting the ATE from validation curve
          ref <- temp.abc - log(result$overall.ate.valid[cv.i])
        } else {
          #  if higher.y=TRUE, lower scores are better, ABC requires subtracting the validation curve from ATE
          ref <- log(result$overall.ate.valid[cv.i]) - temp.abc
        }
        result[[str_replace(name, "score", "ate")]]$abc.valid[cv.i] <-
          auc(x = prop.abc, y = ref, from = prop.abc[1], to = prop.abc[length(prop.abc)], type = "spline")
      } # end of if (abc == TRUE) {}
    } # end of for (name in names(fit.score.train)) {}


    if (verbose == 2){
      if (length(errors.valid) != 0) {
        if (verbose == 2) cat(paste0('     Warning: Error(s) occurred when estimating the ATEs in the nested subgroup using "', paste0(errors.valid, collapse = '", "'), '";\n    return NAs in the corresponding subgroup.'),'\n')
        warning(paste0('Error(s) occurred when estimating the ATEs in the nested subgroup in the validation set using "', paste0(errors.valid, collapse = '", "'), '" in cross-validation iteration ', cv.i, ". NAs are returned for RMTL ratio and HR in the corresponding subgroup; see 'errors/warnings'."))
      }

      if (length(warnings.valid) != 0) {
        if (verbose == 2) cat(paste0('     Warning: Warning(s) occurred when estimating the ATEs in the nested subgroup using "', paste0(warnings.valid, collapse = '", "'), '"'),'\n')
        warning(paste0('Warning(s) occurred when estimating the ATEs in the nested subgroup in the validation set using "', paste0(warnings.valid, collapse = '", "'), '" in cross-validation iteration ', cv.i, "; see 'errors/warnings'."))
      }
    }

    if (verbose >= 1) {
      if ("contrastReg" %in% score.method) {
        if (converge.contrastReg.cv[cv.i]) {
          cat("  Contrast regression converged.\n")
        } else {
          cat("  Contrast regression did not converge.\n")
        }
      }
      cat("  ", date(), "\n")
      setTxtProgressBar(pb, cv.i)
    }
  } # end of cv loop

  if ("contrastReg" %in% score.method) {
    result$ate.contrastReg <- append(result$ate.contrastReg,
                                     list(coefficients = cf$contrastReg,
                                          converge.contrastReg.cv = converge.contrastReg.cv))
    if (any(is.na(converge.contrastReg.cv)) & (verbose == 2)) {
      warning("Contrast regression algorithm did not converge in at least one cross-validation iteration.")
    }
  }

  if ("twoReg" %in% score.method) {
    result$ate.twoReg <- append(result$ate.twoReg,
                                list(coefficients = cf$twoReg))
  }

  if ("poisson" %in% score.method) {
    result$ate.poisson <- append(result$ate.poisson,
                                 list(coefficients = cf$poisson))
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
  result$higher.y <- higher.y
  result$abc <- abc
  result$cv.n <- cv.n
  result$response <- "survival"
  result$formulas <- list(cate.model = cate.model, ps.model = ps.model, trt_labels = out$cat.trt)

  class(result) <- "precmed"

  return(result)
}


#' Cross-validation of the conditional average treatment effect (CATE) score
#' for count outcomes
#'
#' Provides doubly robust estimation of the average treatment effect (ATE) in
#' nested and mutually exclusive subgroups of patients defined by an estimated
#' conditional average treatment effect (CATE) score via cross-validation (CV).
#' The CATE score can be estimated with up to 5 methods among the following:
#' Poisson regression, boosting, two regressions, contrast regression, and
#' negative binomial (see \code{score.method}).
#'
#' @param cate.model A formula describing the outcome model to be fitted.
#' The outcome must appear on the left-hand side.
#' @param ps.model A formula describing the propensity score model to be fitted.
#' The treatment must appear on the left-hand side. The treatment must be a
#' numeric vector coded as 0 or 1. If data are from a randomized trial, specify
#' \code{ps.model} as an intercept-only model.
#' @param data A data frame containing the variables in the outcome and
#' propensity score model; a data frame with \code{n} rows (1 row per
#' observation).
#' @param score.method A vector of one or multiple methods to estimate the CATE
#' score. Allowed values are: \code{'boosting'}, \code{'poisson'},
#' \code{'twoReg'}, \code{'contrastReg'}, and \code{'negBin'}.
#' @param higher.y A logical value indicating whether higher (\code{TRUE}) or
#' lower (\code{FALSE}) values of the outcome are more desirable. Default is
#' \code{TRUE}.
#' @param abc A logical value indicating whether the area between curves (ABC)
#' should be calculated at each cross-validation iterations, for each
#' \code{score.method}. Default is \code{TRUE}.
#' @param prop.cutoff A vector of numerical values between 0 and 1 specifying
#' percentiles of the estimated log CATE scores to define nested subgroups. Each
#' element represents the cutoff to separate observations in nested subgroups
#' (below vs above cutoff). The length of \code{prop.cutoff} is the number of
#' nested subgroups. An equally-spaced sequence of proportions ending with 1 is
#' recommended. Default is \code{seq(0.5, 1, length = 6)}.
#' @param prop.multi A vector of numerical values between 0 and 1 specifying
#' percentiles of the estimated log CATE scores to define mutually exclusive
#' subgroups. It should start with 0, end with 1, and be of
#' \code{length(prop.multi) > 2}. Each element represents the cutoff to separate
#' the observations into \code{length(prop.multi) - 1} mutually exclusive
#' subgroups. Default is \code{c(0, 1/3, 2/3, 1)}.
#' @param ps.method A character value for the method to estimate the propensity
#' score. Allowed values include one of: \code{'glm'} for logistic regression
#' with main effects only (default), or \code{'lasso'} for a logistic regression
#' with main effects and LASSO penalization on two-way interactions (added to
#' the model if interactions are not specified in \code{ps.model}). Relevant
#' only when \code{ps.model} has more than one variable.
#' @param minPS A numerical value between 0 and 1 below which estimated
#' propensity scores should be truncated. Default is \code{0.01}.
#' @param maxPS A numerical value between 0 and 1 above which estimated
#' propensity scores should be truncated. Must be strictly greater than
#' \code{minPS}. Default is \code{0.99}.
#' @param train.prop A numerical value between 0 and 1 indicating the proportion
#' of total data used for training. Default is \code{3/4}.
#' @param cv.n A positive integer value indicating the number of
#' cross-validation iterations. Default is \code{10}.
#' @param error.max A numerical value > 0 indicating the tolerance (maximum
#' value of error) for the largest standardized absolute difference in the
#' covariate distributions or in the doubly robust estimated rate ratios between
#' the training and validation sets. This is used to define a balanced
#' training-validation splitting. Default is \code{0.1}.
#' @param max.iter A positive integer value indicating the maximum number of
#' iterations when searching for a balanced training-validation split. Default
#' is \code{5,000}.
#' @param initial.predictor.method A character vector for the method used to get
#' initial outcome predictions conditional on the covariates in
#' \code{cate.model}. Only applies when \code{score.method} includes
#' \code{'twoReg'} or \code{'contrastReg'}. Allowed values include one of
#' \code{'poisson'} (fastest), \code{'boosting'} and \code{'gam'}. Default is
#' \code{'boosting'}.
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
#' \item{\code{overall.ate.train}: }{A vector of numerical values of length \code{cv.n}.
#' The ith element contains the ATE in the training set of the ith cross-validation
#' iteration, estimated with the doubly robust estimator.}
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
#' @references
#' Yadlowsky, S., Pellegrini, F., Lionetto, F., Braune, S., & Tian, L. (2020).
#' \emph{Estimation and validation of ratio-based conditional average treatment
#' effects using observational data. Journal of the American Statistical
#' Association, 1-18.}
#' \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1772080}
#'
#' @seealso
#' \code{\link{plot.precmed}()}, \code{\link{boxplot.precmed}()},
#' \code{\link{abc}()} methods for \code{"precmed"} objects,
#' and \code{\link{catefitcount}()} function.
#'
#' @examples
#' \donttest{
#' catecv <- catecvcount(data = countExample,
#'                       score.method = "poisson",
#'                       cate.model = y ~ age + female + previous_treatment +
#'                                    previous_cost + previous_number_relapses +
#'                                    offset(log(years)),
#'                       ps.model = trt ~ age + previous_treatment,
#'                       higher.y = FALSE,
#'                       cv.n = 5,
#'                       seed = 999,
#'                       plot.gbmperf = FALSE,
#'                       verbose = 1)
#'
#' plot(catecv, ylab = "Rate ratio of drug1 vs drug0 in each subgroup")
#' boxplot(catecv, ylab = "Rate ratio of drug1 vs drug0 in each subgroup")
#' abc(catecv)
#' }
#' @export
#'
#' @importFrom graphics hist lines
#' @importFrom stats as.formula coef glm median model.frame model.matrix model.offset model.response na.omit optim pchisq predict qnorm quantile sd var
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom gbm gbm gbm.perf
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
                        verbose = 0,
                        ...) {

  # Set seed for reproducibility
  set.seed(seed)

  if (verbose >= 1) t.start <- Sys.time()

  #### CHECK ARGUMENTS ####
  arg.checks(
    fun = "crossv", response = "count", data = data, higher.y = higher.y, score.method = score.method, abc = abc,
    prop.cutoff = prop.cutoff, prop.multi = prop.multi,
    ps.method = ps.method, minPS = minPS, maxPS = maxPS,
    train.prop = train.prop, cv.n = cv.n,
    error.max = error.max, max.iter = max.iter,
    initial.predictor.method = initial.predictor.method,
    tree.depth = tree.depth, n.trees.boosting = n.trees.boosting, B = B, Kfold = Kfold, plot.gbmperf = plot.gbmperf,
    error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune
  )

  #### PRE-PROCESSING ####
  out <- data.preproc(fun = "crossv", cate.model = cate.model, ps.model = ps.model,
                      data = data, prop.cutoff = prop.cutoff,
                      prop.multi = prop.multi, ps.method = ps.method,
                      initial.predictor.method = initial.predictor.method)
  y <- out$y
  trt <- out$trt
  x.ps <- out$x.ps
  x.cate <- out$x.cate
  time <- out$time
  prop.onlyhigh <- out$prop.onlyhigh
  prop.bi <- out$prop.bi
  prop.multi <- out$prop.multi
  initial.predictor.method <- out$initial.predictor.method

  if (abc) prop.abc <- seq(prop.onlyhigh[1], prop.onlyhigh[length(prop.onlyhigh)], length = 100)

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
  for (cv.i in seq(cv.n)) {

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

    if (verbose == 2) {
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
#' @param prop.cutoff A vector of numerical values (in `(0, 1]`) specifying percentiles of the
#' estimated CATE scores to define nested subgroups. Each element represents the cutoff to
#' separate observations in nested subgroups (below vs above cutoff).
#' The length of \code{prop.cutoff} is the number of nested subgroups.
#' An equally-spaced sequence of proportions ending with 1 is recommended.
#' Default is \code{seq(0.5, 1, length = 6)}.
#' @param prop.multi A vector of numerical values (in `[0, 1]`) specifying percentiles of the
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
#' @param minPS A numerical value (in `[0, 1]`) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in `(0, 1]`) above which estimated propensity scores should be
#' truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.
#' @param train.prop A numerical value (in `(0, 1)`) indicating the proportion of total data used
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
#' \code{2} means progress bar, run time, and all errors and warnings. Default is \code{0}.
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
#' \item{\code{ate.randomForest}: }{A list of ATE output measured by the RMTL ratio if
#' \code{score.method} includes \code{'randomForest'}:}
#' \item{\code{ate.gam}: }{A list of results similar to \code{ate.gaussian} output
#' if \code{score.method} includes \code{'gam'}.}
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
#' and \code{\link{catefitmean}()} function.
#'
#' @examples
#' # Not implemented yet!
#'
#' @importFrom graphics hist lines
#' @importFrom stats as.formula coef glm median model.frame model.matrix model.offset model.response na.omit optim pchisq predict qnorm quantile sd var
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom gbm gbm gbm.perf
#' @importFrom stringr str_replace str_extract str_detect
#' @importFrom dplyr %>%
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom MASS glm.nb ginv
#' @importFrom randomForestSRC rfsrc predict.rfsrc

catecvmean <- function(data,
                       score.method,
                       cate.model,
                       ps.model,
                       ps.method = "glm",
                       init.model = NULL,
                       initial.predictor.method = "boosting",
                       minPS = 0.01,
                       maxPS = 0.99,
                       higher.y = TRUE,
                       prop.cutoff = seq(0.5, 1, length = 6),
                       prop.multi = c(0, 1/3, 2/3, 1),
                       abc = TRUE,
                       train.prop = 3/4,
                       cv.n = 10,
                       error.max = 0.1,
                       max.iter = 5000,
                       xvar.smooth.score = NULL,
                       xvar.smooth.init = NULL,
                       tree.depth = 2,
                       n.trees.rf = 1000,
                       n.trees.boosting = 200,
                       B = 3,
                       Kfold = 6,
                       plot.gbmperf = TRUE,
                       error.maxNR = 1e-3,
                       tune = c(0.5, 2),
                       seed = NULL,
                       verbose = 0,
                       ...) {

  stop("This functionality is not implemented yet")

  # TODO: now score.method has no default (mandatory argument)

  # Set seed for reproducibility
  set.seed(seed)

  if (verbose >= 1) t.start <- Sys.time()

  #### CHECK ARGUMENTS ####
  ##TODO: Insert n.trees.rf
  arg.checks(
    fun = "crossv", response = "continuous", data = data, higher.y = higher.y, score.method = score.method, abc = abc,
    prop.cutoff = prop.cutoff, prop.multi = prop.multi,
    ps.method = ps.method, minPS = minPS, maxPS = maxPS,
    train.prop = train.prop, cv.n = cv.n,
    error.max = error.max, max.iter = max.iter,
    initial.predictor.method = initial.predictor.method,
    tree.depth = tree.depth, n.trees.boosting = n.trees.boosting, B = B, Kfold = Kfold, plot.gbmperf = plot.gbmperf,
    error.maxNR = error.maxNR, tune = tune
  )

  #### PRE-PROCESSING ####
  out <- data.preproc.mean(fun = "crossv", cate.model = cate.model, init.model = init.model, ps.model = ps.model,
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
        if (verbose == 2) cat(paste0('    Warning: Error(s) occurred when fitting "', paste0(names(fit$err.fit), collapse = '", "'), '";\n    return NAs in the corresponding parameters'),'\n')
        warning(paste0('Warning: Error(s) occurred when fitting "', paste0(names(fit$err.fit), collapse = '", "'), '" in cross-validation iteration ', cv.i, ". NAs are returned for the corresponding parameters; see 'errors/warnings'."))
        result[['errors/warnings']][[names(fit$err.fit)]][[paste0("cv", cv.i)]]$errors <- fit$err.fit
      }

      if (length(names(fit$warn.fit)) != 0) {
        if (verbose == 2) cat(paste0('    Warning: Warning(s) occurred when fitting "', paste0(names(fit$warn.fit), collapse = '", "'), '"'),'\n')
        warning(paste0('Warning: Warning(s) occurred when fitting "', paste0(names(fit$warn.fit), collapse = '", "'), '" in cross-validation iteration ', cv.i, "; see 'errors/warnings'."))
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
          auc(x = prop.abc, y = ref, from = prop.abc[1], to = prop.abc[length(prop.abc)], type = "spline")
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


#' Compute the area between curves from the \code{"precmed"} object
#'
#' Compute the area between curves (ABC) for each scoring method in the \code{"precmed"} object.
#' This should be run only after results of \code{\link{catecv}()} have been obtained.
#'
#' @param x An object of class \code{"precmed"}.
#' @param ... Additional arguments (currently unused).
#'
#'
#' @return Returns a matrix of numeric values with number of columns equal to the number cross-validation
#' iteration and number of rows equal to the number of scoring methods in \code{x}.
#'
#' @details The ABC is the area between a validation curve and the overall ATE in the validation set.
#' It is calculated for each scoring method separately. Higher ABC values are preferable as they
#' indicate that more treatment effect heterogeneity is captured by the scoring method.
#' Negative values of ABC are possible if segments of the validation curve cross the overall ATE line.
#' The ABC is calculated with the \code{\link{auc}()} in \code{utility.R} with a natural
#' cubic spline interpolation. The calculation of the ABC is always based on validation curves based on
#' 100 proportions equally spaced from \code{min(prop.cutoff)} to \code{max(prop.cutoff)}.
#'
#' The ABC is a metric to help users select the best scoring method in terms of capturing treatment
#' effect heterogeneity in the data. It should be used in complement to the visual inspection of
#' the validation curves in the validation set in \code{\link{plot}()}.
#'
#' @references Zhao, L., Tian, L., Cai, T., Claggett, B., & Wei, L. J. (2013).
#' \emph{Effectively selecting a target population for a future comparative study.
#' Journal of the American Statistical Association, 108(502), 527-539.}
#'
#' @seealso \code{\link{catecv}()} function and \code{\link{plot}()}, \code{\link{boxplot}()} methods for
#' \code{"precmed"} objects.
#'
#' @examples
#' \donttest{
#' # Count outcome
#' cv_count <- catecv(response = "count",
#'                    data = countExample,
#'                    score.method = "poisson",
#'                    cate.model = y ~ age + female + previous_treatment +
#'                                 previous_cost + previous_number_relapses +
#'                                 offset(log(years)),
#'                    ps.model = trt ~ age + previous_treatment,
#'                    higher.y = FALSE, cv.n = 5, verbose = 1)
#'
#' # ABC of the validation curves for each method and each CV iteration
#' abc(cv_count)
#'
#' # Survival outcome
#' library(survival)
#' cv_surv <- catecv(response = "survival",
#'                   data = survivalExample,
#'                   score.method = c("poisson", "randomForest"),
#'                   cate.model = Surv(y, d) ~ age + female + previous_cost +
#'                                previous_number_relapses,
#'                   ps.model = trt ~ age + previous_treatment,
#'                   higher.y = FALSE,
#'                   cv.n = 5)
#'
#' # ABC of the validation curves for each method and each CV iteration
#' abc(cv_surv)
#'
#' }
#'
#' @export
abc <- function(x, ...) UseMethod("abc", x)

#' @export
#' @rdname abc
abc.default <- function(x, ...){
  warning(paste("ABC does not know how to handle object of class ",
                class(x),
                "and can only be used on class precmed"))
}

#' Compute the area between curves from the \code{"precmed"} object
#'
#' Compute the area between curves (ABC) for each scoring method in the \code{"precmed"} object.
#' This should be run only after results of \code{\link{catecv}()} have been obtained.
#'
#' @param x An object of class \code{"precmed"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Returns a matrix of numeric values with number of columns equal to the number cross-validation
#' iteration and number of rows equal to the number of scoring methods in \code{x}.
#'
#' @details The ABC is the area between a validation curve and the overall ATE in the validation set.
#' It is calculated for each scoring method separately. Higher ABC values are preferable as they
#' indicate that more treatment effect heterogeneity is captured by the scoring method.
#' Negative values of ABC are possible if segments of the validation curve cross the overall ATE line.
#' The ABC is calculated with the \code{\link{auc}()} in \code{utility.R} with a natural
#' cubic spline interpolation. The calculation of the ABC is always based on validation curves based on
#' 100 proportions equally spaced from \code{min(prop.cutoff)} to \code{max(prop.cutoff)}.
#'
#' The ABC is a metric to help users select the best scoring method in terms of capturing treatment
#' effect heterogeneity in the data. It should be used in complement to the visual inspection of
#' the validation curves in the validation set in \code{\link{plot}()}.
#'
#' @references Zhao, L., Tian, L., Cai, T., Claggett, B., & Wei, L. J. (2013).
#' \emph{Effectively selecting a target population for a future comparative study.
#' Journal of the American Statistical Association, 108(502), 527-539.}
#'
#' @seealso \code{\link{catecv}()} function and \code{\link{plot}()}, \code{\link{boxplot}()} methods for
#' \code{"precmed"} objects.
#'
#' @examples
#' \donttest{
#' # Count outcome
#' cv_count <- catecv(response = "count",
#'                    data = countExample,
#'                    score.method = "poisson",
#'                    cate.model = y ~ age + female + previous_treatment +
#'                                 previous_cost + previous_number_relapses +
#'                                 offset(log(years)),
#'                    ps.model = trt ~ age + previous_treatment,
#'                    higher.y = FALSE, cv.n = 5, verbose = 1)
#'
#' # ABC of the validation curves for each method and each CV iteration
#' abc(cv_count)
#'
#' # Survival outcome
#' library(survival)
#' cv_surv <- catecv(response = "survival",
#'                   data = survivalExample,
#'                   score.method = c("poisson", "randomForest"),
#'                   cate.model = Surv(y, d) ~ age + female + previous_cost +
#'                                previous_number_relapses,
#'                   ps.model = trt ~ age + previous_treatment,
#'                   higher.y = FALSE,
#'                   cv.n = 5)
#'
#' # ABC of the validation curves for each method and each CV iteration
#' abc(cv_surv)
#'
#' }
#'
#' @export
#'
#' @importFrom stringr str_extract
abc.precmed <- function(x, ...) {

  # Check the value of abc (must have area between curves values)
  if (!x$abc) {
    stop("Area between curves (ABC) must have been calculated in x by setting abc = TRUE in catecv().")
  }

  score.method <- str_extract(names(x), "(?<=^ate\\.).*$")
  score.method <- score.method[!is.na(score.method)]
  m <- length(score.method)

  # Value of ABC
  value.abc <- matrix(NA, nrow = m, ncol = x$cv.n)
  for (i in seq(m)) {
    value.abc[i,] <- x[[paste0("ate.", score.method[i])]]$abc.valid
  }
  rownames(value.abc) <- score.method
  colnames(value.abc) <- paste0("cv", 1:(x$cv.n))

  return(value.abc)
}

#' Two side-by-side line plots of validation curves from the \code{"precmed"} object
#'
#' Provides validation curves in two side-by-side plots, visualizing the estimated ATEs in a series
#' of nested subgroups in the training set and validation set separately, where each line represents
#' one scoring method specified in \code{\link{catecv}()} or \code{\link{catecvmean}()}. This should be run
#' only after results of \code{\link{catecv}()} or \code{\link{catecvmean}()} have been obtained.
#'
#' @param x An object of class \code{"precmed"}.
#' @param cv.i A positive integer indicating the index of the CV iteration results to be plotted.
#' Allowed values are: a positive integer \eqn{<=} \code{cv.n} in \code{\link{catecv}()} or
#' \code{NULL}. If \code{cv.i = NULL}, the results across all CV iterations are combined according
#' to \code{combine} and then plotted. Default is \code{NULL}.
#' @param combine A character value indicating how to combine the estimated ATEs across all CV
#' iterations into a validation curve for each nested subgroup, separately for the training and
#' validation results. Allowed values are: \code{'mean'} or \code{'median'}. Used only if
#' \code{cv.i = NULL}. Default is \code{'mean'}.
#' @param show.abc A logical value indicating whether to show the ABC statistics in the validation set. Used
#' only if \code{x$abc = TRUE} and \code{xlim} is not limited to a smaller range (i.e., \code{xlim = NULL} or
#' equal to the entire \code{x$prop.onlyhigh} range). If \code{cv.i} is NULL, ABC statistics will be based
#' on the combined CV iterations. If \code{cv.i} is an integer, ABC statistics will be based solely on that
#' CV iteration. Default is \code{TRUE}.
#' @param valid.only A logical value indicating whether only the validation curves in the validation set
#' should be plotted (\code{TRUE}). Otherwise, the validation curves in both the training and validation
#' sets are plotted side-by-side (\code{FALSE}). Default is \code{FALSE}.
#' @param plot.hr A logical value indicating whether the hazard ratios should be plotted in the
#' validation curves (\code{TRUE}). Otherwise, the restricted mean time lost is plotted (\code{FALSE}).
#' This argument is only applicable to survival outcomes. Default is \code{FALSE}.
#' @param ylab A character value for the y-axis label to describe what the ATE is. Default is \code{NULL},
#' which creates a default y-axis label based on available data.
#' @param legend.position A character value for the legend position argument to be passed to \code{ggplot}
#' object. Default is \code{'bottom'}.
#' @param xlim A numeric value for the range of the x-axis. Default is \code{NULL}, which means there is no
#' range specified.
#' @param title The text for the title
#' @param theme Defaults to \code{theme_classic()}. Other options include \code{theme_grey()}, \code{theme_bw()}, \code{theme_light()}, \code{theme_dark()}, and \code{theme_void()}
#' @param ... Other parameters
#'
#' @return Returns two side-by-side line plots, one of which shows the validation curves of the training
#' sets and the other the validation curves in the validation sets. A gray horizontal dashed line of
#' overall ATE is included as a reference. ABC statistics will be added to the legend if
#' \code{show.abc = TRUE}.
#'
#' @details \code{\link{plot}()} takes in outputs from \code{\link{catecv}()} and generates two plots
#' of validation curves side-by-side, one for the training set and one for validation set.
#' Separate validation curves are produced for each scoring method specified via \code{score.method}
#' in \code{\link{catecv}()} or \code{\link{catecvmean}()}.
#'
#' The validation curves (and ABC statistics, if applicable) can help compare the performance of
#' different scoring methods in terms of discerning potential treatment heterogeneity in subgroups
#' with internal validation. Steeper validation curves in the validation set suggest presence of
#' treatment effect heterogeneity (and the ability of the scoring methods to capture it) while flat
#' validation curves indicate absence of treatment effect heterogeneity (or inability of the scoring method
#' to capture it).
#'
#' @references Yadlowsky, S., Pellegrini, F., Lionetto, F., Braune, S., & Tian, L. (2020).
#' \emph{Estimation and validation of ratio-based conditional average treatment effects using
#' observational data. Journal of the American Statistical Association, 1-18.}
#' \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1772080}
#'
#' @seealso \code{\link{abc}()} and \code{\link{boxplot}()} for \code{"precmed"} objects.
#'
#' @examples
#' \donttest{
#' # Count outcome
#' eval_1 <- catecv(response = "count",
#'                  data = countExample,
#'                  score.method = "poisson",
#'                  cate.model = y ~ age + female + previous_treatment +
#'                                   previous_cost + previous_number_relapses + offset(log(years)),
#'                  ps.model = trt ~ age + previous_treatment,
#'                  higher.y = FALSE,
#'                  cv.n = 5)
#'
#' # default setting
#' plot(eval_1)
#'
#' # turn off ABC annotation
#' plot(eval_1, show.abc = FALSE)
#'
#' # use a different theme
#' plot(eval_1, theme = ggplot2::theme_bw())
#'
#' # plot the validation curves from the 2nd CV iteration instead of the mean
#' # of all validation curves
#' plot(eval_1, cv.i = 2)
#'
#' # median of the validation curves
#' plot(eval_1, combine = "median")
#'
#' # plot validation curves in validation set only
#' plot(eval_1, valid.only = TRUE)
#'
#' # Survival outcome
#' library(survival)
#' tau0 <- with(survivalExample,
#'              min(quantile(y[trt == "drug1"], 0.95), quantile(y[trt == "drug0"], 0.95)))
#' eval_2 <- catecv(response = "survival",
#'                  data = survivalExample,
#'                  score.method = c("poisson", "randomForest"),
#'                  cate.model = Surv(y, d) ~ age + female + previous_cost +
#'                                            previous_number_relapses,
#'                  ps.model = trt ~ age + previous_treatment,
#'                  initial.predictor.method = "randomForest",
#'                  ipcw.model = ~ age + previous_cost + previous_treatment,
#'                  tau0 = tau0,
#'                  cv.n = 5,
#'                  seed = 999)
#'
#'
#' # default setting, plot RMTL ratios in both training and validation sets
#' plot(eval_2)
#'
#' # plot hazard ratio
#' plot(eval_2, plot.hr = TRUE)
#'
#'}
#'
#' @export
#'
#' @importFrom rlang .data
#' @importFrom stringr str_extract
#' @importFrom dplyr filter select mutate group_by ungroup
#' @importFrom ggplot2 ggplot aes facet_wrap geom_boxplot geom_hline geom_line
#' labs scale_color_manual scale_linetype_manual scale_y_continuous
#' theme theme_classic waiver
#' @importFrom tidyr fill

plot.precmed <- function(x,
                         cv.i = NULL,
                         combine = "mean",
                         show.abc = TRUE,
                         valid.only = FALSE,
                         plot.hr = FALSE,
                         ylab = NULL,
                         legend.position = "bottom",
                         xlim = NULL,
                         title = waiver(),
                         theme = theme_classic(), ...) {

  # Must be count or survival data for this function
  stopifnot(x$response %in% c("count", "survival", "continuous"))

  # Check the value of show.abc, combine, valid.only, plot.hr
  if (!(show.abc %in% c(TRUE, FALSE))) stop("show.abc has to be boolean.")
  if (!(combine %in% c("mean", "median"))) stop("combine should be: 'mean' or 'median'.")
  if (!(valid.only %in% c(TRUE, FALSE))) stop("valid.only has to be boolean.")

  # Check the value of cv.i
  if (is.null(cv.i) == FALSE && (cv.i > x$cv.n)) {
    stop("cv.i provided is too large. Must be <= the number of cross-validation folds x$cv.n.")
  } else if (is.null(cv.i) == FALSE && (cv.i %% 1 != 0)) {
    stop("cv.i must be an integer or NULL.")
  }

  # Check plot.hr and only available for survival data
  if (plot.hr == TRUE & x$response != "survival") stop("Hazard ratio plot is only available for precmed object with x$response equals to 'survival'.")
  if (!(plot.hr %in% c(TRUE, FALSE))) stop("plot.hr has to be boolean.")

  # Retrieve proportion
  prop.onlyhigh <- round(x$props$prop.onlyhigh, 5)
  if (is.null(xlim) == TRUE) xlim <- c(min(prop.onlyhigh), max(prop.onlyhigh))
  xlim <- round(xlim, 5)
  xmin.idx <- min(which(prop.onlyhigh >= xlim[1]))
  xmax.idx <- max(which(prop.onlyhigh <= xlim[2]))
  l <- length(prop.onlyhigh)

  # Argument checks for HR with survival outcomes
  if (show.abc == TRUE) {
    # Must have abc results if want to show abc
    if (x$abc == FALSE) {
      warning("ABC will not be shown because x$abc is FALSE.")
    } else if (plot.hr == TRUE) {
      warning("ABC will not be shown in hazard ratio plot.")
    } else if (is.null(xlim) == FALSE & (xlim[1] != min(prop.onlyhigh)) | (xlim[2] != max(prop.onlyhigh))) {
      warning("ABC will not be shown when xlim does not cover the entire prop.cutoff range.")
    }
  }

  # Get name of treatment variable
  trt <- all.vars(x$formulas$ps.model)[1]

  # Define y-axis if default is NULL
  if (is.null(ylab)) {

    if (x$response == "count") {
      plot.ratio <- "Rate ratio"
    } else if (x$response == "survival") {
      if (plot.hr) {
        plot.ratio <- "Hazard ratio"
      } else {
        plot.ratio <- "RMTL ratio"
      }
    } else if (x$response == "continuous") {
      plot.ratio <- "Mean difference"
    }

    if (is.null(x$formulas$trt_labels) & x$response == "continuous") {
      ylab <- paste0(plot.ratio, " ", trt, "=1 - ", trt, "=0\nin each subgroup")

    } else if (is.null(x$formulas$trt_labels) & x$response != "continuous") {
      ylab <- paste0(plot.ratio, " ", trt, "=1 - ", trt, "=0\nin each subgroup")
    } else {
      ylab <- paste0(plot.ratio, " ", trt, "=", x$formulas$trt_labels[2], " - ", trt, "=", x$formulas$trt_labels[1], "\nin each subgroup")
    }
  }

  # Define x-axis label
  if ((x$higher.y == TRUE & x$response == "survival") | (x$higher.y == FALSE & x$response == "count") | (x$higher.y == FALSE & x$response == "continuous")) {
    xlab <- "Subgroup size\n(proportion of subjects with lowest estimated CATE)"
  } else {
    xlab <- "Subgroup size\n(proportion of subjects with highest estimated CATE)"
  }

  # Detect the number of score methods
  map <- c("randomForest" = "Random Forest",
           "boosting" = "Boosting",
           "poisson" = "Naive Poisson",
           "twoReg" = "Two Regressions",
           "contrastReg" = "Contrast Regression",
           "negBin" = "Negative Binomial",
           "gaussian" = "Linear Regression",
           "gam" = "Generalized Additive Models")
  score.method <- str_extract(names(x), "(?<=^ate\\.).*$")
  score.method <- score.method[!is.na(score.method)]
  m <- length(score.method)

  # Extract overall ATE
  if (plot.hr == FALSE | (x$response == "count" | x$response == "continuous")) {
    if (is.null(cv.i) == TRUE & combine == "mean") {
      overall.hline.train <- overall.hline.train <- mean(x$overall.ate.train)
      overall.hline.valid <- mean(x$overall.ate.valid)
    } else if (is.null(cv.i) == TRUE & combine == "median") {
      overall.hline.train <- median(x$overall.ate.train)
      overall.hline.valid <- median(x$overall.ate.valid)
    } else { # if cv.i = integer
      overall.hline.train <- x$overall.ate.train[cv.i]
      overall.hline.valid <- x$overall.ate.valid[cv.i]
    }
  } else { # if we look at hr
    if (is.null(cv.i) == TRUE & combine == "mean") {
      overall.hline.train <- mean(x$overall.hr.train)
      overall.hline.valid <- mean(x$overall.hr.valid)
    } else if (is.null(cv.i) == TRUE & combine == "median") {
      overall.hline.train <- median(x$overall.hr.train)
      overall.hline.valid <- median(x$overall.hr.valid)
    } else { # if cv.i = integer
      overall.hline.train <- x$overall.hr.train[cv.i]
      overall.hline.valid <- x$overall.hr.valid[cv.i]
    }
  }

  # Retrieve results needed to generate the plot
  results <- rep(NA, m * l * 2)
  abc <- rep(NA, m * l * 2)

  if (plot.hr == FALSE | (x$response == "count" | x$response == "continuous")) {
    prefix <- "ate."
  } else {
    prefix <- "hr."
    show.abc <- FALSE
  }

  for (i in seq(m)) {
    ate <- x[[paste0(prefix, score.method[i])]]
    if (is.null(cv.i) == TRUE) {
      if (combine == "mean") {
        # Plot the average of all CVs
        results[((i - 1) * l + 1):(i * l)] <- rowMeans(ate[[paste0(prefix, "est.train.high.cv")]], na.rm = TRUE)
        results[((i - 1) * l + 1 + m * l):(i * l + m * l)] <- temp <- rowMeans(ate[[paste0(prefix, "est.valid.high.cv")]], na.rm = TRUE)
      } else {
        # Plot the median of all CVs
        results[((i - 1) * l + 1):(i * l)] <- apply(ate[[paste0(prefix, "est.train.high.cv")]], 1, median, na.rm = TRUE)
        results[((i - 1) * l + 1 + m * l):(i * l + m * l)] <- temp <- apply(ate[[paste0(prefix, "est.valid.high.cv")]], 1, median, na.rm = TRUE)
      }

      # Calculate ABC based on average validation curve
      # abc[m * l + (i - 1) * l + 1] <- abc(x)
      if (show.abc == TRUE & (xlim[1] == min(prop.onlyhigh)) & (xlim[2] == max(prop.onlyhigh))) {
        if ((x$higher.y == FALSE & x$response == "survival") | (x$higher.y == TRUE & x$response == "count")) {
          abc[m * l + (i - 1) * l + 1] <- auc(x = prop.onlyhigh, y = log(temp) - log(overall.hline.valid),
                                              from = prop.onlyhigh[xmin.idx], to = prop.onlyhigh[xmax.idx],
                                              type = "spline")

        } else if ((x$higher.y == TRUE & x$response == "survival") | (x$higher.y == FALSE & x$response == "count")){
          abc[m * l + (i - 1) * l + 1] <- auc(x = prop.onlyhigh, y = log(overall.hline.valid) - log(temp),
                                              from = prop.onlyhigh[xmin.idx], to = prop.onlyhigh[xmax.idx],
                                              type = "spline")

        } else if (x$higher.y == TRUE & x$response == "continuous") {
          abc[m * l + (i - 1) * l + 1] <- auc(x = prop.onlyhigh, y = temp - overall.hline.valid,
                                              from = prop.onlyhigh[xmin.idx], to = prop.onlyhigh[xmax.idx],
                                              type = "spline")

        } else if (x$higher.y == FALSE & x$response == "continuous"){
          abc[m * l + (i - 1) * l + 1] <- auc(x = prop.onlyhigh, y = overall.hline.valid - temp,
                                              from = prop.onlyhigh[xmin.idx], to = prop.onlyhigh[xmax.idx],
                                              type = "spline")
        }
      }
    } else {
      # Plot each individual cv.i
      results[((i - 1) * l + 1):(i * l)] <- ate[[paste0(prefix, "est.train.high.cv")]][, cv.i]
      results[((i - 1) * l + 1 + m * l):(i * l + m * l)] <- ate[[paste0(prefix, "est.valid.high.cv")]][, cv.i]

      # Calculate ABC based on validation curve
      if (show.abc == TRUE) abc[m * l + (i - 1) * l + 1] <- ate$abc.valid[cv.i]
    }
  }

  # Gather data to be plotted
  mapped.method <- map[score.method]
  if (x$response == "count") abc.gray.value <- ifelse(sprintf('%.3f', abc) == "NA", "", paste(as.character(mapped.method), "ABC:", sprintf('%.3f', abc)))
  if (x$response == "survival") abc.gray.value <- ifelse(sprintf('%.3f', abc) == "NA", "", paste(rep(rep(mapped.method, each = l), 2), "ABC:", sprintf('%.3f', abc)))
  if (x$response == "continuous") abc.gray.value <- ifelse(sprintf('%.3f', abc) == "NA", "", paste(as.character(mapped.method), "ABC:", sprintf('%.3f', abc)))

  dat <- data.frame(ate.cv = results,
                    set = factor(c(rep(paste0("Training set CV", cv.i), l * m),
                                   rep(paste0("Validation set CV", cv.i), l * m)),
                                 levels = c(paste0("Training set CV", cv.i),
                                            paste0("Validation set CV", cv.i))),
                    prop = rep(prop.onlyhigh, times = 2 * m),
                    score = factor(rep(rep(mapped.method, each = l), 2),
                                   levels = mapped.method)) %>%
    mutate(abc.gray = abc.gray.value,
           abc = ifelse(sprintf('%.3f', abc) == "NA", "", paste("ABC:", sprintf('%.3f', abc))),
           score.abc = ifelse(abc == "", NA, paste0(.data$score, "\nValidation ", abc))) %>%
    group_by(.data$score) %>%
    tidyr::fill(.data$score.abc, .direction = "downup") %>%
    ungroup() %>%
    filter(.data$prop >= xlim[1] & .data$prop <= xlim[2])

  if (valid.only == TRUE) {
    dat <- dat %>% filter(.data$set == paste0("Validation set CV", cv.i))
  }

  lowery <- min(c(dat$ate.cv, overall.hline.train, overall.hline.valid))
  uppery <- max(c(dat$ate.cv, overall.hline.train, overall.hline.valid))

  if (show.abc == TRUE & (xlim[1] == min(prop.onlyhigh)) & (xlim[2] == max(prop.onlyhigh))) {
    dat$score <- factor(dat$score.abc, levels = unique(dat$score.abc))
  }
  scorevars <- unique(dat$score)

  map.color <- c("randomForest" = "magenta4",
                 "boosting" = "orangered2",
                 "poisson" = "gray30",
                 "gaussian" = "gray30", #Only one of poisson or gaussian is mapped depending on the data type
                 "twoReg" = "springgreen3",
                 "contrastReg" = "dodgerblue",
                 "negBin" = "magenta4",
                 "gam" = "red4")

  colors <- map.color[score.method]
  attr(colors, "names") <- sapply(1:m, function(x) grep(mapped.method[x], scorevars, value = TRUE))

  p <- dat %>%
    ggplot(aes(x = .data$prop, y = .data$ate.cv, color = .data$score)) +
    geom_line(size = 1.2) +
    facet_wrap(~ .data$set, ncol = 2) +
    scale_y_continuous(limits = c(lowery, uppery)) +
    labs(y = ylab, x = xlab, title = title) +
    theme +
    theme(legend.position = legend.position) +
    scale_color_manual(name = "Method", values = colors)

  # Add overall ATE
  abline <- data.frame(yintercept = c(overall.hline.train, overall.hline.valid),
                       set = factor(c(paste0("Training set CV", cv.i),
                                      paste0("Validation set CV", cv.i))))
  if (valid.only == TRUE) {
    abline <- abline %>% filter(.data$set == paste0("Validation set CV", cv.i))
  }

  p <- p + geom_hline(data = abline, aes(yintercept = .data$yintercept),
                      linetype = 2, color = "gray90", size = 0.8)

  return(p)
}


#' A set of box plots of estimated ATEs from the \code{"precmed"} object
#'
#' Provides box plots which depict distributions of estimated ATEs for each multi-category subgroup in
#' the validation set across all cross-validation iterations. The subgroups are mutually exclusive and
#' are categorized by the CATE score percentiles (\code{prop.multi} specified in \code{\link{catecv}()} or
#' \code{\link{catecvmean}()}). Box plots of mutually exclusive subgroups are constructed separately by scoring
#' method specified in \code{\link{catecv}()}. This should be run only after results of \code{\link{catecv}()} or
#' \code{\link{catecvmean}()}) have been obtained.
#'
#' @param x An object of class \code{"precmed"}.
#' @param ylab A character value for the y-axis label to describe what the ATE is. Default is \code{NULL},
#' which creates a default y-axis label based on available data.
#' @param plot.hr A logical value indicating whether the hazard ratios should be plotted in the
#' validation curves (\code{TRUE}). Otherwise, the restricted mean time lost is plotted (\code{FALSE}).
#' This argument is only applicable to survival outcomes. Default is \code{FALSE}.
#' @param title The text for the title
#' @param theme Defaults to \code{theme_classic()}. Other options include \code{theme_grey()}, \code{theme_bw()}, \code{theme_light()}, \code{theme_dark()}, and \code{theme_void()}
#' @param ... Other parameters
#'
#' @return Returns sets of box plots, one set for each scoring method, over each of the multi-category
#' subgroups. A gray horizontal dashed line of the overall ATE is included as a reference.
#'
#' @details \code{\link{boxplot}()} takes in outputs from \code{\link{catecv}()} and generates
#' the box plots of estimated ATEs for multi-category subgroups of the validation set.
#' The box plots together with the overall ATE reference line can help compare the scoring methods'
#' ability to distinguish subgroups of patients with different treatment effects.
#'
#' For a given scoring method, box plots showing increasing or decreasing trends across the
#' multi-category subgroups indicate presence of treatment effect heterogeneity
#' (and the ability of the scoring method to capture it). On the contrary, box plots which
#' are relatively aligned across the multi-category subgroups indicate absence of treatment
#' effect heterogeneity (or the inability of the scoring method to capture it).
#'
#' @references Yadlowsky, S., Pellegrini, F., Lionetto, F., Braune, S., & Tian, L. (2020).
#' \emph{Estimation and validation of ratio-based conditional average treatment effects using
#' observational data. Journal of the American Statistical Association, 1-18.}
#' \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1772080}
#'
#' @seealso \code{\link{plot}} and \code{\link{abc}()} for \code{"precmed"} objects.
#'
#' @examples
#' \donttest{
#' # Count outcome
#' eval_1 <- catecv(response = "count",
#'                  data = countExample,
#'                  score.method = "poisson",
#'                  cate.model = y ~ age + female + previous_treatment +
#'                                   previous_cost + previous_number_relapses +
#'                                   offset(log(years)),
#'                  ps.model = trt ~ age + previous_treatment,
#'                  higher.y = FALSE,
#'                  cv.n = 5)
#'
#' boxplot(eval_1, ylab = "Rate ratio of drug1 vs drug0 in each subgroup")
#'
#' # Survival outcome
#' library(survival)
#' tau0 <- with(survivalExample,
#'              min(quantile(y[trt == "drug1"], 0.95), quantile(y[trt == "drug0"], 0.95)))
#' eval_2 <- catecv(response = "survival",
#'                  data = survivalExample,
#'                  score.method = c("poisson", "randomForest"),
#'                  cate.model = Surv(y, d) ~ age + female + previous_cost +
#'                                            previous_number_relapses,
#'                  ps.model = trt ~ age + previous_treatment,
#'                  initial.predictor.method = "randomForest",
#'                  ipcw.model = ~ age + previous_cost + previous_treatment,
#'                  tau0 = tau0,
#'                  higher.y = TRUE,
#'                  cv.n = 5,
#'                  seed = 999)
#'
#' boxplot(eval_2, ylab = "RMTL ratio of drug1 vs drug0 in each subgroup")
#'}
#'
#' @export
#'
#' @importFrom graphics boxplot
#' @importFrom dplyr mutate_at vars contains
#' @importFrom ggplot2 ggplot aes facet_wrap geom_boxplot geom_hline geom_line
#' labs scale_color_manual scale_linetype_manual scale_y_continuous theme
#' theme_classic waiver scale_fill_hue
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data


boxplot.precmed <- function(x, ylab = NULL,
                            plot.hr = FALSE,
                            title = waiver(),
                            theme = theme_classic(), ...) {

  # Must be count or survival data for this function
  stopifnot(x$response %in% c("count", "survival", "continuous"))

  # Argument checks on plot.hr
  if (!(plot.hr %in% c(TRUE, FALSE))) stop("plot.hr has to be boolean.")
  if (plot.hr & x$response != "survival") stop("Hazard ratio plot is only available for precmed objects with x$response is \"survival\".")

  # Grab number of CV iterations
  cv.n <- ncol(x[[1]]$ate.est.valid.group.cv)

  # Retrieve proportion
  prop.multi <- x$props$prop.multi
  l <- length(prop.multi) - 1

  # Define x-axis
  if ((x$higher.y == TRUE & x$response == "survival") | (x$higher.y == FALSE & x$response == "count") | (x$higher.y == FALSE & x$response == "continuous")) {
    xlab <- paste0("Subgroup \n (from lowest [0%] to highest [100%] estimated CATE)") # lower score = better responder to trt=1
  } else {
    xlab <- paste0("Subgroup \n (from highest [0%] to lowest [100%] estimated CATE)") # higher score = better responder to trt=1
  }

  # Get name of treatment variable
  trt <- all.vars(x$formulas$ps.model)[1]

  # Define y-axis if default is NULL
  if (is.null(ylab)) {
    if (x$response == "count") {
      plot.ratio = "Rate ratio"
    } else if (x$response == "survival") {
      if (plot.hr) {
        plot.ratio = "Hazard ratio"
      } else {
        plot.ratio = "RMTL ratio"
      }
    } else if (x$response == "continuous") {
      plot.ratio = "Mean difference"
    }

    if (is.null(x$formulas$trt_labels) == TRUE & x$response == "continuous") {
      ylab <- paste0(plot.ratio, " ", trt, "=1 - ", trt, "=0\nin each subgroup")

    } else if (is.null(x$formulas$trt_labels) == TRUE & x$response != "continuous") {
      ylab <- paste0(plot.ratio, " ", trt, "=1 - ", trt, "=0\nin each subgroup")
    } else {
      ylab <- paste0(plot.ratio, " ", trt, "=", x$formulas$trt_labels[2], " - ", trt, "=", x$formulas$trt_labels[1], "\nin each subgroup")
    }
  }

  # Detect the number of score methods
  map <- c("randomForest" = "Random Forest",
           "boosting" = "Boosting",
           "poisson" = "Naive Poisson",
           "gaussian" = "Linear Regression",
           "negBin" = "Negative Binomial",
           "twoReg" = "Two Regressions",
           "contrastReg" = "Contrast Regression",
           "gam" = "Generalized Additive Models")
  score.method <- str_extract(names(x), "(?<=^ate\\.).*$")
  score.method <- score.method[!is.na(score.method)]
  m <- length(score.method)
  mapped.method <- map[score.method]

  # Retrieve results needed to generate the plot
  results <- matrix(NA, nrow = cv.n * m, ncol = l + 2)
  prefix <- ifelse(plot.hr, "hr.", "ate.")

  for (i in seq(m)) {
    ate <- x[[paste0(prefix, score.method[i])]]
    if (cv.n == 1) {
      results[((i - 1) * cv.n + 1):((i - 1) * cv.n + cv.n),] <- c(
        t(ate[[paste0(prefix, "est.valid.group.cv")]])[l:1],
        1,
        unname(mapped.method[i]))
    } else {
      results[((i - 1) * cv.n + 1):((i - 1) * cv.n + cv.n),] <- cbind(
        t(ate[[paste0(prefix, "est.valid.group.cv")]])[,l:1],
        1:cv.n,
        unname(mapped.method[i]))
    }
  }

  # Extract overall ATE
  overall.hline <- ifelse(plot.hr, mean(x$overall.hr.valid), mean(x$overall.ate.valid))


  results <- as.data.frame(results)
  colnames(results) <- c(1:l, "CV", "Method") # Subgroup 1 = lowest responder to drug1
  colnames(results)[1:l] <- subgroup_label <- paste(round(prop.multi[1:l] * 100), paste0(round(prop.multi[2:(l + 1)] * 100), "%"), sep = "-")


  p <- results %>%
    mutate_at(c(1:l), as.numeric) %>%
    pivot_longer(-(.data$CV:.data$Method), names_to = "Subgroup", values_to = "ate") %>%
    mutate(Method = factor(.data$Method, levels = mapped.method),
           Subgroup = factor(.data$Subgroup, levels = subgroup_label)) %>%
    ggplot(aes(x = .data$Subgroup, y = .data$ate)) +
    geom_boxplot(aes(fill = .data$Subgroup)) +
    geom_hline(yintercept = overall.hline, linetype = "dashed", color = "grey70") +
    facet_wrap(~ .data$Method, nrow = 2, scales = "fixed") +
    scale_fill_hue(l = 45) +
    labs(y = ylab, x = xlab, title = title) +
    theme + theme(legend.position = "none")

  return(p)
}


#' Split the given dataset into balanced training and validation sets
#' (within a pre-specified tolerance)
#' Balanced means 1) The ratio of treated and controls is maintained in the training and validation sets
#'                2) The covariate distributions are balanced between the training and validation sets
#'
#' @param y Observed outcome; vector of size \code{n} (observations)
#' @param trt Treatment received; vector of size \code{n} with treatment coded as 0/1
#' @param x.cate Matrix of \code{p.cate} baseline covariates; dimension \code{n} by \code{p.cate}
#' (covariates in the outcome model)
#' @param x.ps Matrix of \code{p.ps} baseline covariates (plus a leading column of 1 for the intercept);
#' dimension \code{n} by \code{p.ps + 1} (covariates in the propensity score model plus intercept)
#' @param time Log-transformed person-years of follow-up; vector of size \code{n}
#' @param minPS A numerical value (in `[0, 1]`) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in `(0, 1]`) above which estimated propensity scores should be
#' truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.
#' @param train.prop A numerical value (in `(0, 1)`) indicating the proportion of total data used
#' for training. Default is \code{3/4}.
#' @param error.max A numerical value > 0 indicating the tolerance (maximum value of error)
#' for the largest standardized absolute difference in the covariate distributions or in the
#' doubly robust estimated rate ratios between the training and validation sets. This is used
#' to define a balanced training-validation splitting. Default is \code{0.1}.
#' @param max.iter A positive integer value indicating the maximum number of iterations when
#' searching for a balanced training-validation split. Default is \code{5,000}.
#'
#' @return A list of 10 objects, 5 training and 5 validation of y, trt, x.cate, x.ps, time:
#'             y.train          - observed outcome in the training set; vector of size \code{m} (observations in the training set)
#'             trt.train        - treatment received in the training set; vector of size \code{m} coded as 0/1
#'             x.cate.train     - baseline covariates for the outcome model in the training set; matrix of dimension \code{m} by \code{p.cate}
#'             x.ps.train       - baseline covariates (plus intercept) for the propensity score model in the training set; matrix of dimension \code{m} by \code{p.ps + 1}
#'             time.train       - log-transformed person-years of follow-up in the training set; vector of size \code{m}
#'             y.valid          - observed outcome in the validation set; vector of size \code{n-m}
#'             trt.valid        - treatment received in the validation set; vector of size \code{n-m} coded as 0/1
#'             x.cate.valid     - baseline covariates for the outcome model in the validation set; matrix of dimension \code{n-m} by \code{p.cate}
#'             x.ps.valid       - baseline covariates (plus intercept) for the propensity score model in the validation set; matrix of dimension \code{n-m} by \code{p.ps + 1}
#'             time.valid       - log-transformed person-years of follow-up in the validation set; vector of size \code{n-m}
balance.split <- function(y, trt, x.cate, x.ps, time,
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
    time.valid <- time[id.valid]

    y.train <- y[-id.valid]
    trt.train <- trt[-id.valid]
    x.train <- x[-id.valid, , drop = FALSE]
    x.cate.train <- x.cate[-id.valid, , drop = FALSE]
    x.ps.train <- x.ps[-id.valid, , drop = FALSE]
    time.train <- time[-id.valid]

    diffx <- colMeans(x.train) - colMeans(x.valid)
    errorx <- max(abs(diffx / sdx))

    est.valid  <- drcount(y = y.valid, x.cate = x.cate.valid, x.ps = x.ps.valid, trt = trt.valid, time = time.valid,
                          ps.method = "glm", minPS = minPS, maxPS = maxPS, interactions = TRUE)$log.rate.ratio
    est.train <- drcount(y = y.train, x.cate = x.train, x.ps = x.ps.train, trt = trt.train, time = time.train,
                         ps.method = "glm", minPS = minPS, maxPS = maxPS, interactions = TRUE)$log.rate.ratio

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
    time.valid <- time[bestid.valid]

    y.train <- y[-bestid.valid]
    trt.train <- trt[-bestid.valid]
    x.cate.train <- x.cate[-bestid.valid, , drop = FALSE]
    x.ps.train <- x.ps[-bestid.valid, , drop = FALSE]
    time.train <- time[-bestid.valid]
    warning(paste("Maximum iteration reached and the SMD between training and validation set is still greater than error.max (error=", round(errorbest, 4), "). Consider increasing max.iter, decreasing error.max, or increasing sample size.", sep = ""))
  }

  return(list(y.train = y.train, trt.train = trt.train, x.cate.train = x.cate.train, x.ps.train = x.ps.train, time.train = time.train,
              y.valid = y.valid, trt.valid = trt.valid, x.cate.valid = x.cate.valid, x.ps.valid = x.ps.valid, time.valid = time.valid))
}


