# ------------------------------------------------------------------
#
# Project: Precision Medicine MS (precmed) - Comprehensive R package
#
# Purpose: Main output functions for precmed package
#
# Platform: Windows
# R Version: 4.1.0
#



#' Cross-validation of the conditional average treatment effect (CATE) score for count, survival or continuous outcomes
#'
#' Provides (doubly robust) estimation of the average treatment effect (ATE) for count, survival or continuous
#' outcomes in nested and mutually exclusive subgroups of patients defined by an estimated conditional
#' average treatment effect (CATE) score via cross-validation (CV).
#'
#' @param response A string describing the type of outcome in the data. Allowed values include
#' "count" (see \code{\link{cvcount}()}), "survival" (see \code{\link{cvsurv}()}) and "continuous" (see \code{\link{cvmean}()}) .
#' @param cate.model A formula describing the outcome model to be fitted.
#' The outcome must appear on the left-hand side. For survival outcomes, a \code{Surv} object
#' must be used to describe the outcome.
#' @param init.model A formula describing the initial predictor model. The outcome must appear on the left-hand side.
#' It must be specified when \code{score.method = contrastReg} or \code{twoReg}.
#' @param ps.model A formula describing the propensity score (PS) model to be fitted. The treatment must
#' appear on the left-hand side. The treatment must be a numeric vector coded as 0/1.
#' If data are from a randomized controlled trial, specify \code{ps.model = ~1} as an intercept-only model.
#' @param data A data frame containing the variables in the outcome, propensity score, and inverse
#' probability of censoring models (if specified); a data frame with \code{n} rows (1 row per observation).
#' @param score.method A vector of one or multiple methods to estimate the CATE score.
#' Allowed values are: \code{'boosting'}, \code{'twoReg'}, \code{'contrastReg'}, \code{'poisson'} (count and survival outcomes only),
#' \code{'randomForest'} (survival, continuous outcomes only), \code{negBin} (count outcomes only), \code{'gam'} (continuous outcomes only),
#' \code{'gaussian'} (continuous outcomes only).
#' @param ipcw.model A formula describing the inverse probability of censoring weighting (IPCW)
#' model to be fitted. The left-hand side must be empty. Only applies for survival outcomes.
#' Default is \code{NULL}, which corresponds to specifying the IPCW with the same covariates
#' as the outcome model \code{cate.model}, plus the treatment.
#' @param followup.time A column name in \code{data} specifying the maximum follow-up time,
#' interpreted as the potential censoring time. Only applies for survival outcomes.
#' Default is \code{NULL}, which corresponds to unknown potential censoring time.
#' @param tau0 The truncation time for defining restricted mean time lost. Only applies for
#' survival outcomes. Default is \code{NULL}, which corresponds to setting the truncation time as the
#' maximum survival time in the data.
#' @param surv.min Lower truncation limit for the probability of being censored.
#' It must be a positive value and should be chosen close to 0. Only applies for survival
#' outcomes. Default is \code{0.025}.
#' @param ipcw.method A character value for the censoring model. Only applies for survival
#' outcomes. Allowed values are: \code{'breslow'} (Cox regression with Breslow estimator of t
#' he baseline survivor function), \code{'aft (exponential)'}, \code{'aft (weibull)'},
#' \code{'aft (lognormal)'} or \code{'aft (loglogistic)'} (accelerated failure time model
#' with different distributions for y variable). Default is \code{'breslow'}.
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
#' outcome predictions conditional on the covariates specified in \code{cate.model}. Only applies
#' when \code{score.method} includes \code{'twoReg'} or \code{'contrastReg'}. Allowed values include one of
#' \code{'randomForest'} (survival outcomes only), \code{'boosting'}, \code{'logistic'}
#' (survival outcomes only, fast), \code{'poisson'} (count outcomes only, fast), \code{'gaussian'} (continuous outcomes only) and
#' \code{'gam'} (count and continuous outcomes only). Default is \code{NULL}, which assigns \code{'boosting'}
#' for count outcomes and \code{'randomForest'} for survival outcomes.
#' @param xvar.smooth.init A vector of characters indicating the name of the variables used as
#' the smooth terms if \code{initial.predictor.method = 'gam'}. The variables must be selected
#' from the variables listed in \code{init.model}.
#' Default is \code{NULL}, which uses all variables in \code{init.model}.
#' @param xvar.smooth.score A vector of characters indicating the name of the variables used as
#' the smooth terms if \code{score.method = 'gam'}. The variables must be selected
#' from the variables listed in \code{cate.model}.
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
#' \code{2} means progress bar, run time, and all errors and warnings. Default is \code{2}.
#'
#' @return For count response, see description of outputs in \code{\link{cvcount}()}.
#' For survival response, see description of outputs in \code{\link{cvsurv}()}.
#' For continuous response, see description of outputs in \code{\link{cvmean}()}.
#'
#' @details For count response, see details in \code{\link{cvcount}()}.
#' For survival response, see details in \code{\link{cvsurv}()}.
#' For continuous response, see details in \code{\link{cvmean}()}.
#'
#' @references Yadlowsky, S., Pellegrini, F., Lionetto, F., Braune, S., & Tian, L. (2020).
#' \emph{Estimation and validation of ratio-based conditional average treatment effects using
#' observational data. Journal of the American Statistical Association, 1-18.}
#' \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1772080}
#'
#' @seealso \code{\link{pm}()} function and \code{\link{boxplot}()}, \code{\link{abc}} methods for
#' \code{"precmed"} objects.
#'
#' @examples
#' \dontrun{
#' output_cv <- cv(response = "count",
#'                 cate.model = y ~ age + female + previous_treatment +
#'                                      previous_cost + previous_number_relapses + offset(log(years)),
#'                 ps.model = trt ~ age + previous_treatment,
#'                 data = countExample,
#'                 higher.y = FALSE,
#'                 score.method = "poisson",
#'                 cv.n = 5,
#'                 plot.gbmperf = FALSE,
#'                 seed = 999)
#'
#' # Try:
#' plot(x = output_cv, ylab = "RMTL ratio of drug1 vs drug0 in each subgroup")
#' boxplot(x = output_cv, ylab = "RMTL ratio of drug1 vs drug0 in each subgroup")
#' abc(x = output_cv)
#'
#' # Survival outcome
#' tau0 <- with(survivalExample,
#'              min(quantile(y[trt == "drug1"], 0.95), quantile(y[trt == "drug0"], 0.95)))
#'
#' output_cv2 <- cv(response = "survival",
#'                  cate.model = survival::Surv(y, d) ~ age +
#'                                                      female +
#'                                                      previous_cost +
#'                                                      previous_number_relapses,
#'                  ps.model = trt ~ age + previous_treatment,
#'                  ipcw.model = ~ age + previous_cost + previous_treatment,
#'                  data = survivalExample,
#'                  score.method = c("poisson", "randomForest"),
#'                  followup.time = NULL,
#'                  tau0 = tau0,
#'                  surv.min = 0.025,
#'                  higher.y = TRUE,
#'                  cv.n = 5,
#'                  initial.predictor.method = "randomForest",
#'                  plot.gbmperf = FALSE,
#'                  seed = 999)
#'
#' # Try:
#' plot(x = output_cv2, ylab = "RMTL ratio of drug1 vs drug0 in each subgroup")
#' boxplot(x = output_cv2, ylab = "RMTL ratio of drug1 vs drug0 in each subgroup")
#' abc(x = output_cv2)
#'
#'
#' # Continuous outcome
#' output_cv3 <- cv(response = "continuous",
#'                 cate.model = y ~ age +
#'                                  previous_treatment +
#'                                  previous_cost +
#'                                  previous_status_measure,
#'                 init.model = y ~ age +
#'                                  previous_treatment +
#'                                  previous_cost +
#'                                  previous_status_measure,
#'                 ps.model = trt ~ previous_status_measure,
#'                 data = meanExample,
#'                 higher.y = FALSE,
#'                 score.method = c("gaussian", "twoReg"),
#'                 xvar.smooth.score = c("age", "previous_cost"),
#'                 initial.predictor.method = "gaussian",
#'                 cv.n = 5,
#'                 plot.gbmperf = FALSE,
#'                 seed = 999)
#'
#' # Try:
#' plot(x = output_cv3)
#' boxplot(x = output_cv3)
#' abc(x = output_cv3)
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
#' @importFrom MESS auc
#' @importFrom randomForestSRC rfsrc predict.rfsrc
#' @importFrom stats as.formula coef glm median model.frame model.matrix model.offset model.response na.omit optim pchisq predict qnorm quantile sd var
#' @importFrom stringr str_replace str_extract str_detect
#' @importFrom survival Surv coxph coxph.detail survreg
#' @importFrom utils setTxtProgressBar txtProgressBar

cv <- function(response, cate.model, ps.model, data, score.method,  # Mandatory arguments (count & survival & continuous)
               init.model = NULL, # This is used in continuous data, when contrast/two regression is used
               ipcw.model = NULL, followup.time = NULL, tau0 = NULL, # Non-mandatory arguments survival only
               surv.min = 0.025, ipcw.method = "breslow",
               higher.y = TRUE, abc = TRUE, # Non-mandatory arguments survival & count & continuous (except xvar.smooth)
               prop.cutoff = seq(0.5, 1, length = 6), prop.multi = c(0, 1/3, 2/3, 1),
               ps.method = "glm", minPS = 0.01, maxPS = 0.99,
               train.prop = 3/4, cv.n = 10, error.max = 0.1, max.iter = 5000,
               initial.predictor.method = NULL, xvar.smooth.score = NULL, xvar.smooth.init = NULL,
               tree.depth = 2, n.trees.rf = 1000, n.trees.boosting = 200, B = 3, Kfold = 5,
               error.maxNR = 1e-3, max.iterNR = 150, tune = c(0.5, 2),
               seed = NULL, plot.gbmperf = TRUE, verbose = 2) {

  stopifnot("`response` must be either `count`, `survival`, or `continuous`." = any(response== c("count", "survival", "continuous")))

  if (response == "count"){
    cvout <- cvcount(cate.model = cate.model, ps.model = ps.model, data = data, score.method = score.method,
                     higher.y = higher.y, abc = abc,
                     prop.cutoff = prop.cutoff, prop.multi = prop.multi,
                     ps.method = ps.method, minPS = minPS, maxPS = maxPS,
                     train.prop = train.prop, cv.n = cv.n, error.max = error.max, max.iter = max.iter,
                     initial.predictor.method = initial.predictor.method, xvar.smooth = xvar.smooth.score,
                     tree.depth = tree.depth, n.trees.boosting = n.trees.boosting, B = B, Kfold = Kfold,
                     error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune,
                     seed = seed, plot.gbmperf = plot.gbmperf, verbose = verbose)
  }
  if (response == "survival"){
    cvout <- cvsurv(cate.model = cate.model, ps.model = ps.model, data = data, score.method = score.method,
                    ipcw.model = ipcw.model, followup.time = followup.time, tau0 = tau0,
                    surv.min = surv.min, ipcw.method = ipcw.method,
                    higher.y = higher.y, abc = abc,
                    prop.cutoff = prop.cutoff, prop.multi = prop.multi,
                    ps.method = ps.method, minPS = minPS, maxPS = maxPS,
                    train.prop = train.prop, cv.n = cv.n, error.max = error.max, max.iter = max.iter,
                    initial.predictor.method = initial.predictor.method,
                    tree.depth = tree.depth, n.trees.rf = n.trees.rf, n.trees.boosting = n.trees.boosting, B = B, Kfold = Kfold,
                    error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune,
                    seed = seed, plot.gbmperf = plot.gbmperf, verbose = verbose)
  }
  if (response == "continuous"){
    cvout <- cvmean(cate.model = cate.model, init.model = init.model, ps.model = ps.model, data = data, score.method = score.method,
                    ipcw.model = ipcw.model,
                    higher.y = higher.y, abc = abc,
                    prop.cutoff = prop.cutoff, prop.multi = prop.multi,
                    ps.method = ps.method, minPS = minPS, maxPS = maxPS,
                    train.prop = train.prop, cv.n = cv.n, error.max = error.max, max.iter = max.iter,
                    initial.predictor.method = initial.predictor.method,
                    xvar.smooth.score = xvar.smooth.score, xvar.smooth.init = xvar.smooth.init,
                    tree.depth = tree.depth, n.trees.rf = n.trees.rf, n.trees.boosting = n.trees.boosting, B = B, Kfold = Kfold,
                    error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune,
                    seed = seed, plot.gbmperf = plot.gbmperf, verbose = verbose)
  }

  return(cvout)
}

#' Estimation of the conditional average treatment effect (CATE) score for count, survival and continuous data
#'
#' Provides singly robust and doubly robust estimation of CATE score for count, survival and continuous data with
#' up the following scoring methods among the following: Random forest (survival, continuous only), boosting,
#' poisson regression (count, survival only), two regressions, contrast regression, negative binomial regression (count only),
#' linear regression (continuous only), and generalized additive model (continuous only).
#'
#' @param response A string describing the type of outcome in the data. Allowed values include
#' "count" (see \code{\link{cvcount}()}), "survival" (see \code{\link{cvsurv}()}) and "continuous" (see \code{\link{cvmean}()}) .
#' @param cate.model A formula describing the outcome model to be fitted.
#' The outcome must appear on the left-hand side. For survival outcomes, a \code{Surv} object
#' must be used to describe the outcome.
#'@param init.model A formula describing the initial predictor model. The outcome must appear on the left-hand side.
#' It must be specified when \code{score.method = contrastReg} or \code{twoReg}.
#' @param ps.model A formula describing the propensity score (PS) model to be fitted. The treatment must
#' appear on the left-hand side. The treatment must be a numeric vector coded as 0/1.
#' If data are from a randomized controlled trial, specify \code{ps.model = ~1} as an intercept-only model.
#' @param data A data frame containing the variables in the outcome, propensity score, and inverse
#' probability of censoring models (if specified); a data frame with \code{n} rows (1 row per observation).
#' @param score.method A vector of one or multiple methods to estimate the CATE score.
#' Allowed values are: \code{'boosting'}, \code{'twoReg'}, \code{'contrastReg'}, \code{'poisson'} (count and survival outcomes only),
#' \code{'randomForest'} (survival, continuous outcomes only), \code{negBin} (count outcomes only), \code{'gam'} (continuous outcomes only),
#' \code{'gaussian'} (continuous outcomes only).
#' @param ipcw.model A formula describing the inverse probability of censoring weighting (IPCW)
#' model to be fitted. The left-hand side must be empty. Only applies for survival outcomes.
#' Default is \code{NULL}, which corresponds to specifying the IPCW with the same covariates
#' as the outcome model \code{cate.model}, plus the treatment.
#' @param followup.time A column name in \code{data} specifying the maximum follow-up time,
#' interpreted as the potential censoring time. Only applies for survival outcomes.
#' Default is \code{NULL}, which corresponds to unknown potential censoring time.
#' @param tau0 The truncation time for defining restricted mean time lost. Only applies for
#' survival outcomes. Default is \code{NULL}, which corresponds to setting the truncation time as the
#' maximum survival time in the data.
#' @param surv.min Lower truncation limit for the probability of being censored.
#' It must be a positive value and should be chosen close to 0. Only applies for survival
#' outcomes. Default is \code{0.025}.
#' @param ipcw.method A character value for the censoring model. Only applies for survival
#' outcomes. Allowed values are: \code{'breslow'} (Cox regression with Breslow estimator of t
#' he baseline survivor function), \code{'aft (exponential)'}, \code{'aft (weibull)'},
#' \code{'aft (lognormal)'} or \code{'aft (loglogistic)'} (accelerated failure time model
#' with different distributions for y variable). Default is \code{'breslow'}.
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
#' outcome predictions conditional on the covariates specified in \code{cate.model}. Only applies
#' when \code{score.method} includes \code{'twoReg'} or \code{'contrastReg'}.Allowed values include one of
#' \code{'randomForest'} (survival outcomes only), \code{'boosting'}, \code{'logistic'}
#' (survival outcomes only, fast), \code{'poisson'} (count outcomes only, fast), \code{'gaussian'} (continuous outcomes only) and
#' \code{'gam'} (count and continuous outcomes only). Default is \code{NULL}, which assigns \code{'boosting'}
#' for count outcomes and \code{'randomForest'} for survival outcomes.
#' @param xvar.smooth.init A vector of characters indicating the name of the variables used as
#' the smooth terms if \code{initial.predictor.method = 'gam'}. The variables must be selected
#' from the variables listed in \code{init.model}.
#' Default is \code{NULL}, which uses all variables in \code{init.model}.
#' @param xvar.smooth.score A vector of characters indicating the name of the variables used as
#' the smooth terms if \code{score.method = 'gam'}. The variables must be selected
#' from the variables listed in \code{cate.model}.
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
#' @param verbose An integer value indicating whether intermediate progress messages and histograms should
#' be printed. \code{1} indicates messages are printed and \code{0} otherwise. Default is \code{1}.
#'
#' @return For count response, see description of outputs in \code{\link{pmcount}()}.
#' For survival response, see description of outputs in \code{\link{pmsurv}()}.
#'
#' @details For count response, see details in \code{\link{pmcount}()}.
#' For survival response, see details in \code{\link{pmsurv}()}.
#'
#' @references Yadlowsky, S., Pellegrini, F., Lionetto, F., Braune, S., & Tian, L. (2020).
#' \emph{Estimation and validation of ratio-based conditional average treatment effects using
#' observational data. Journal of the American Statistical Association, 1-18.}
#' \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1772080}
#'
#' @seealso \code{\link{cv}()}
#'
#' @examples
#' \dontrun{
#' # Count outcome
#' output_pm <- pm(response = "count",
#'                 cate.model = y ~ age + female + previous_treatment +
#'                                      previous_cost + previous_number_relapses + offset(log(years)),
#'                 ps.model = trt ~ age + previous_treatment,
#'                 data = countExample,
#'                 higher.y = TRUE,
#'                 score.method = "poisson",
#'                 seed = 999)
#'
#' # Try:
#' coef(output_pm)
#'
#' # Survival outcome
#' tau0 <- with(survivalExample,
#'                  min(quantile(y[trt == "drug1"], 0.95), quantile(y[trt == "drug0"], 0.95)))
#' output_pm2 <- pm(response = "survival",
#'                  cate.model = survival::Surv(y, d) ~ age + female
#'                  + previous_cost + previous_number_relapses,
#'                  ps.model = trt ~ age + previous_treatment,
#'                  ipcw.model = ~ age + previous_cost + previous_treatment,
#'                  data = survivalExample,
#'                  higher.y = TRUE,
#'                  score.method = c("poisson", "boosting", "randomForest"),
#'                  followup.time = NULL,
#'                  tau0 = tau0,
#'                  surv.min = 0.025,
#'                  initial.predictor.method = "logistic",
#'                  seed = 999,
#'                  plot.gbmperf = FALSE)
#'
#' # Try:
#' coef(output_pm2)
#'
#'
#' # Continuous outcome
#' output_pm3 <- pm(response = "continuous",
#'                 cate.model = y ~ age +
#'                                  previous_treatment +
#'                                  previous_cost +
#'                                  previous_status_measure,
#'                 init.model = y ~ age +
#'                                  previous_treatment +
#'                                  previous_cost +
#'                                  previous_status_measure,
#'                 ps.model = trt ~ previous_treatment,
#'                 data = meanExample,
#'                 higher.y = FALSE,
#'                 score.method = c("gaussian", "randomForest", "twoReg", "contrastReg"),
#'                 initial.predictor.method = "gaussian",
#'                 seed = 999)
#'
#' # Try:
#' coef(output_pm3)
#'
#' }
#' @export

pm <- function(response, cate.model, ps.model, data, score.method,
               init.model = NULL, # This is used in continuous data, when contrast/two regression is used
               ipcw.model = NULL, followup.time = NULL, tau0 = NULL,
               surv.min = 0.025, ipcw.method = "breslow",
               higher.y = TRUE, prop.cutoff = seq(0.5, 1, length = 6),
               ps.method = "glm", minPS = 0.01, maxPS = 0.99,
               initial.predictor.method = NULL,  xvar.smooth.score = NULL, xvar.smooth.init = NULL,
               tree.depth = 2, n.trees.rf = 1000, n.trees.boosting = 200, B = 3, Kfold = 5,
               error.maxNR = 1e-3, max.iterNR = 150, tune = c(0.5, 2),
               seed = NULL, plot.gbmperf = TRUE, verbose = 1) {

  stopifnot("`response` must be either `count`, `survival`, or `continuous`." = any(response== c("count", "survival", "continuous")))

  if (response == "count") {
    pmout <- pmcount(cate.model = cate.model, ps.model = ps.model, data = data, score.method = score.method,
                     higher.y = higher.y,
                     prop.cutoff = prop.cutoff,
                     ps.method = ps.method, minPS = minPS, maxPS = maxPS,
                     initial.predictor.method = initial.predictor.method, xvar.smooth = xvar.smooth.score,
                     tree.depth = tree.depth, n.trees.boosting = n.trees.boosting, B = B, Kfold = Kfold,
                     error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune,
                     seed = seed, plot.gbmperf = plot.gbmperf)
  }
  if (response == "survival") {
    pmout <- pmsurv(cate.model = cate.model, ps.model = ps.model, data = data, score.method = score.method,
                    ipcw.model = ipcw.model, followup.time = followup.time, tau0 = tau0,
                    surv.min = surv.min, ipcw.method = ipcw.method,
                    higher.y = higher.y,
                    prop.cutoff = prop.cutoff,
                    ps.method = ps.method, minPS = minPS, maxPS = maxPS,
                    initial.predictor.method = initial.predictor.method,
                    tree.depth = tree.depth, n.trees.rf = n.trees.rf, n.trees.boosting = n.trees.boosting, B = B, Kfold = Kfold,
                    error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune,
                    seed = seed, plot.gbmperf = plot.gbmperf)
  }
  if (response == "continuous") {
    pmout <- pmmean(cate.model = cate.model, init.model = init.model, ps.model = ps.model, data = data, score.method = score.method,
                    ipcw.model = ipcw.model,
                    higher.y = higher.y, abc = abc,
                    prop.cutoff = prop.cutoff,
                    ps.method = ps.method, minPS = minPS, maxPS = maxPS,
                    initial.predictor.method = initial.predictor.method,
                    xvar.smooth.score = xvar.smooth.score, xvar.smooth.init = xvar.smooth.init,
                    tree.depth = tree.depth, n.trees.rf = n.trees.rf, n.trees.boosting = n.trees.boosting, B = B, Kfold = Kfold,
                    error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune,
                    seed = seed, plot.gbmperf = plot.gbmperf, verbose = verbose)
  }
  return(pmout)
}

#' Doubly robust estimator of and inference for the average treatment effect for count, survival and continuous data
#'
#' Doubly robust estimator of the average treatment effect between two treatments, which is the rate ratio
#' for count outcomes, the restricted mean time lost ratio for survival outcomes and the mean difference for continuous outcome. Bootstrap is used for
#' inference.
#'
#' @param response A string describing the type of outcome in the data. Allowed values include
#' "count" (see \code{\link{cvcount}()}), "survival" (see \code{\link{cvsurv}()}) and "continuous" (see \code{\link{cvmean}()}).
#' @param cate.model A formula describing the outcome model to be fitted.
#' The outcome must appear on the left-hand side. For survival outcomes, a \code{Surv} object
#' must be used to describe the outcome.
#' @param ps.model A formula describing the propensity score (PS) model to be fitted. The treatment must
#' appear on the left-hand side. The treatment must be a numeric vector coded as 0/1.
#' If data are from a randomized controlled trial, specify \code{ps.model = ~1} as an intercept-only model.
#' @param data A data frame containing the variables in the outcome, propensity score, and inverse
#' probability of censoring models (if specified); a data frame with \code{n} rows (1 row per observation).
#' @param ipcw.model A formula describing the inverse probability of censoring weighting (IPCW)
#' model to be fitted. The left-hand side must be empty. Only applies for survival outcomes.
#' Default is \code{NULL}, which corresponds to specifying the IPCW with the same covariates
#' as the outcome model \code{cate.model}, plus the treatment.
#' @param followup.time A column name in \code{data} specifying the maximum follow-up time,
#' interpreted as the potential censoring time. Only applies for survival outcomes.
#' Default is \code{NULL}, which corresponds to unknown potential censoring time.
#' @param tau0 The truncation time for defining restricted mean time lost. Only applies for
#' survival outcomes. Default is \code{NULL}, which corresponds to setting the truncation time as the
#' maximum survival time in the data.
#' @param surv.min Lower truncation limit for the probability of being censored.
#' It must be a positive value and should be chosen close to 0. Only applies for survival
#' outcomes. Default is \code{0.025}.
#' @param ipcw.method A character value for the censoring model. Only applies for survival
#' outcomes. Allowed values are: \code{'breslow'} (Cox regression with Breslow estimator of t
#' he baseline survivor function), \code{'aft (exponential)'}, \code{'aft (weibull)'},
#' \code{'aft (lognormal)'} or \code{'aft (loglogistic)'} (accelerated failure time model
#' with different distributions for y variable). Default is \code{'breslow'}.
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
#' between \code{x} and \code{trt}. Applies only to count outcomes. If \code{TRUE}, interactions will
#' be assumed only if at least 10 patients received each treatment option. Default is \code{TRUE}.
#' @param n.boot A numeric value indicating the number of bootstrap samples used. Default is \code{500}.
#' @param seed An optional integer specifying an initial randomization seed for reproducibility.
#' Default is \code{NULL}, corresponding to no seed.
#' @param verbose An integer value indicating whether intermediate progress messages and histograms should
#' be printed. \code{1} indicates messages are printed and \code{0} otherwise. Default is \code{1}.
#' @param plot.boot A logical value indicating whether histograms of the bootstrapped log(rate ratio)
#' (for count outcomes) log(restricted mean time lost ratio) (for survival outcomes) should be produced at
#' every \code{n.boot/10}-th iteration and whether the final histogram should be outputted. This argument is
#' only taken into account if \code{verbose = 1}. Default is \code{FALSE}.
#'
#' @return For count response, see description of outputs in \code{\link{drcount.inference}()}.
#' For survival response, see description of outputs in \code{\link{drsurv.inference}()}.
#'
#' @details For count response, see details in \code{\link{drcount.inference}()}.
#' For survival response, see details in \code{\link{drsurv.inference}()}.
#'
#' @examples
#' \dontrun{
#' # Count outcome
#' output <- dr.inference(response = "count",
#'                        cate.model = y ~ age + female + previous_treatment +
#'                        previous_cost + previous_number_relapses + offset(log(years)),
#'                        ps.model = trt ~ age + previous_treatment,
#'                        data = countExample,
#'                        plot.boot = TRUE,
#'                        seed = 999)
#' print(output)
#' output$plot
#'
#' # Survival outcome
#' tau0 <- with(survivalExample,
#'                  min(quantile(y[trt == "drug1"], 0.95), quantile(y[trt == "drug0"], 0.95)))
#' output2 <- dr.inference(response = "survival",
#'                         cate.model = survival::Surv(y, d) ~ age + female +
#'                         previous_cost + previous_number_relapses,
#'                         ps.model = trt ~ age + previous_treatment,
#'                         data = survivalExample,
#'                         tau0 = tau0,
#'                         plot.boot = TRUE,
#'                         seed = 999)
#' print(output2)
#' output2$plot
#'
#'# Continuous outcome
#' output3 <- dr.inference(response = "continuous",
#'                        cate.model = y ~ age +
#'                                         previous_treatment +
#'                                         previous_cost +
#'                                         previous_status_measure,
#'                        previous_cost + previous_number_relapses + offset(log(years)),
#'                        ps.model = trt ~ previous_treatment,
#'                        data = meanExample,
#'                        plot.boot = TRUE,
#'                        seed = 999)
#' print(output3)
#' output$plot
#'}
#' @export
#'
#' @importFrom ggplot2 ggplot geom_histogram geom_vline
#' @importFrom dplyr mutate
#' @importFrom tidyr gather

dr.inference <- function(response, cate.model, ps.model, data,
                         ipcw.model = NULL, followup.time = NULL, tau0 = NULL, surv.min = 0.025, ipcw.method = "breslow",
                         ps.method = "glm", minPS = 0.01, maxPS = 0.99, interactions = TRUE,
                         n.boot = 500, seed = NULL, verbose = 1, plot.boot = FALSE) {

  stopifnot("`response` must be either `count`, `survival`, or `continuous`." = any(response== c("count", "survival", "continuous")))

  if (response == "count"){
    drout <- drcount.inference(cate.model = cate.model, ps.model = ps.model, data = data,
                               ps.method = ps.method, minPS = minPS, maxPS = maxPS, interactions = interactions,
                               n.boot = n.boot, verbose = verbose, plot.boot = plot.boot, seed = seed)
  }
  if (response == "survival"){
    drout <- drsurv.inference(cate.model = cate.model, ps.model = ps.model, data = data,
                              ipcw.model = ipcw.model, followup.time = followup.time, tau0 = tau0,
                              surv.min = surv.min, ipcw.method = ipcw.method,
                              ps.method = ps.method, minPS = minPS, maxPS = maxPS,
                              n.boot = n.boot, verbose = verbose, plot.boot = plot.boot, seed = seed)
  }
  if (response == "continuous"){
    drout <- drmean.inference(cate.model = cate.model, ps.model = ps.model, data = data,
                              ps.method = ps.method, minPS = minPS, maxPS = maxPS, interactions = interactions,
                              n.boot = n.boot, verbose = verbose, plot.boot = plot.boot, seed = seed)
  }
  return(drout)

}


#' Compute the area between curves from the \code{"precmed"} object
#'
#' Compute the area between curves (ABC) for each scoring method in the \code{"precmed"} object.
#' This should be run only after results of \code{\link{cv}()} have been obtained.
#'
#' @param x An object of class \code{"precmed"}.
#'
#' @return Returns a matrix of numeric values with number of columns equal to the number cross-validation
#' iteration and number of rows equal to the number of scoring methods in \code{x}.
#'
#' @details The ABC is the area between a validation curve and the overall ATE in the validation set.
#' It is calculated for each scoring method separately. Higher ABC values are preferable as they
#' indicate that more treatment effect heterogeneity is captured by the scoring method.
#' Negative values of ABC are possible if segments of the validation curve cross the overall ATE line.
#' The ABC is calculated with the \code{\link[MESS]{auc}()} in the \code{MESS} package with a natural
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
#' \url{https://www.jstor.org/stable/24246461?seq=1}
#'
#' @seealso \code{\link{cv}()} function and \code{\link{plot}()}, \code{\link{boxplot}()} methods for
#' \code{"precmed"} objects.
#'
#' @examples
#' \dontrun{
#' # Count outcome
#' cv_count <- cv(response = "count",
#'                cate.model = y ~ age + female + previous_treatment +
#'                                 previous_cost + previous_number_relapses + offset(log(years)),
#'                ps.model = trt ~ age + previous_treatment,
#'                data = countExample,
#'                higher.y = FALSE,
#'                score.method = "poisson",
#'                cv.n = 5,
#'                plot.gbmperf = FALSE)
#'
#' abc(x = cv_count) # ABC of the validation curves for each method and each CV iteration
#'
#' # Survival outcome
#' cv_surv <- cv(response = "survival",
#'               cate.model = survival::Surv(y, d) ~ age + female
#'                                                   + previous_cost + previous_number_relapses,
#'               ps.model = trt ~ age + previous_treatment,
#'               data = survivalExample,
#'               score.method = c("poisson", "randomForest"),
#'               higher.y = FALSE,
#'               cv.n = 5,
#'               plot.gbmperf = FALSE)
#'
#' abc(x = cv_surv) # ABC of the validation curves for each method and each CV iteration
#'
#' # Continuous outcome
#' cv.mean <- cvmean(cate.model = y ~ age +
#'                                    previous_treatment +
#'                                    previous_cost +
#'                                    previous_status_measure,
#'                   init.model = y ~ age +
#'                                    previous_treatment +
#'                                    previous_cost +
#'                                    previous_status_measure,
#'                   ps.model = trt ~ previous_treatment,
#'                   data = meanExample,
#'                   higher.y = FALSE,
#'                   score.method = "gaussian",
#'                   cv.n = 5,
#'                   plot.gbmperf = FALSE)
#'
#'
#' abc(x = cv.mean) # ABC of the validation curves for each method and each CV iteration
#' }
#'
#' @export
#'
#' @importFrom stringr str_extract

abc <- function(x) {

  # Check the value of abc (must have area between curves values)
  if (x$abc == FALSE) {
    stop("Area between curves (ABC) must have been calculated in x by setting abc = TRUE in cv().")
  }

  # Detect the number of score methods
  map <- c("boosting" = "Boosting",
           "poisson" = "Naive Poisson",
           "twoReg" = "Two Regressions",
           "contrastReg" = "Contrast Regression",
           "negBin" = "Negative Binomial",
           "randomForest" = "Random Forest")
  score.method <- str_extract(names(x), "(?<=^ate\\.).*$")
  score.method <- score.method[!is.na(score.method)]
  m <- length(score.method)

  # Value of ABC
  value.abc <- matrix(NA, nrow = m, ncol = x$cv.n)
  for (i in 1:m){
    value.abc[i,] <- x[[paste0("ate.", score.method[i])]]$abc.valid
  }
  rownames(value.abc) <- score.method
  colnames(value.abc) <- paste0("cv", 1:(x$cv.n))

  return(value.abc)
}
