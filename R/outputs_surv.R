# ------------------------------------------------------------------
# Project: Precision Medicine MS - Comprehensive R package
#
# Purpose: Output functions for Survival outcomes
#
# Platform: Windows
# R Version: 4.1.0
#
#   Modifications:
#
#   Date			By			Description
# --------		--------	-----------------------------
#  14JUN2021  dl      Start the script
#  09JUL2021  dl      Add cvsurv()
#  12JUL2021  dl      Add pmsurv()
#  13JUL2021  dl      Change cvsurv argument: 1) y, d, trt, x to cate.model and ps.model. y, d, trt, x.ps, x.cate are extracted
#                     within the function using data.preproc. 2) prop.onlyhigh and prop.bi to prop.cutoff so that prop.onlyhigh
#                     and prop.bi are also extracted from data.preproc;
#                     Change pmsurv argument: 1) y, d, trt, x to cate.model and ps.model. y, d, trt, x.ps, x.cate are extracted
#                     within the function using data.preproc. 2) prop to prop.cutoff so that data.preproc will extract prop.onlyhigh;
#  13JUL2021  gs      Review cvsurv() with minor edits and pmsurv()
#  14JUL2021  dl      Add an argument higher.y to cvsurv() and pmsurv();
#                     Add an argument abc to cvsurv()
#  15JUL2021  dl      Add plotsurv.PrecMed() and boxplotsurv.PrecMed();
#                     Add drsurv.inference() as a separate function
#  16JUL2021  dl      Add ipcw.model argument to cvsurv and pmsurv;
#                     Change argument yf to followup.time;
#                     Add arg.check.surv() inside cvsurv() and pmsurv() to check argument
#  20JUL2021  dl      Removed tau0 from balancesurv.split since unused
#  28JUL2021  dl      Add warning to cvsurv and pmsurv if tau0 is less than the median of the time (either observed or censored);
#                     Minor updates in drsurv.inference
#  04AUG2021  dl      Added error-handling code for subgroup ATE estimation
#  05AUG2021  dl      Added plot.hr, valid.only, and legend.position option to plotsurv.PrecMed
#  13AUG2021  dl      Added cvsurv2 function (see improvement document);
#                     Added plotsurv.PrecMed2 function to handle cvsurv2 output;
#                     Modified boxplotsurv.PrecMed and added boxplot.bi.surv.PrecMed
#  26AUG2021  gs      Add initial.predictor.method in arg.check
#  09SEP2021  gs      Add best.iter in cvsurv(), cvsurv2() and pmsurv()
#                     Minor edits in cvsurv2() and plot2() (syntax, warning messages, line break)
#  29SEP2021  pj      Remove argument seed.cf
#  05OCT2021  gs      Add ipcw.method argument to survprd -> inherited in cvsurv, pmsurv, drsurv.inference
#                     Add seed argument to pmsurv (like in pmcount)
#  06OCT2021  pj      Remove set.seed functions inside and leave only one set.seed in cvsurv/pmsurv/drsurv.inference once
#                     Change all seeds to NULL to be consistent
#  13OCT2021  pj      Reward the warning messages (iteration [0-9]) and fix typos
#  15OCT2021  gs      Lower case ipcw.method argument names
#  20OCT2021  pj      Add error/warning results to outputed list in pmsurv
#                     Calculate the ATE when prop = 1 only once so estsurv calculates prop.no1 in pmsurv to save time
#  22OCT2021  gs      Add initial.predictor logistic regression
#  03NOV2021  pj      Replace est.train with est.high in pm errors/warnings
#  05NOV2021  gs      Fix ATE calculation for higher.y = TRUE vs FALSE in plot()
#  12NOV2021  gs      Fix higher.y in abc calculation and in cv() documentation
#                     Update documentation for cvsurv()
#  19NOV2021  gs      Add formulas and original treatment labels in output of cvcount() for more meaningful y-axis in plot
#                     Add more meaningful y-axis in plot() and boxplot() if default is used
#  01DEC2021  pj      Update documentation of cvsurve() and pmsurv()
#  03DEC2021  gs      Update prop.multi documentation
#                     Fix y-axis title in plot when plot.hr = TRUE
#                     Clarify subgroup label in boxplot
#  09DEC2021  pj      Update plot documentation and minor changes
#  13DEC2021  gs      Fix axes in boxplot, remove legend
#                     Always extract overall ATE (RMTL or HR) by cv iterations in cvsurv() - to be used as overall ATE in boxplot
#  15DEC2021  pj      Minor changes
#  17DEC2021  gs      Always plot overall ATE in lineplot
#  07JAN2022  pj      Revise plotsurv.PrecMed() so it accepts both count and survival outcomes
#  19JAN2022  pj      Revise documentation for plotsurv.PrecMed() for both count and survival outcomes
#                     Revise boxplotsurv.PrecMed() so it accepts both count and survival outcomes
#  20JAN2022  gs      Revise plotsurv.PrecMed() to ensure that the ATE line is accurate
#  27JAN2022  pj      Revise ylab, xlim, and abc results as well as documentation in plot and boxplot
#  01FEB2022  gs      Output cv.n from cvsurv; update error message in plotsurv.PrecMed() when cv.i > cv.n
#                     Location of plot.hr + count outcome argument check fixed in plotsurv.PrecMed()
#  04FEB2022  pj      Add missing argument in documentation and fix \link and example in documentation
#  07FEB2022  pj      Revise documentation examples to pass testing (fewer score method, small cv.n, rename plots)
#  18FEB2022  gs      Edits in boxplot after testing
#  22FEB2022  gs      Minor edits in boxplot() documentation
#  23FEB2022  gs      Build skeleton of cv() based on cvsurv arguments
#  03MAR2022  pj      Merge plot from surv and count to one plot.PrecMed()
#  04MAR2022  gs      Skeleton pm() and dr.inference()
#                     Update arguments and documentation cv(), pm(), dr.inference()
#                     Change all references to gbm to boosting (except plot.gbmperf)
#  09MAR2022  pj      Merge boxplot from surv and count to one boxplot.PrecMed()
#                     Add response argument to cv/pm and fill in the cv()/pm() function
#  11MAR2022  gs      Update return, detail, example in wrapper functions
#                     Update cvsurv, pmsurv, drsurv.inf arguments & documentation to match wrapper
#  16MAR2022  pj      Complete dr.inference wrapper
#  24MAR2022  pj      Create wrapper for arg.checks() to distinguish common and specific args based on outcome type
#  31MAR2022  gs      Reverse x-axis on boxplot
#  13APR2022  pj      Minor edits in doc, add more examples to doc
#  14APR2022  gs      Remove follow.up line of code at beginning of cvsurv
#  28APR2022  gs      Change "ITR score" to "Estimated CATE" in the plots
#  03MAY2022  pj      Revise cv and pm examples in doc
#  06MAY2022  gs      Check argument & data.preproc with initial.predictor.method with default NULL
#                     Check argument & data.preproc with tau0 default NULL
#  10MAY2022  gs      Add tau0 argument to data.preproc
#  18MAY2022  gs      Remove follow.up line of code at beginning of drsurv.inference
#  25MAY2022  gs      Minor changes in dr.inference (output warning, change output)
#  25MAY2022  pj      Revise verbose argument in cv() and dr.inference() to integer values
#  30MAY2022  gs      Add survival example in dr.inference
#  01JUn2022  pj      Fix dr.inference() verbose tp control only text outputs
#  11JUL2022  gs      Minor edits to cvsurv() output
#  13JUL2022  gs      Check pm() arguments
#  05JUL2022  sk      Included continuous case in plot
#  22JUL2022  gs      Update drinf function description
#  02AUG2022  sk      Update cv() function, need to consider whether to include init.model (it is included for continuous case for now)
#  15AUG2022  gs      Transfer abc() from output.R, documentation update
# ------------------------------------------------------------------

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
#' \code{"PrecMed"} objects.
#'
#' @examples
#' # Count outcome
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
#'
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
#' output_pm$coefficients
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
#'                  score.method = c("poisson", "boosting", "randomForest", "twoReg", "contrastReg"),
#'                  followup.time = NULL,
#'                  tau0 = tau0,
#'                  surv.min = 0.025,
#'                  initial.predictor.method = "logistic",
#'                  seed = 999,
#'                  plot.gbmperf = FALSE)
#'
#' # Try:
#' output_pm2$coefficients
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
#' output_pm3$coefficients
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
               seed = NULL, plot.gbmperf = TRUE) {

  if (response == "count"){
    pmout <- pmcount(cate.model = cate.model, ps.model = ps.model, data = data, score.method = score.method,
                     higher.y = higher.y,
                     prop.cutoff = prop.cutoff,
                     ps.method = ps.method, minPS = minPS, maxPS = maxPS,
                     initial.predictor.method = initial.predictor.method, xvar.smooth = xvar.smooth.score,
                     tree.depth = tree.depth, n.trees.boosting = n.trees.boosting, B = B, Kfold = Kfold,
                     error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune,
                     seed = seed, plot.gbmperf = plot.gbmperf)
  }
  if (response == "survival"){
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
  if (response == "continuous"){
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
#'
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
#'
#' @export
#'
#' @importFrom ggplot2 ggplot geom_histogram geom_vline
#' @importFrom dplyr mutate
#' @importFrom tidyr gather

dr.inference <- function(response, cate.model, ps.model, data,
                         ipcw.model = NULL, followup.time = NULL, tau0 = NULL, surv.min = 0.025, ipcw.method = "breslow",
                         ps.method = "glm", minPS = 0.01, maxPS = 0.99, interactions = TRUE,
                         n.boot = 500, seed = NULL, verbose = 1, plot.boot = FALSE) {

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


#' Compute the area between curves from the \code{"PrecMed"} object
#'
#' Compute the area between curves (ABC) for each scoring method in the \code{"PrecMed"} object.
#' This should be run only after results of \code{\link{cv}()} have been obtained.
#'
#' @param x An object of class \code{"PrecMed"}.
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
#' \code{"PrecMed"} objects.
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
#' abc(x = cv.mean)
#' }
#'
#' @export
#'
#' @importFrom stringr str_extract

abc <- function(x) {

  # Check the value of abc (must have sbc values)
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
#' \code{2} means progress bar, run time, and all errors and warnings. Default is \code{2}.
#'
#' @return Returns a list containing the following components saved as a \code{"PrecMed"} object:
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
#'  \item{\code{overall.hr.valid}: }{A vector of numerical values of length \code{cv.n}.
#'  The ith element contains the ATE (HR) in the validation set of the ith cross-validation
#'  iteration.}
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
#' @seealso \code{\link{pmsurv}()} function and \code{\link{boxplot}()}, \code{\link{abc}} methods for
#' \code{"PrecMed"} objects.
#'
#' @examples
#'
#' tau0 <- with(survivalExample,
#'              min(quantile(y[trt == "drug1"], 0.95), quantile(y[trt == "drug0"], 0.95)))
#'
#' cv <- cvsurv(cate.model = survival::Surv(y, d) ~ age + female
#'                                                      + previous_cost + previous_number_relapses,
#'              ps.model = trt ~ age + previous_treatment,
#'              ipcw.model = ~ age + previous_cost + previous_treatment,
#'              data = survivalExample,
#'              followup.time = NULL,
#'              tau0 = tau0,
#'              surv.min = 0.025,
#'              higher.y = TRUE,
#'              score.method = c("poisson"),
#'              cv.n = 5,
#'              initial.predictor.method = "randomForest",
#'              plot.gbmperf = FALSE,
#'              seed = 999)
#'
#' # Try:
#' plot(x = cv, ylab = "RMTL ratio of drug1 vs drug0 in each subgroup")
#' boxplot(x = cv, ylab = "RMTL ratio of drug1 vs drug0 in each subgroup")
#' abc(x = cv)
#'
#' @export
#'
#' @importFrom stats as.formula coef glm median model.frame model.matrix model.response na.omit optim pchisq predict qnorm quantile sd var
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom survival Surv coxph coxph.detail survreg
#' @importFrom randomForestSRC rfsrc predict.rfsrc
#' @importFrom gbm gbm gbm.perf predict.gbm
#' @importFrom gam gam
#' @importFrom stringr str_replace
#' @importFrom MESS auc
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom dplyr %>%

cvsurv <- function(cate.model, ps.model, data, score.method,
                   ipcw.model = NULL, tau0 = NULL, followup.time = NULL,
                   surv.min = 0.025, ipcw.method = "breslow",
                   higher.y = TRUE, abc = TRUE,
                   prop.cutoff = seq(0.5, 1, length = 6), prop.multi = c(0, 1/3, 2/3, 1),
                   ps.method = "glm", minPS = 0.01, maxPS = 0.99,
                   train.prop = 3/4, cv.n = 10, error.max = 0.1, max.iter = 5000,
                   initial.predictor.method = "randomForest",
                   tree.depth = 2, n.trees.rf = 1000, n.trees.boosting = 200, B = 3, Kfold = 5,
                   error.maxNR = 1e-3, max.iterNR = 150, tune = c(0.5, 2),
                   seed = NULL, plot.gbmperf = TRUE, verbose = 2) {

  # Set seed for reproducibility
  set.seed(seed)

  if (verbose >= 1) t.start <- Sys.time()

  #### CHECK ARGUMENTS ####
  arg.checks(
    fun = "cv", response = "survival", data = data, followup.time = followup.time, tau0 = tau0, surv.min = surv.min,
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


  #### PRE-PROCESSING ####
  out <- data.preproc.surv(fun = "cv", cate.model = cate.model, ps.model = ps.model, ipcw.model = ipcw.model, tau0 = tau0,
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
        if (verbose == 2) cat(paste0('    Error(s) occurred when estimating the ATEs in the nested subgroup using "', paste0(errors.train, collapse = '", "'), '";\n    return NAs in the corresponding subgroup.'),'\n')
        warning(paste0('Error(s) occurred when estimating the ATEs in the nested subgroup in the training set using "', paste0(errors.train, collapse = '", "'), '" in cross-validation iteration ', cv.i, ". NAs are returned for RMTL ratio and HR in the corresponding subgroup; see 'errors/warnings'."))
      }

      if (length(warnings.train) != 0) {
        if (verbose == 2) cat(paste0('    Warning(s) occurred when estimating the ATEs in the nested subgroup using "', paste0(warnings.train, collapse = '", "'), '"'),'\n')
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
        if (verbose == 2) cat(paste0('    Error(s) occurred when estimating the ATEs in the nested subgroup using "', paste0(errors.valid, collapse = '", "'), '";\n    return NAs in the corresponding subgroup.'),'\n')
        warning(paste0('Error(s) occurred when estimating the ATEs in the nested subgroup in the validation set using "', paste0(errors.valid, collapse = '", "'), '" in cross-validation iteration ', cv.i, ". NAs are returned for RMTL ratio and HR in the corresponding subgroup; see 'errors/warnings'."))
      }

      if (length(warnings.valid) != 0) {
        if (verbose == 2) cat(paste0('    Warning(s) occurred when estimating the ATEs in the nested subgroup using "', paste0(warnings.valid, collapse = '", "'), '"'),'\n')
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
    t.diff <- difftime(t.end, t.start)
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

  class(result) <- "PrecMed"

  return(result)
}


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
#'
#' @return Returns a list containing the following components:
#' \itemize{
#'  \item{\code{ate.randomForest}: }{A vector of numerical values of length \code{prop.cutoff}
#'  containing the estimated ATE by the RMTL ratio in nested subgroups (defined by \code{prop.cutoff})
#'  constructed based on the estimated CATE scores with random forest method.
#'  Only provided if \code{score.method} includes \code{'randomForest'}.}
#'  \item{\code{ate.boosting}: }{Same as \code{ate.randomForest}, but with the nested subgroups based
#'  the estimated CATE scores with boosting. Only provided if \code{score.method}
#'  includes \code{'boosting'}.}
#'  \item{\code{ate.poisson}: }{Same as \code{ate.randomForest}, but with the nested subgroups based
#'  the estimated CATE scores with poisson regression.
#'  Only provided if \code{score.method} includes \code{'poisson'}.}
#'  \item{\code{ate.twoReg}: }{Same as \code{ate.randomForest}, but with the nested subgroups based
#'  the estimated CATE scores with two regressions.
#'  Only provided if \code{score.method} includes \code{'twoReg'}.}
#'  \item{\code{ate.contrastReg}: }{Same as \code{ate.randomForest}, but with the nested subgroups based
#'  the estimated CATE scores with contrast regression.
#'  Only provided if \code{score.method} includes \code{'contrastReg'}.}
#'  \item{\code{hr.randomForest}: }{A vector of numerical values of length \code{prop.cutoff}
#'  containing the adjusted hazard ratio in nested subgroups (defined by \code{prop.cutoff})
#'  constructed based on the estimated CATE scores with random forest method.
#'  Only provided if \code{score.method} includes \code{'randomForest'}.}
#'  \item{\code{hr.boosting}: }{Same as \code{hr.randomForest}, but with the nested subgroups based
#'  the estimated CATE scores with boosting. Only provided if \code{score.method}
#'  includes \code{'boosting'}.}
#'  \item{\code{hr.poisson}: }{Same as \code{hr.randomForest}, but with the nested subgroups based
#'  the estimated CATE scores with poisson regression.
#'  Only provided if \code{score.method} includes \code{'poisson'}.}
#'  \item{\code{hr.twoReg}: }{Same as \code{hr.randomForest}, but with the nested subgroups based
#'  the estimated CATE scores with two regressions.
#'  Only provided if \code{score.method} includes \code{'twoReg'}.}
#'  \item{\code{hr.contrastReg}: }{Same as \code{hr.randomForest}, but with the nested subgroups based
#'  the estimated CATE scores with contrast regression.
#'  Only provided if \code{score.method} includes \code{'contrastReg'}.}
#'  \item{\code{score.randomForest}: }{A vector of numerical values of length n
#'  (number of observations in \code{data}) containing the estimated log-CATE scores
#'  according to random forest. Only provided if \code{score.method}
#'  includes \code{'randomForest'}.}
#'  \item{\code{score.boosting}: }{Same as \code{score.randomForest}, but with estimated log-CATE score
#'  according to boosting. Only provided if \code{score.method} includes
#'  \code{'boosting'}.}
#'  \item{\code{score.poisson}: }{Same as \code{score.randomForest}, but with estimated log-CATE score
#'  according to the Poisson regression. Only provided if \code{score.method}
#'  includes \code{'poisson'}.}
#'  \item{\code{score.twoReg}: }{Same as \code{score.randomForest}, but with estimated log-CATE score
#'  according to two regressions. Only provided if \code{score.method} includes
#'  \code{'twoReg'}.}
#'  \item{\code{score.contrastReg}: }{Same as \code{score.randomForest}, but with estimated log-CATE score
#'  according to contrast regression. Only provided if \code{score.method} includes
#'  \code{'contrastReg'}.}
#'  \item{\code{fit}: }{Additional details on model fitting if \code{score.method}
#'  includes 'randomForest', 'boosting' or 'contrastReg':}
#'  \itemize{
#'    \item{\code{result.randomForest}: }{Details on the random forest model fitted to observations
#'    with treatment = 0 \code{($fit0.rf)} and to observations with treatment = 1 \code{($fit1.rf)}.
#'    Only provided if \code{score.method} includes \code{'randomForest'}.}
#'    \item{\code{result.boosting}: }{Details on the boosting model fitted to observations
#'    with treatment = 0, \code{($fit0.boosting)} and \code{($fit0.gam)} and to observations with treatment = 1,
#'    \code{($fit1.boosting)} and \code{($fit1.gam)}.
#'    Only provided if \code{score.method} includes \code{'boosting'}.}
#'    \item{\code{result.contrastReg$converge.contrastReg}: }{Whether the contrast regression algorithm converged
#'    or not. Only provided if \code{score.method} includes \code{'contrastReg'}.}
#'  }
#'  \item{\code{coefficients}: }{A data frame with the coefficients of the estimated log-CATE
#'  score by \code{score.method}. The data frame has number of rows equal to the number of
#'  covariates in \code{cate.model} and number of columns equal to \code{length(score.method)}.
#'  If \code{score.method} includes \code{'contrastReg'}, the data frame has an additional
#'  column containing the standard errors of the coefficients estimated with contrast regression.
#'  \code{'randomForest'} and \code{'boosting'} do not have coefficient results because
#'  tree-based methods typically do not express the log-CATE as a linear combination of coefficients
#'  and covariates.}
#'  \item{\code{errors/warnings}: }{A nested list of errors and warnings that were wrapped during the
#'  calculation of ATE. Errors and warnings are organized by \code{score.method}.}
#' }
#'
#' @details The CATE score represents an individual-level treatment effect for survival data,
#' estimated with random forest, boosting, Poisson regression, and the doubly
#' robust estimator (two regressions, Yadlowsky, 2020) applied separately by treatment group
#' or with the other doubly robust estimators (contrast regression, Yadlowsky, 2020) applied
#' to the entire data set.
#'
#' \code{\link{pmsurv}()} provides the coefficients of the CATE score for each scoring method requested
#' through \code{score.method}. Currently, contrast regression is the only method which allows
#' for inference of the CATE coefficients by providing standard errors of the coefficients.
#' The coefficients can be used to learn the effect size of each variable and predict the
#' CATE score for a new observation.
#'
#' \code{\link{pmsurv}()} also provides the predicted CATE score of each observation in the data set,
#' for each scoring method. The predictions allow ranking the observations from potentially
#' high responders to the treatment to potentially low or standard responders.
#'
#' The estimated ATE among nested subgroups of high responders are also provided by scoring method.
#' Note that the ATEs in \code{\link{pmsurv}()} are derived based on the CATE score which is estimated
#' using the same data sample. Therefore, overfitting may be an issue. \code{\link{cvsurv}()} is more
#' suitable to inspect the estimated ATEs across scoring methods as it implements internal cross
#' validation to reduce optimism.
#'
#' @references Yadlowsky, S., Pellegrini, F., Lionetto, F., Braune, S., & Tian, L. (2020).
#' \emph{Estimation and validation of ratio-based conditional average treatment effects using
#' observational data. Journal of the American Statistical Association, 1-18.}
#' \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1772080}
#'
#' @seealso \code{\link{cvsurv}()}
#'
#' @examples
#'
#' tau0 <- with(survivalExample, min(quantile(y[trt == 1], 0.95), quantile(y[trt == 0], 0.95)))
#'
#' pm <- pmsurv(cate.model = survival::Surv(y, d) ~ age +
#'                                                  female +
#'                                                  previous_cost +
#'                                                  previous_number_relapses,
#'                           ps.model = trt ~ age + previous_treatment,
#'                           ipcw.model = ~ age + previous_cost + previous_treatment,
#'                           data = survivalExample,
#'                           tau0 = tau0,
#'                           score.method = "randomForest",
#'                           seed = 999)
#'
#'
#' @export

pmsurv <- function(cate.model, ps.model, score.method, data,
                   ipcw.model = NULL, followup.time = NULL, tau0 = NULL,
                   surv.min = 0.025, ipcw.method = "breslow",
                   higher.y = TRUE,
                   prop.cutoff = seq(0.5, 1, length = 6),
                   ps.method = "glm", minPS = 0.01, maxPS = 0.99,
                   initial.predictor.method = "randomForest",
                   tree.depth = 2, n.trees.rf = 1000, n.trees.boosting = 200, B = 3, Kfold = 5, plot.gbmperf = TRUE,
                   error.maxNR = 1e-3, max.iterNR = 100, tune = c(0.5, 2), seed = NULL) {

  # Set seed for reproducibility
  set.seed(seed)

  t.start <- Sys.time()
  # followup.time <- as.list(match.call()[-1])$followup.time

  #### CHECK ARGUMENTS ####
  arg.checks(
    fun = "pm", response = "survival", data = data, followup.time = followup.time, tau0 = tau0, surv.min = surv.min,
    higher.y = higher.y, score.method = score.method, abc = abc,
    prop.cutoff = prop.cutoff,
    ps.method = ps.method, minPS = minPS, maxPS = maxPS, ipcw.method = ipcw.method,
    train.prop = 0.5, cv.n = 1,
    error.max = 1, max.iter = 2,
    initial.predictor.method = initial.predictor.method,
    tree.depth = tree.depth, n.trees.rf = n.trees.rf, n.trees.boosting = n.trees.boosting,
    B = B, Kfold = Kfold, plot.gbmperf = plot.gbmperf,
    error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune
  )

  #### PRE-PROCESSING ####
  out <- data.preproc.surv(fun = "pm", cate.model = cate.model, ps.model = ps.model, ipcw.model = ipcw.model, tau0 = tau0,
                           data = data, prop.cutoff = prop.cutoff,
                           ps.method = ps.method, response = "survival")
  y <- out$y
  d <- out$d
  trt <- out$trt
  x.ps <- out$x.ps
  x.cate <- out$x.cate
  x.ipcw <- out$x.ipcw
  prop <- out$prop
  prop.no1 <- out$prop.no1

  if (is.null(followup.time)) {
    yf <- NULL
  } else {
    yf <- data[[followup.time]]
  }

  # Check if tau0 is large enough, i.e., if tau0 is larger than the 50% quantile of the observed or censoring time.
  if (tau0 < median(y)) warning("It is recommended to increase tau0 close to the largest observed or censored time.")

  #### FUNCTION STARTS HERE ####
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

  ###### Fit the interaction model --------------------------------------------------------------
  fit <- intxsurv(y = y, d = d, trt = trt, x.cate = x.cate, x.ps = x.ps, x.ipcw = x.ipcw, yf = yf, tau0 = tau0, surv.min = surv.min,
                  score.method = score.method,
                  ps.method = ps.method, minPS = minPS, maxPS = maxPS, ipcw.method = ipcw.method,
                  initial.predictor.method = initial.predictor.method,
                  tree.depth = tree.depth, n.trees.rf = n.trees.rf, n.trees.boosting = n.trees.boosting,
                  B = B, Kfold = Kfold, plot.gbmperf = plot.gbmperf,
                  error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune)

  if (initial.predictor.method == "boosting" & fit$best.iter == n.trees.boosting) {
    warning(paste("The best boosting iteration was iteration number", n.trees.boosting, " out of ", n.trees.boosting, ". Consider increasing the maximum number of trees and turning on boosting performance plot (plot.gbmperf = TRUE).", sep = ""))
  }

  ####### Construct the score in the whole dataset (no training or validation for pmsurv) ------------------------------------------
  fit.score <- scoresurv(fit = fit, x.cate = x.cate, tau0 = tau0, score.method = score.method)

  ####### Estimate the treatment effect in the whole dataset --------------------------------
  errors <- warnings <- c()
  if (prop[length(prop)] == 1){
    est.prop1 <- survCatch(drsurv(y = y, d = d, x.cate = x.cate, x.ps = x.ps, x.ipcw = x.ipcw, trt = trt, yf = yf,
                                  tau0 = tau0, surv.min = surv.min, ps.method = ps.method, minPS = minPS, maxPS = maxPS,
                                  ipcw.method = ipcw.method))
    est.prop1$ate.rmtl.high <- exp(est.prop1$log.rmtl.ratio)
    est.prop1$hr.high <- exp(est.prop1$log.hazard.ratio)
  } else {
    est.prop1 <- NULL
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
    cat(paste0('    Error(s) occurred when estimating the ATEs in the nested subgroup using "', paste0(errors, collapse = '", "'), '";\n    return NAs in the corresponding subgroup.'),'\n')
    warning(paste0('Error(s) occurred when estimating the ATEs in the nested subgroup using "', paste0(errors, collapse = '", "'), "\". NAs are returned for RMTL ratio and HR in the corresponding subgroup; see 'errors/warnings'."))
  }

  if (length(warnings) != 0) {
    cat(paste0('    Warning(s) occurred when estimating the ATEs in the nested subgroup using "', paste0(warnings, collapse = '", "'), '"'),'\n')
    warning(paste0('Warning(s) occurred when estimating the ATEs in the nested subgroup using "', paste0(warnings, collapse = '", "'), "\"; see 'errors/warnings'."))
  }

  if (sum(score.method %in% c("poisson", "twoReg", "contrastReg")) > 0) {
    cf <- data.frame(matrix(NA, nrow = ncol(x.cate) + 1, ncol = 3))
    colnames(cf) <- c("poisson", "twoReg", "contrastReg")
    rownames(cf) <- c("(Intercept)", colnames(x.cate))

    if ("poisson" %in% score.method) cf$poisson <- fit$result.poisson

    if("twoReg" %in% score.method) cf$twoReg <- fit$result.twoReg

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

  t.end <- Sys.time()
  t.diff <- difftime(t.end, t.start)
  cat('Total runtime :',as.numeric(t.diff), attributes(t.diff)$units, '\n')
  return(result)
}



#' Two side-by-side line plots of validation curves from the \code{"PrecMed"} object
#'
#' Provides validation curves in two side-by-side plots, visualizing the estimated ATEs in a series
#' of nested subgroups in the training set and validation set separately, where each line represents
#' one scoring method specified in \code{\link{cv}()} or \code{\link{cvmean}()}. This should be run
#' only after results of \code{\link{cv}()} or \code{\link{cvmean}()} have been obtained.
#'
#' @param x An object of class \code{"PrecMed"}.
#' @param cv.i A positive integer indicating the index of the CV iteration results to be plotted.
#' Allowed values are: a positive integer \eqn{<=} \code{cv.n} in \code{\link{cv}()} or
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
#' @param grayscale A logical value indicating grayscale plots (\code{TRUE}) or colored plots (\code{FALSE}).
#' Default is \code{FALSE}.
#' @param xlim A numeric value for the range of the x-axis. Default is \code{NULL}, which means there is no
#' range specified.
#' @param ... Other parameters
#'
#' @return Returns two side-by-side line plots, one of which shows the validation curves of the training
#' sets and the other the validation curves in the validation sets. A gray horizontal dashed line of
#' overall ATE is included as a reference. ABC statistics will be added to the legend if
#' \code{show.abc = TRUE}.
#'
#' @details \code{\link{plot}()} takes in outputs from \code{\link{cv}()} and generates two plots
#' of validation curves side-by-side, one for the training set and one for validation set.
#' Separate validation curves are produced for each scoring method specified via \code{score.method}
#' in \code{\link{cv}()} or \code{\link{cvmean}()}.
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
#' @seealso \code{\link{abc}()} and \code{\link{boxplot}()} for \code{"PrecMed"} objects.
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
#' # default setting
#' plot(x = cv_count)
#'
#' # turn off ABC annotation
#' plot(x = cv_count, show.abc = FALSE)
#'
#' # grayscale
#' plot(x = cv_count, grayscale = TRUE)
#'
#' # plot the validation curves from the 2nd CV iteration instead of the mean of all validation curves
#' plot(x = cv_count, cv.i = 2)
#'
#' # median of the validation curves
#' plot(x = cv_count, combine = "median")
#'
#' # plot validation curves in validation set only
#' plot(x = cv_count, valid.only = TRUE)
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
#' # default setting, plot RMTL ratios in both training and validation sets
#' plot(x = cv_surv)
#'
#' # plot hazard ratio
#' plot(x = cv_surv, plot.hr = TRUE)
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
#' plot(x = cv.mean)
#'}
#'
#' @export
#'
#' @importFrom rlang .data
#' @importFrom stringr str_extract
#' @importFrom MESS auc
#' @importFrom dplyr filter select mutate group_by ungroup
#' @importFrom ggplot2 ggplot aes facet_wrap geom_boxplot geom_hline geom_line labs scale_color_manual scale_linetype_manual scale_y_continuous theme theme_classic
#' @importFrom tidyr fill

plot.PrecMed <- function(x,
                         cv.i = NULL,
                         combine = "mean",
                         show.abc = TRUE,
                         valid.only = FALSE,
                         plot.hr = FALSE,
                         ylab = NULL,
                         legend.position = "bottom",
                         grayscale = FALSE,
                         xlim = NULL, ...) {

  # Must be count or survival data for this function
  stopifnot(x$response == "count" | x$response == "survival" | x$response == "continuous")

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
  if (plot.hr == TRUE & x$response != "survival") stop("Hazard ratio plot is only available for PrecMed object with x$response equals to 'survival'.")
  if (!(plot.hr %in% c(TRUE, FALSE))) stop("plot.hr has to be boolean.")

  # Retrieve proportion
  prop.onlyhigh <- round(x$props$prop.onlyhigh, 5)
  if(is.null(xlim) == TRUE) xlim <- c(min(prop.onlyhigh), max(prop.onlyhigh))
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
  if (is.null(ylab)){
    if (x$response == "count"){
      what.ratio <- "Rate ratio"
    } else if (x$response == "survival"){
      if (plot.hr == TRUE){
        what.ratio <- "Hazard ratio"
      } else {
        what.ratio <- "RMTL ratio"
      }
    } else if (x$response == "continuous"){
      ##TODO: what.ratio? correct name?
      what.ratio <- "Mean difference"
    }

    if (is.null(x$formulas$trt_labels) == TRUE & x$response == "continuous"){
      ylab <- paste0(what.ratio, " ", trt, "=1 - ", trt, "=0\nin each subgroup")

    } else if (is.null(x$formulas$trt_labels) == TRUE & x$response != "continuous") {
      ylab <- paste0(what.ratio, " ", trt, "=1 - ", trt, "=0\nin each subgroup")
    } else {
      ylab <- paste0(what.ratio, " ", trt, "=", x$formulas$trt_labels[2], " - ", trt, "=", x$formulas$trt_labels[1], "\nin each subgroup")
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

  for (i in 1:m) {
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


  # Plot
  if (grayscale == TRUE){
    map.color <- c("randomForest" = "#969696",
                   "boosting" = "#969696",
                   "poisson" = "#969696",
                   "gaussian" = "#969696",
                   "twoReg" = "black",
                   "contrastReg" = "black",
                   "negBin" = "#969696",
                   "gam" = "#969696")

    map.linetype <- c("randomForest" = "dotdash",
                      "boosting" = "dotted",
                      "poisson" = "solid",
                      "gaussian" = "solid", #Only one of poisson or gaussian is mapped depending on the data type
                      "twoReg" = "dotdash",
                      "contrastReg" = "solid",
                      "negBin" = "dotdash",
                      "gam" = "twodash")

    colors <- map.color[score.method]
    linetypes <- map.linetype[score.method]
    attr(colors, "names") <- attr(linetypes, "names") <- sapply(1:m, function(x) grep(mapped.method[x], scorevars, value = TRUE))

    p <- dat %>%
      ggplot(aes(x = .data$prop, y = .data$ate.cv, color = .data$score, linetype = .data$score)) +
      geom_line(size = 1.2) +
      facet_wrap(~ .data$set, ncol = 2) +
      scale_y_continuous(limits = c(lowery, uppery)) +
      labs(y = ylab, x = xlab, color = "Method", linetype = "Method") +
      theme_classic() +
      theme(legend.position = legend.position) +
      scale_linetype_manual(values = linetypes) +
      scale_color_manual(values = colors)

  } else {

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
      labs(y = ylab, x = xlab) +
      theme_classic() +
      theme(legend.position = legend.position) +
      scale_color_manual(name = "Method", values = colors)
  }

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


#' A set of box plots of estimated ATEs from the \code{"PrecMed"} object
#'
#' Provides box plots which depict distributions of estimated ATEs for each multi-category subgroup in
#' the validation set across all cross-validation iterations. The subgroups are mutually exclusive and
#' are categorized by the CATE score percentiles (\code{prop.multi} specified in \code{\link{cv}()} or
#' \code{\link{cvmean}()}). Box plots of mutually exclusive subgroups are constructed separately by scoring
#' method specified in \code{\link{cv}()}. This should be run only after results of \code{\link{cv}()} or
#' \code{\link{cvmean}()}) have been obtained.
#'
#' @param x An object of class \code{"PrecMed"}.
#' @param ylab A character value for the y-axis label to describe what the ATE is. Default is \code{NULL},
#' which creates a default y-axis label based on available data.
#' @param plot.hr A logical value indicating whether the hazard ratios should be plotted in the
#' validation curves (\code{TRUE}). Otherwise, the restricted mean time lost is plotted (\code{FALSE}).
#' This argument is only applicable to survival outcomes. Default is \code{FALSE}.
#' @param grayscale A logical value indicating grayscale plots (\code{TRUE}) or colored plots (\code{FALSE}).
#' Default is \code{FALSE}.
#' @param ... Other parameters
#'
#' @return Returns sets of box plots, one set for each scoring method, over each of the multi-category
#' subgroups. A gray horizontal dashed line of the overall ATE is included as a reference.
#'
#' @details \code{\link{boxplot}()} takes in outputs from \code{\link{cv}()} and generates
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
#' @seealso \code{\link{plot}} and \code{\link{abc}()} for \code{"PrecMed"} objects.
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
#' boxplot(x = cv_count, ylab = "Rate ratio of drug1 vs drug0 in each subgroup")
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
#' boxplot(x = cv_surv, ylab = "RMTL ratio of drug1 vs drug0 in each subgroup")
#'
#'# Continuous outcome
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
#' boxplot(x = cv.mean)
#'}
#'
#' @export
#'
#' @importFrom graphics boxplot
#' @importFrom dplyr mutate_at vars contains
#' @importFrom ggplot2 ggplot aes facet_wrap geom_boxplot geom_hline geom_line labs scale_color_manual scale_linetype_manual scale_y_continuous theme theme_classic
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data

# Note: additional arguments ggplot scale, legend.position, nrow.plot (what does this mean?)

boxplot.PrecMed <- function(x, ylab = NULL,
                            plot.hr = FALSE,
                            grayscale = FALSE, ...) {

  # Must be count or survival data for this function
  stopifnot(x$response == "count" | x$response == "survival" | x$response == "continuous")

  # Argument checks on plot.hr
  if (!(plot.hr %in% c(TRUE, FALSE))) stop("plot.hr has to be boolean.")
  if (plot.hr == TRUE & x$response != "survival") stop("Hazard ratio plot is only available for PrecMed objects with x$response is \"survival\".")

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
  if (is.null(ylab)){
    if (x$response == "count"){
      what.ratio = "Rate ratio"
    } else if (x$response == "survival"){
      if (plot.hr == TRUE){
        what.ratio = "Hazard ratio"
      } else {
        what.ratio = "RMTL ratio"
      }
    } else if (x$response == "continuous"){
      ##TODO: what.ratio? correct name?
      what.ratio = "Mean difference"
    }

    if (is.null(x$formulas$trt_labels) == TRUE & x$response == "continuous"){
      ylab <- paste0(what.ratio, " ", trt, "=1 - ", trt, "=0\nin each subgroup")

    } else if (is.null(x$formulas$trt_labels) == TRUE & x$response != "continuous") {
      ylab <- paste0(what.ratio, " ", trt, "=1 - ", trt, "=0\nin each subgroup")
    } else {
      ylab <- paste0(what.ratio, " ", trt, "=", x$formulas$trt_labels[2], " - ", trt, "=", x$formulas$trt_labels[1], "\nin each subgroup")
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
  if (plot.hr == FALSE) {
    prefix <- "ate."
  } else {
    prefix <- "hr."
  }

  for (i in 1:m) {
    ate <- x[[paste0(prefix, score.method[i])]]
    if (cv.n == 1) { # TODO: Temporary solutions as cbind not working when cv.n = 1
      results[((i - 1) * cv.n + 1):((i - 1) * cv.n + cv.n),] <- c(t(ate[[paste0(prefix, "est.valid.group.cv")]])[l:1],
                                                                  1,
                                                                  unname(mapped.method[i]))
    } else {
      results[((i - 1) * cv.n + 1):((i - 1) * cv.n + cv.n),] <- cbind(t(ate[[paste0(prefix, "est.valid.group.cv")]])[,l:1],
                                                                      1:cv.n,
                                                                      unname(mapped.method[i]))
    }
  }

  # Extract overall ATE
  if (plot.hr == FALSE) {
    overall.hline <- mean(x$overall.ate.valid)
  } else {
    overall.hline <- mean(x$overall.hr.valid)
  }

  results <- as.data.frame(results)
  colnames(results) <- c(1:l, "CV", "Method") # Subgroup 1 = lowest responder to drug1
  colnames(results)[1:l] <- subgroup_label <- paste(round(prop.multi[1:l] * 100), paste0(round(prop.multi[2:(l + 1)] * 100), "%"), sep = "-")

  if (grayscale == TRUE) {
    p <- results %>%
      mutate_at(c(1:l), as.numeric) %>%
      pivot_longer(-(.data$CV:.data$Method), names_to = "Subgroup", values_to = "ate") %>%
      mutate(Method = factor(.data$Method, levels = mapped.method),
             Subgroup = factor(.data$Subgroup, levels = subgroup_label)) %>%
      ggplot(aes(x = .data$Subgroup, y = .data$ate)) +
      geom_boxplot() +
      geom_hline(yintercept = overall.hline, linetype = "dashed", color= "grey70") +
      facet_wrap(~ .data$Method, nrow = 2, scales = "fixed") +
      labs(y = ylab, x = xlab) +
      theme_classic()
  } else {
    p <- results %>%
      mutate_at(c(1:l), as.numeric) %>%
      pivot_longer(-(.data$CV:.data$Method), names_to = "Subgroup", values_to = "ate") %>%
      mutate(Method = factor(.data$Method, levels = mapped.method),
             Subgroup = factor(.data$Subgroup, levels = subgroup_label)) %>%
      ggplot(aes(x = .data$Subgroup, y = .data$ate)) +
      geom_boxplot(aes(fill = .data$Subgroup)) +
      geom_hline(yintercept = overall.hline, linetype = "dashed", color = "grey70") +
      facet_wrap(~ .data$Method, nrow = 2, scales = "fixed") +
      ggplot2::scale_fill_hue(l = 45) +
      labs(y = ylab, x = xlab) +
      theme_classic() + theme(legend.position = "none")
  }

  return(p)
}



#' Doubly robust estimator of and inference for the average treatment effect for survival data
#'
#' Doubly robust estimator of the average treatment effect between two treatments, which is the
#' restricted mean time lost ratio for survival outcomes. Bootstrap is used for inference.
#'
#'
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
#' @return Return a list of 6 elements:
#' \itemize{
#'   \item{\code{rmst1}: } A vector of numeric values of the estimated RMST, bootstrap standard error,
#'   lower and upper limits of 95% confidence interval, and the p-value in the group \code{trt=1}.
#'   \item{\code{rmst0}: } A vector of numeric values of the estimated RMST, bootstrap standard error,
#'   lower and upper limits of 95% confidence interval, and the p-value in the group \code{trt=0}.
#'   \item{\code{log.rmtl.ratio}: } A vector of numeric values of the estimated log RMTL ratio of
#'   \code{trt=1} over \code{trt=0}, bootstrap standard error, lower and upper limits of 95% confidence
#'   interval, and the p-value.
#'   \item{\code{log.hazard.ratio}: } A vector of numeric values of the estimated adjusted log hazard ratio
#'   of \code{trt=1} over \code{trt=0}, bootstrap standard error, lower and upper limits of 95% confidence
#'   interval, and the p-value.
#'   \item{\code{warning}: } A warning message produced if the treatment variable was not coded as 0/1.
#'   The key to map the original coding of the variable to a 0/1 key is displayed in the warning to facilitate
#'   the interpretation of the remaining of the output.
#'   \item{\code{plot}: } If \code{plot.boot} is \code{TRUE}, a histogram displaying the distribution of the
#'   bootstrapped rmst1, rmst0, log.rmtl.ratio and log.hazard.ratio. The red vertical reference line in the
#'   histogram represents the estimates.
#' }
#'
#' @details This helper function estimates the average treatment effect (ATE) for survival data between two
#' treatment groups in a given dataset. The ATE is estimated with a doubly robust estimator that accounts for
#' imbalances in covariate distributions between the two treatment groups with inverse probability treatment and
#' censoring weighting. For survival outcomes, the estimated ATE is the estimated by RMTL ratio between treatment
#' 1 versus treatment 0. The log-transformed ATEs and log-transformed adjusted hazard ratios are returned, as well
#'  as the estimated RMST in either treatment group. The variability of the estimated RMTL ratio is calculated
#'  using bootstrap. Additional outputs include standard error of the log RMTL ratio, 95% confidence interval,
#'  p-value, and a histogram of the bootstrap estimates.
#'
#' @examples
#' tau0 <- with(survivalExample,
#'              min(quantile(y[trt == "drug1"], 0.95), quantile(y[trt == "drug0"], 0.95)))
#'
#' output <- drsurv.inference(cate.model = survival::Surv(y, d) ~ age +
#'                                                                female +
#'                                                                previous_cost +
#'                                                                previous_number_relapses,
#'                            ps.model = trt ~ age + previous_treatment,
#'                            data = survivalExample,
#'                            tau0 = tau0,
#'                            plot.boot = TRUE,
#'                            seed = 999)
#' print(output)
#' output$plot
#'
#' @export
#'
#' @importFrom ggplot2 ggplot geom_histogram geom_vline
#' @importFrom dplyr mutate
#' @importFrom tidyr gather

drsurv.inference <- function(cate.model, ps.model, data, ipcw.model = NULL, followup.time = NULL, tau0 = NULL, surv.min = 0.025, ipcw.method = "breslow",
                             ps.method = "glm", minPS = 0.01, maxPS = 0.99, n.boot = 500, seed = NULL, verbose = 1, plot.boot = FALSE) {

  # Set seed once for reproducibility
  set.seed(seed)

  #### CHECK ARGUMENTS ####
  arg.checks(fun = "drinf", response = "survival", data = data, followup.time = followup.time, tau0 = tau0, surv.min = surv.min,
             ps.method = ps.method, minPS = minPS, maxPS = maxPS, ipcw.method = ipcw.method, n.boot = n.boot, plot.boot = plot.boot)

  #### PRE-PROCESSING ####
  preproc <- data.preproc.surv(fun = "drinf", cate.model = cate.model, ps.model = ps.model, ipcw.model = ipcw.model, tau0 = tau0,
                               data = data, ps.method = ps.method, response = "survival")

  y <- preproc$y
  d <- preproc$d
  trt <- preproc$trt
  x.ps <- preproc$x.ps
  x.cate <- preproc$x.cate
  x.ipcw <- preproc$x.ipcw

  # Check if tau0 is large enough, i.e., if tau0 is larger than the 50% quantile of the observed or censoring time.
  if (tau0 < median(y)) warning("It is recommended to increase tau0 close to the largest observed or censored time.")

  if (is.null(followup.time)) {
    yf <- NULL
  } else {
    yf <- data[[followup.time]]
  }

  #### FUNCTION STARTS HERE ####
  # Point estimate
  est <- drsurv(y = y, d = d, x.cate = x.cate, x.ps = x.ps, x.ipcw = x.ipcw, trt = trt, yf = yf,
                tau0 = tau0, surv.min = surv.min, ps.method = ps.method, minPS = minPS, maxPS = maxPS, ipcw.method = ipcw.method)
  if (plot.boot == TRUE) {
    est.df <- data.frame(est) %>%
      gather() %>%
      mutate(key = factor(.data$key, levels = c("rmst1", "rmst0", "log.rmtl.ratio", "log.hazard.ratio")))
  }
  est <- unlist(est)

  # Apply bootstrap
  n <- length(y)
  save.boot <- matrix(NA, nrow = n.boot, ncol = length(est))
  colnames(save.boot) <- c("rmst1", "rmst0", "log.rmtl.ratio", "log.hazard.ratio")
  for (i in 1:n.boot) {
    idsub.boot <- sample(n, size = n, replace  = TRUE)
    est.boot <- drsurv(y = y[idsub.boot],
                       d = d[idsub.boot],
                       x.cate = x.cate[idsub.boot, ],
                       x.ps = x.ps[idsub.boot, ],
                       x.ipcw = x.ipcw[idsub.boot, ],
                       trt = trt[idsub.boot],
                       yf = yf[idsub.boot],
                       tau0 = tau0, surv.min = surv.min, ps.method = ps.method, minPS = minPS, maxPS = maxPS, ipcw.method = ipcw.method)
    save.boot[i,] <- unlist(est.boot)

    if (i %% (n.boot %/% 10) == 0) {
      if (verbose == 1) cat("Bootstrap iteration", i, "\n")
      if (plot.boot == TRUE) {
        p <- save.boot %>% as.data.frame() %>% gather(key = "key", value = "value") %>% na.omit() %>%
          mutate(key = factor(.data$key, levels = c("rmst1", "rmst0", "log.rmtl.ratio", "log.hazard.ratio"))) %>%
          ggplot(aes(x = .data$value)) + geom_histogram(bins = 50, alpha = 0.7) +
          facet_wrap(. ~ .data$key, scales = "free") +
          geom_vline(data = est.df, aes(xintercept = .data$value), linetype = 1, color = "red") +
          theme_classic() +
          labs(x = "Bootstrap values", y = "Frequency", title = paste0(i, " bootstrap iterations"))
        print(p)

      } # end of if (plot.boot == TRUE)
    } # end of if (i %% (n.boot %/% 10) == 0)
  } # end of for (i in 1:n.boot)

  se.est <- apply(save.boot, 2, sd)
  cil.est <- est - qnorm(0.975) * se.est
  ciu.est <- est + qnorm(0.975) * se.est
  p.est  <- 1 - pchisq(est^2 / se.est^2, df = 1)

  est.all <- data.frame(estimate = est, SE = se.est, CI.lower = cil.est, CI.upper = ciu.est, pvalue = p.est)

  out <- c()
  out$rmst1 <- data.frame(est.all[1, ])
  out$rmst0 <- data.frame(est.all[2, ])
  out$log.rmtl.ratio <- data.frame(est.all[3, ])
  out$log.hazard.ratio <- data.frame(est.all[4, ])
  out$warning <- preproc$warning

  if (plot.boot == TRUE) {
    plot <- save.boot %>% as.data.frame() %>% gather(key = "key", value = "value") %>% na.omit() %>%
      mutate(key = factor(.data$key, levels = c("rmst1", "rmst0", "log.rmtl.ratio", "log.hazard.ratio"))) %>%
      ggplot(aes(x = .data$value)) + geom_histogram(bins = 50, alpha = 0.7) +
      facet_wrap(. ~ .data$key, scales = "free") +
      geom_vline(data = est.df, aes(xintercept = .data$value), linetype = 1, color = "red") +
      theme_classic() +
      labs(x = "Bootstrap values", y = "Frequency", title = paste0(n.boot, " bootstrap iterations"))
    out$plot <- plot
  }

  return(out)
}


