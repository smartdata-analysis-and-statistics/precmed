#' Estimation of the conditional average treatment effect (CATE) score for count, survival and continuous data
#'
#' Provides singly robust and doubly robust estimation of CATE score for count, survival and continuous data with
#' up the following scoring methods among the following: Random forest (survival, continuous only), boosting,
#' poisson regression (count, survival only), two regressions, contrast regression, negative binomial regression (count only),
#' linear regression (continuous only), and generalized additive model (continuous only).
#'
#' @param response A string describing the type of outcome in the data. Allowed values include
#' "count" (see \code{\link{catecvcount}()}), "survival" (see \code{\link{catecvsurv}()}) and "continuous" (see \code{\link{catecvmean}()}).
#' @param data A data frame containing the variables in the outcome, propensity score, and inverse
#' probability of censoring models (if specified); a data frame with \code{n} rows (1 row per observation).
#' @param score.method A vector of one or multiple methods to estimate the CATE score.
#' Allowed values are: \code{'boosting'}, \code{'twoReg'}, \code{'contrastReg'}, \code{'poisson'} (count and survival outcomes only),
#' \code{'randomForest'} (survival, continuous outcomes only), \code{negBin} (count outcomes only), \code{'gam'} (continuous outcomes only),
#' \code{'gaussian'} (continuous outcomes only).
#' @param cate.model A formula describing the outcome model to be fitted.
#' The outcome must appear on the left-hand side. For survival outcomes, a \code{Surv} object
#' must be used to describe the outcome.
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
#' when \code{score.method} includes \code{'twoReg'} or \code{'contrastReg'}.Allowed values include one of
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
#' @param minPS A numerical value (in [0, 1]) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in (0, 1]) above which estimated propensity scores should be
#' truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.
#' @param followup.time A column name in \code{data} specifying the maximum follow-up time,
#' interpreted as the potential censoring time. Only applies for survival outcomes.
#' Default is \code{NULL}, which corresponds to unknown potential censoring time.
#' @param tau0 The truncation time for defining restricted mean time lost. Only applies for
#' survival outcomes. Default is \code{NULL}, which corresponds to setting the truncation time as the
#' maximum survival time in the data.
#' @param higher.y A logical value indicating whether higher (\code{TRUE}) or
#' lower (\code{FALSE}) values of the outcome are more desirable. Default is \code{TRUE}.
#' @param prop.cutoff A vector of numerical values (in (0, 1]) specifying percentiles of the
#' estimated log CATE scores to define nested subgroups. Each element represents the cutoff to
#' separate observations in nested subgroups (below vs above cutoff).
#' The length of \code{prop.cutoff} is the number of nested subgroups.
#' An equally-spaced sequence of proportions ending with 1 is recommended.
#' Default is \code{seq(0.5, 1, length = 6)}.
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
#' @param verbose An integer value indicating whether intermediate progress messages and histograms should
#' be printed. \code{1} indicates messages are printed and \code{0} otherwise. Default is \code{0}.
#'
#' @return For count response, see description of outputs in \code{\link{catefitcount}()}.
#' For survival response, see description of outputs in \code{\link{catefitsurv}()}.
#'
#' @details For count response, see details in \code{\link{catefitcount}()}.
#' For survival response, see details in \code{\link{catefitsurv}()}.
#'
#' @references Yadlowsky, S., Pellegrini, F., Lionetto, F., Braune, S., & Tian, L. (2020).
#' \emph{Estimation and validation of ratio-based conditional average treatment effects using
#' observational data. Journal of the American Statistical Association, 1-18.}
#' \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1772080}
#'
#' @seealso \code{\link{catecv}()}
#'
#' @examples
#' # Count outcome
#' fit_1 <- catefit(response = "count",
#'                  data = countExample,
#'                  score.method = "poisson",
#'                  cate.model = y ~ age + female + previous_treatment +
#'                                   previous_cost + previous_number_relapses +
#'                                   offset(log(years)),
#'                  ps.model = trt ~ age + previous_treatment,
#'                  higher.y = TRUE,
#'                  seed = 999)
#'
#' coef(fit_1)
#'
#' \dontrun{
#' # Survival outcome
#' library(survival)
#' tau0 <- with(survivalExample,
#'                  min(quantile(y[trt == "drug1"], 0.95), quantile(y[trt == "drug0"], 0.95)))
#'
#' fit_2 <- catefit(response = "survival",
#'                  data = survivalExample,
#'                  score.method = c("poisson", "boosting", "randomForest"),
#'                  cate.model = Surv(y, d) ~ age + female + previous_cost +
#'                                            previous_number_relapses,
#'                  ps.model = trt ~ age + previous_treatment,
#'                  initial.predictor.method = "logistic",
#'                  ipcw.model = ~ age + previous_cost + previous_treatment,
#'                  tau0 = tau0,
#'                  higher.y = TRUE,
#'                  seed = 999)
#'
#' coef(fit_2)
#'
#'
#' # Continuous outcome
#' fit_3 <- catefit(response = "continuous",
#'                  data = meanExample,
#'                  score.method = c("gaussian", "randomForest", "twoReg", "contrastReg"),
#'                  cate.model = y ~ age + previous_treatment + previous_cost +
#'                                   previous_status_measure,
#'                  ps.model = trt ~ age + previous_treatment,
#'                  init.model = y ~ age + previous_treatment + previous_cost +
#'                                   previous_status_measure,
#'                  initial.predictor.method = "gaussian",
#'                  higher.y = FALSE,
#'                  seed = 999)
#'
#' coef(fit_3)
#'
#' }
#' @export

catefit <- function(response,
                    data,
                    score.method,
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
                    tau0 = NULL,
                    higher.y = TRUE,
                    prop.cutoff = seq(0.5, 1, length = 6),
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
           fitout <- do.call(catefitcount, .args)
         },
         survival = {
           fitout <- do.call(catefitsurv, .args)
         },
         continuous = {
           fitout <- do.call(catefitmean, .args)
         }
  )
  return(fitout)
}


#' Print function for atefit
#'
#' @param x An object of class \code{"catefit"}.
#' @param ... Other parameters
#'
#' @details Display the estimated treatment effects for survival outcomes (log
#' restricted mean time lost ratio and log hazard ratio) and count outcomes
#' (the log rate ratio).
#'
#' @value No return value
#'
#' @author Thomas Debray
#'
#' @export
#'
#' @importFrom dplyr %>%
#'
print.catefit <- function(x, ...) {

  cat("Estimated ATE in the nested subgroups:\n\n ")

  evalm <- c("boosting", "twoReg", "contrastReg", "poisson", "randomForest",
             "negBin", "gam", "gaussian")

  snames <- c()
  out <- NULL

  for (evalm_i in evalm) {
    if (evalm_i %in% x$score.method) {
      out <- rbind(out, x[[paste0("ate.",evalm_i)]])
      snames <- c(snames, evalm_i)
    }
  }

  rownames(out) <- snames

  print(out)
  cat("\n")

  print(coef(x))

  if (!is.null(x$warning)) {
    cat("\n")
    warning(x$warning)
  }
}



