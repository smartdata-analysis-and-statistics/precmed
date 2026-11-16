#' Estimate conditional treatment effect scores
#'
#' Estimate CATE scores for binary, count and survival outcomes. Binary scores
#' include pooled logistic regression, doubly robust two-regression estimation and
#' boosting. Count and survival methods include both outcome-regression and
#' doubly robust scores. Continuous outcomes are not yet implemented.
#'
#' @param response A string describing the type of outcome in the data. Allowed values include
#' "binary", "count" (see \code{\link{catecvcount}()}), "survival" (see \code{\link{catecvsurv}()}) and "continuous" (see \code{\link{catecvmean}()}).
#' @param data A data frame containing the variables in the outcome, treatment,
#' and censoring models (if specified); a data frame with \code{n} rows (1 row per observation).
#' @param score.method A vector of one or multiple methods to estimate the CATE score.
#' Allowed values are: \code{'boosting'}, \code{'twoReg'}, \code{'contrastReg'}, \code{'poisson'} (count and survival outcomes only),
#' \code{'randomForest'} (survival, continuous outcomes only), \code{negBin} (count outcomes only), \code{'gam'} (continuous outcomes only),
#' \code{'gaussian'} (continuous outcomes only). For binary outcomes, use
#' \code{'logistic'}, \code{'twoReg'}, or \code{'boosting'}.
#' @param cate.model A formula describing the outcome model to be fitted.
#' The outcome must appear on the left-hand side. For survival outcomes, a \code{Surv} object
#' must be used to describe the outcome.
#' @param ps.model A formula describing the treatment model (propensity-score model) to be fitted. The treatment must
#' appear on the left-hand side. For binary outcomes, logical treatment values
#' \code{FALSE} and \code{TRUE} are mapped to 0 and 1, respectively.
#' Numeric 0/1 values or two distinct treatment labels are also accepted.
#' If data are from a randomized controlled trial, specify \code{ps.model = trt ~ 1} as an intercept-only model.
#' @param ps.method A character value for the method to estimate the propensity score.
#' Allowed values include one of:
#' \code{'glm'} for logistic regression with main effects only (default), or
#' \code{'lasso'} for a logistic regression with main effects and LASSO penalization on
#' two-way interactions (added to the model if interactions are not specified in \code{ps.model}).
#' Relevant only when \code{ps.model} has more than one variable.
#' @param init.model A formula describing the initial predictor model. The outcome must appear on the left-hand side.
#' It must be specified when \code{score.method = contrastReg} or \code{twoReg}, except for binary outcomes, which use \code{cate.model}.
#' @param initial.predictor.method A character vector for the method used to get initial
#' outcome predictions conditional on the covariates specified in \code{cate.model}. Only applies
#' when \code{score.method} includes \code{'twoReg'} or \code{'contrastReg'}.Allowed values include one of
#' \code{'randomForest'} (survival outcomes only), \code{'boosting'}, \code{'logistic'}
#' (survival outcomes only, fast), \code{'poisson'} (count outcomes only, fast), \code{'gaussian'} (continuous outcomes only) and
#' \code{'gam'} (count and continuous outcomes only). Default is \code{NULL}, which assigns \code{'boosting'}
#' for count outcomes and \code{'randomForest'} for survival outcomes.
#' @param ipcw.model A formula describing the censoring model used for inverse
#' probability of censoring weighting (IPCW). The left-hand side must be empty. Only applies for survival outcomes.
#' Default is \code{NULL}, which uses a censoring model with the same covariates
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
#' @param rf.control Named list of randomForestSRC options: \code{mtry}, \code{nodesize},
#' \code{nsplit}, \code{nodedepth}, and \code{splitrule}. Applies to both treatment arms
#' and survival initial-predictor forests. Use \code{n.trees.rf} for the number of trees.
#' @param n.trees.rf A positive integer specifying the number of trees in random forest.
#' Used if \code{score.method = 'randomForest'} or if \code{initial.predictor.method = 'randomForest'}
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
#' @param ... Additional boosting arguments; see \code{\link{boosting-control}}.
#'
#' For \code{response = "binary"}, \code{score.method = "logistic"}
#' fits a logistic model with treatment interactions. \code{"twoReg"} cross-fits
#' separate logistic outcome models and a treatment model, then projects the
#' augmented inverse-propensity outcomes for each arm onto the covariates in
#' \code{cate.model} using linear regression. The score is the difference of
#' these projections, averaged over \code{B} repetitions of \code{Kfold} folds.
#' It estimates the linear projection of the conditional risk difference if
#' either the treatment model or both outcome models are correct, under
#' identification and positivity assumptions. Propensity score truncation can
#' affect this property. The projection need not lie between -1 and 1 and
#' equals the conditional risk difference only if that difference is linear in
#' the specified design. The difference between arm coefficients is returned in
#' \code{coefficients$twoReg} and by \code{coef()}; the separate arm coefficients
#' are stored in \code{fit$twoReg$coefficients}. Logistic model coefficients are
#' on the log-odds scale and are stored in \code{fit$logistic}; they are not
#' linear risk-difference score coefficients. Boosting does not have a linear
#' score coefficient vector. With neither \code{twoReg} nor another linear
#' binary score, \code{coef()} returns \code{NULL}.
#' \code{"boosting"} fits
#' separate Bernoulli gradient-boosting models in each treatment arm, with
#' the number of trees selected by five-fold cross-validation. All three scores
#' are predicted risk differences; subgroup effects are doubly robust risk
#' differences. Binary boosting requires at least one outcome-model covariate
#' and both outcome values in each arm. Its fitted models and selected tree
#' counts are stored in \code{fit$boosting}. Scores and subgroup effects are
#' estimated on the same data; they are not independent validation. For extra boosting controls,
#' see \code{\link{boosting-control}}.
#' Subgroup fits without both arms or with rank-deficient models return NA.
#' @return Fitted objects retain \code{treatment.labels} in treatment 0/1 order
#' for score plot labels. For count response, see description of other outputs in \code{\link{catefitcount}()}.
#' For survival response, see description of outputs in \code{\link{catefitsurv}()}.
#'
#' @details For count response, see details in \code{\link{catefitcount}()}.
#' For survival response, see details in \code{\link{catefitsurv}()}.
#'
#' @references Yadlowsky, S., Pellegrini, F., Lionetto, F., Braune, S., & Tian, L. (2020).
#' \emph{Estimation and validation of ratio-based conditional average treatment effects using
#' observational data. Journal of the American Statistical Association, 1-18.} DOI: 10.1080/01621459.2020.1772080.
#'
#' Hoogland, J., Debray, T. P. A., Desachy, G., & Nguyen, T.-L. (2026).
#' Introduction to methods for personalizing medicine. In T. P. A. Debray,
#' T.-L. Nguyen, & R. W. Platt (Eds.),
#' \emph{Comparative Effectiveness and Personalized Medicine Research Using Real-World Data}
#' (Chapter 17, pp. 421--463). Chapman & Hall/CRC.
#' \doi{10.1201/9781003300809-17}
#'
#' Hoogland, J., IntHout, J., Belias, M., et al. (2021). A tutorial on
#' individualized treatment effect prediction from randomized trials with a
#' binary endpoint. \emph{Statistics in Medicine}, 40(26), 5961--5981.
#' \doi{10.1002/sim.9154}
#'
#' Jiang, X., Simoneau, G., Pellegrini, F., Shen, C., de Moor, C., & Tian, L. (2026).
#' Visualization and interpretation of individualized treatment rule results. In
#' T. P. A. Debray, T.-L. Nguyen, & R. W. Platt (Eds.),
#' \emph{Comparative Effectiveness and Personalized Medicine Research Using Real-World Data}
#' (Chapter 20, pp. 513--536). Chapman & Hall/CRC.
#' \doi{10.1201/9781003300809-20}
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
#' \donttest{
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
#'                  tau0 = tau0, higher.y = TRUE, seed = 999, n.cores = 1)
#'
#' coef(fit_2)
#'
#' }
#' @param contrast.control A named list of options for count or survival contrast
#' regression. Supported names are \code{solver} (\code{"C"}, the default, or
#' \code{"R"}), \code{error.maxNR}, \code{max.iterNR}, \code{tune}, and
#' \code{optim.control} (a named list passed to \code{stats::optim} on fallback).
#' List entries override the corresponding top-level NR arguments. Options are
#' kept separate from \code{...}, which is forwarded to boosting. Count and survival
#' outcomes support this list; continuous outcomes retain their existing solver.
#' @param rmst.control Named options for Cox restricted mean survival integration.
#' Supported entry: \code{solver = "C"} (default) or \code{solver = "R"}.
#' The C backend accumulates each patient's integral without allocating a
#' patient-by-event-time matrix. Cox model fitting is unchanged. Available only
#' for survival outcomes; independent of \code{contrast.control}.
#' @param ps.control Named options for LASSO propensity-score interaction screening.
#' Supported entry: \code{solver = "C"} (default) or \code{solver = "R"}.
#' Only affects screening when \code{ps.method = "lasso"}; the LASSO fit and
#' interaction selection rule are unchanged. Kept separate from boosting's
#' \code{...} and other numerical controls.
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
                    verbose = 0,
                    contrast.control = list(),
                    rmst.control = list(), ps.control = list(), rf.control = list(), ...) {
  check_response(response)
  check_score_method(response, score.method)
  resolve_ps_control(ps.control)
  resolve_rf_control(rf.control)
  if (length(rf.control) && response == "count") stop("rf.control is supported only for survival or continuous outcomes.")


  if (length(contrast.control) && response == "continuous")
    stop("contrast.control is supported only for count or survival outcomes.")
  if (length(rmst.control) && response != "survival")
    stop("rmst.control is currently supported only for survival outcomes.")
  .args <- as.list(match.call())[-1]
  .args$response <- NULL
  if (response == "count") .args$rf.control <- NULL
  if (response != "survival") .args$rmst.control <- NULL
  if (response == "continuous") .args$contrast.control <- NULL
  switch(response,
         binary = {
           fitout <- do.call(catefitbinary, .args, envir = parent.frame())
         },
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
#' @details Display subgroup treatment effects and available score coefficients.
#' Binary effects and \code{twoReg} coefficients are on the risk-difference
#' scale; count and survival score coefficients are on log-ratio scales.
#'
#' @return No return value
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
             "negBin", "gam", "gaussian", "logistic")

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

  if (identical(x$response, "survival")) {
    chr.methods <- evalm[vapply(evalm, function(method) !is.null(x[[paste0("chr.", method)]]), logical(1))]
    if (length(chr.methods)) {
      chr.out <- do.call(rbind, lapply(chr.methods, function(method) x[[paste0("chr.", method)]]))
      rownames(chr.out) <- chr.methods
      cat("Cumulative hazard ratio in the nested subgroups:\n\n")
      print(chr.out)
      cat("\n")
    }
  }
  if (!identical(x$response, "binary") || !is.null(coef(x))) print(coef(x))

  if (!is.null(x$warning)) {
    cat("\n")
    message(x$warning)
  }
}



