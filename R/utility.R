#' Check arguments that are common to all types of outcome
#' USed inside \code{arg.checks()}
#'
#' @param fun A function for which argument check is needed; "pm" for \code{pmcount()}, "cv" for \code{cvcount()},
#' and "drinf" for \code{drcount.inference()}. No default.
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
#' @param higher.y A logical value indicating whether higher (\code{TRUE}) or lower (\code{FALSE})
#' values of the outcome are more desirable. Default is \code{TRUE}.
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
#' @param B A positive integer specifying the number of time cross-fitting is repeated in
#' \code{score.method = 'twoReg'} and \code{'contrastReg'}. Default is \code{3}.
#' @param Kfold A positive integer specifying the number of folds (parts) used in cross-fitting
#' to partition the data in \code{score.method = 'twoReg'} and \code{'contrastReg'}.
#' Default is \code{6}.
#' @param plot.gbmperf A logical value indicating whether to plot the performance measures in
#' boosting. Used only if \code{score.method = 'boosting'} or if \code{score.method = 'twoReg'}
#' or \code{'contrastReg'} and \code{initial.predictor.method = 'boosting'}. Default is \code{TRUE}.
#' @param tree.depth A positive integer specifying the depth of individual trees in boosting
#' (usually 2-3). Used only if \code{score.method = 'boosting'} or
#' if \code{score.method = 'twoReg'} or \code{'contrastReg'} and
#' \code{initial.predictor.method = 'boosting'}. Default is \code{2}.
#' @param n.trees.boosting A positive integer specifying the maximum number of trees in boosting
#' (usually 100-1000). Used only if \code{score.method = 'boosting'} or
#' if \code{score.method = 'twoReg'} or \code{'contrastReg'} and
#' \code{initial.predictor.method = 'boosting'}. Default is \code{200}.
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
#' @param n.boot A numeric value indicating the number of bootstrap samples used. This is only relevant
#' if \code{inference = TRUE}. Default is \code{500}.
#' @param plot.boot A logic value indicating whether histograms of the bootstrapped log(rate ratio) should
#' be produced at every \code{n.boot/10}-th iteration and whether the final histogram should be outputted.
#' Default is \code{FALSE}.
#'
#' @return Nothing. Will stop if arguments are incorrect.

arg.checks.common <- function(fun,
                              ps.method, minPS, maxPS,
                              higher.y = NULL, abc = NULL,
                              prop.cutoff = NULL, prop.multi = NULL,
                              B = NULL, Kfold = NULL, plot.gbmperf = NULL,
                              tree.depth = NULL, n.trees.boosting = NULL,
                              error.maxNR = NULL, max.iterNR = NULL, tune = NULL,
                              train.prop = NULL, cv.n = NULL,
                              error.max = NULL, max.iter = NULL,
                              n.boot = NULL, plot.boot = NULL) {


  #### Check common arguments ####
  # Check values of ps.method
  if (!(ps.method %in% c("glm", "lasso"))) stop("ps.method must be either 'glm' or 'lasso'.")
  if (!(fun %in% c("catefit", "crossv", "drinf"))) stop("fun must be either 'catefit', 'crossv' or 'drinf'.")

  # Check values of minPS, maxPS
  if (minPS > 1 | minPS < 0) stop("minPS must be a number between 0 and 1, inclusive.")
  if (maxPS > 1 | maxPS < 0) stop("maxPS must be a number between 0 and 1, inclusive.")
  if (minPS >= maxPS) stop("minPS must be smaller than maxPS.")


  if (fun %in% c("crossv", "catefit")) {
    # Check values of higher.y
    if (!(higher.y %in% c(TRUE, FALSE))) stop("higher.y has to be boolean.")

    # Check values of prop
    if (any(prop.cutoff > 1) | any(prop.cutoff < 0)) stop("prop.cutoff must be > 0 and <= 1.")

    # Check cross-fitting arguments
    if (any(c(B, Kfold, tree.depth, n.trees.boosting) <= 0)) stop("B, Kfold, tree.depth, or n.trees.boosting must be > 0.")
    if (!(plot.gbmperf %in% c(TRUE, FALSE))) stop("plot.gbmperf has to be boolean.")

    # Check control arguments
    if (length(tune) != 2) stop("tune must be a vector with two elements.")
    if (any(c(error.maxNR, max.iterNR, tune[1], tune[2]) <= 0)) stop("error.maxNR, max.iterNR and tune must be > 0.")

    if (fun == "crossv") {
      # Check values of abc
      if (!(abc %in% c(TRUE, FALSE))) stop("abc has to be boolean.")

      # Check values of prop.multi
      if (any(prop.multi > 1) | any(prop.multi < 0)) stop("prop.multi must be >= 0 and <= 1.")
      if (length(prop.multi) == 2 & sum(prop.multi == c(0, 1)) == 2) stop("prop.multi must defined at least 2 subgroups.")

      # Check values of CV
      if (train.prop >= 1 | train.prop <= 0) stop("train.prop must be a number between 0 and 1, exclusive.")

      # Check if cv.n is a positive integer
      check_positive_integer(cv.n, var_name = "cv.n")


      # Check control values for balance.split
      if (any(c(error.max, max.iter) <= 0)) stop("error.max and max.iterNR must be > 0.")
    }
  } else if (fun == "drinf") {
    # Check if n.boot is a positive integer
    check_positive_integer(n.boot, var_name = "n.boot")

    # Check if plot.boot is a boolean
    check_boolean(plot.boot, "plot.boot")
  }
}


#' Check arguments
#' Catered to all types of outcome
#' Apply at the beginning of \code{pmcount()}, \code{cvcount()}, \code{drcount.inference()}, \code{catefitsurv()}, \code{catecvsurv()}, and \code{drsurv.inference()}
#'
#' @param fun A function for which argument check is needed; "catefit" for \code{catefitcount()} and \code{catefitsurv()}, "crossv" for \code{catecvcount()} and \code{catecvsurv()},
#' and "drinf" for \code{drcount.inference()} and \code{drsurv.inference()}. No default.
#' @param response The type of response. Always 'survival' for this function.
#' @param data A data frame containing the variables in the outcome and propensity score models;
#' a data frame with \code{n} rows (1 row per observation).
#' @param followup.time Follow-up time, interpreted as the potential censoring time. If the potential censoring time is known,
#' followup.time is the name of a corresponding column in the data. Otherwise, set \code{followup.time == NULL}.
#' @param tau0 The truncation time for defining restricted mean time lost.
#' @param surv.min Lower truncation limit for probability of being censored (positive and very close to 0).
#' @param ipcw.method The censoring model. Allowed values are: \code{'breslow'} (Cox regression with Breslow estimator of the baseline survivor function),
#' \code{'aft (exponential)'}, \code{'aft (weibull)'}, \code{'aft (lognormal)'} or \code{'aft (loglogistic)'}. Default is \code{'breslow'}.
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
#' @param higher.y A logical value indicating whether higher (\code{TRUE}) or lower (\code{FALSE})
#' values of the outcome are more desirable. Default is \code{TRUE}.
#' @param score.method A vector of one or multiple methods to estimate the CATE score.
#' Allowed values are: \code{'boosting'}, \code{'poisson'}, \code{'twoReg'}, \code{'contrastReg'},
#' \code{'negBin'}. Default specifies all 5 methods.
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
#' one of \code{'randomForest'}, \code{'boosting'} and \code{'logistic'} (fastest). Default is \code{'randomForest'}.
#' @param tree.depth A positive integer specifying the depth of individual trees in boosting
#' (usually 2-3). Used only if \code{score.method = 'boosting'} or
#' if \code{score.method = 'twoReg'} or \code{'contrastReg'} and
#' \code{initial.predictor.method = 'boosting'}. Default is \code{2}.
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
#' Default is \code{6}.
#' @param plot.gbmperf A logical value indicating whether to plot the performance measures in
#' boosting. Used only if \code{score.method = 'boosting'} or if \code{score.method = 'twoReg'}
#' or \code{'contrastReg'} and \code{initial.predictor.method = 'boosting'}. Default is \code{TRUE}.
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
#' @param n.boot A numeric value indicating the number of bootstrap samples used. This is only relevant
#' if \code{inference = TRUE}. Default is \code{500}.
#' @param plot.boot A logic value indicating whether histograms of the bootstrapped log(rate ratio) should
#' be produced at every \code{n.boot/10}-th iteration and whether the final histogram should be outputted.
#' Default is \code{FALSE}.
#' @param interactions A logical value indicating whether the outcome model should assume interactions
#' between \code{x} and \code{trt}. If \code{TRUE}, interactions will be assumed only if at least 10 patients
#' received each treatment option. Default is \code{TRUE}.
#'
#' @return Nothing. Will stop if arguments are incorrect.

arg.checks <- function(fun, response, data,
                       followup.time = NULL,
                       tau0 = NULL, surv.min = NULL,
                       ipcw.method = NULL,
                       ps.method, minPS, maxPS,
                       higher.y = NULL, score.method = NULL, abc = NULL,
                       prop.cutoff = NULL, prop.multi = NULL,
                       train.prop = NULL, cv.n = NULL,
                       error.max = NULL, max.iter = NULL,
                       initial.predictor.method = NULL,
                       tree.depth = NULL, n.trees.rf = NULL, n.trees.boosting = NULL,
                       B = NULL, Kfold = NULL, plot.gbmperf = NULL,
                       error.maxNR = NULL, max.iterNR = NULL, tune = NULL,
                       n.boot = NULL, plot.boot = NULL, interactions = NULL){

  #### Check common arguments ####
  arg.checks.common(fun = fun,
                    ps.method = ps.method, minPS = minPS, maxPS = maxPS,
                    higher.y = higher.y,
                    prop.cutoff = prop.cutoff, prop.multi = prop.multi,
                    B = B, Kfold = Kfold, plot.gbmperf = plot.gbmperf,
                    tree.depth = tree.depth, n.trees.boosting = n.trees.boosting,
                    error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune,
                    abc = abc,
                    train.prop = train.prop, cv.n = cv.n,
                    error.max = error.max, max.iter = max.iter,
                    n.boot = n.boot, plot.boot = plot.boot)

  #### Check argument only used for count outcome ####
  if (response == "count") {
    # Check initial predictor method
    if (is.null(initial.predictor.method) == FALSE) {
      if (!(initial.predictor.method %in% c("poisson", "boosting", "gam")))
        stop("initial.predictor.method must be 'poisson', 'boosting' or 'gam'.")
    }

    # Check values of score.method
    if (any(!(score.method %in% c("boosting", "poisson", "twoReg", "contrastReg", "negBin"))))
      stop("Elements of score.method must come from: 'boosting', 'poisson', 'twoReg', 'contrastReg', 'negBin'.")

    if (fun == "drinf") {
      check_boolean(interactions, "interactions")
    }
  }

  #### Check argument only used for survival outcome ####
  if (response == "survival") {

    # Check if followup.time is either NULL or a valid column name in the data
    if (!is.null(followup.time) && !as.character(followup.time) %in% names(data)) {
      stop("followup.time must be either NULL or a valid column name in the data.")
    }

    # Check if tau0 > 0
    if (!is.null(tau0) && tau0 <= 0) {
      stop("tau0 must be positive.")
    }

    # Check if surv.min is positive and within the acceptable range (0, 0.1)
    if (surv.min <= 0 || surv.min >= 0.1)
      stop("surv.min must be positive and less than 0.1.")


    # Check n.trees.rf
    if (any(c(n.trees.rf) <= 0)) stop("n.trees.rf must be > 0.")

    # Check values for ipcw.method
    if (!(ipcw.method %in% c("breslow", "aft (exponential)", "aft (weibull)", "aft (lognormal)", "aft (loglogistic)"))) stop("ipcw.method must be either 'breslow', 'aft (exponential)', 'aft (weibull)', 'aft (lognormal)', 'aft (loglogistic)'.")

    # Check initial predictor method
    if (is.null(initial.predictor.method) == FALSE) {
      if (!(initial.predictor.method %in% c("randomForest", "boosting", "logistic"))) stop("initial.predictor.method must be 'randomForest', 'boosting' or 'logistic'.")
    }

    # Check values of score.method
    if (any(!(score.method %in% c("boosting", "poisson", "twoReg", "contrastReg", "randomForest")))) stop("Elements of score.method must come from: 'boosting', 'poisson', 'twoReg', 'contrastReg', 'randomForest'.")

  }
}

#' Data preprocessing
#' Apply at the beginning of \code{pmcount()} and \code{cvcount()}, after \code{arg.checks()}
#'
#' @param fun A function for which argument check is needed; "pm" for \code{pmcount()}, "cv" for \code{cvcount()},
#' and "drinf" for \code{drcount.inference()}. No default.
#' @param cate.model A formula describing the outcome model to be fitted.
#' The outcome must appear on the left-hand side.
#' @param ps.model A formula describing the propensity score model to be fitted.
#' The treatment must appear on the left-hand side. The treatment must be a numeric vector
#' coded as 0/1. If data are from a RCT, specify \code{ps.model} as an intercept-only model.
#' @param data A data frame containing the variables in the outcome and propensity score models;
#' a data frame with \code{n} rows (1 row per observation).
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
#' @param initial.predictor.method A character vector for the method used to get initial
#' outcome predictions conditional on the covariates. Only applies when \code{score.method}
#' includes \code{'twoReg'} or \code{'contrastReg'}. Allowed values include one of
#' \code{'boosting'}, \code{'poisson'} (fast), and \code{'gam'}. Default is \code{NULL}, which assigns
#' \code{'boosting'} for count outcomes.
#'
#' @return A list of 6 elements:
#'            - y: outcome; vector of length \code{n} (observations)
#'            - trt: binary treatment; vector of length \code{n}
#'            - x.ps: matrix of \code{p.ps} baseline covariates (plus intercept); dimension \code{n} by \code{p.ps + 1}
#'            - x.cate: matrix of \code{p.cate} baseline covariates; dimension \code{n} by \code{p.cate}
#'            - time: offset; vector of length \code{n}
#'            - if \code{fun = "pm"}:
#'                - prop: formatted \code{prop.cutoff}
#'            - if \code{fun = "cv"}
#'                - prop.onlyhigh: formatted \code{prop.cutoff} with 0 removed if applicable
#'                - prop.bi; formatted \code{prop.cutoff} with 0 and 1 removed if applicable
#'                - prop.multi: formatted \code{prop.multi}, starting with 0 and ending with 1

data.preproc <- function(fun, cate.model, ps.model, data, prop.cutoff = NULL,
                         prop.multi = NULL, ps.method, initial.predictor.method = NULL) {

  ## cate.model and ps.model as formulas
  cate.model <- as.formula(cate.model)
  ps.model <- as.formula(ps.model)

  ## Extraction step: extract y, trt, x.cate, x.ps, time in matrix form from cate.model and ps.model
  cate.mat <- model.frame(cate.model, data, na.action = 'na.pass')
  ps.mat <- model.frame(ps.model, data, na.action = 'na.pass')

  # y
  y <- model.response(cate.mat)
  if (is.null(y) == TRUE) stop("Outcome must be supplied on the left-hand side of the cate.model formula.")

  # trt
  trt <- model.response(ps.mat)
  if (is.null(trt) == TRUE) stop("Treatment must be supplied on the left-hand side of the ps.model formula.")

  # Covariate matrices
  x.cate <- model.matrix(cate.model, cate.mat)[, -1, drop = FALSE]
  x.ps <- model.matrix(ps.model, ps.mat)
  if (ncol(x.ps) == 1 & ps.method == "lasso") stop("LASSO penalization irrelevant when ps.model specified as a function of an intercept only. Consider setting ps.method='glm'.")

  # time
  time <- model.offset(cate.mat)
  if (is.null(time) == TRUE) { # Eventually some outcomes will not need an offset
    time <- rep(0, nrow(data))
    warning("No offset supplied. Offset set to 0.")
  }

  ## Check missing data
  if (any(is.na(y)) | any(is.na(trt)) | any(is.na(x.cate)) | any(is.na(x.ps)) | any(is.na(time))) stop("Missing data not allowed in cate.model, ps.model or offset.")

  ## Check negative count
  if (any(y < 0)) stop("Negative y values not allowed for count outcomes.")

  ## Check treatment binary numeric and coded as 0/1
  cat <- warn.dr <- NULL
  if (as.numeric(length(unique(trt))) != 2) {
    stop("trt should describe two distinct treatments.")
  } else if (!all(unique(trt) %in% c(0, 1))) {
    cat <- sort(unique(trt))
    trt <- ifelse(trt == cat[1], 0, 1)
    warn.dr <- warning(paste0("Variable trt was recoded to 0/1 with ", cat[1], "->0 and ", cat[2], "->1.\n"))
  } else if (is.factor(trt) == TRUE) {
    trt <- as.numeric(trt == 1)
  }

  ## Assign default to initial.predictor.method if NULL
  if (is.null(initial.predictor.method) == TRUE & fun %in% c("pm", "cv")) initial.predictor.method <- "boosting"

  if (fun == "catefit") {
    ## Check values of prop
    prop <- sort(prop.cutoff) # sort the proportions from small to large
    if (prop[1] == 0) {
      prop <- prop[-1] # if first element is 0, remove it because this means we leave out 0% of individuals
      warning("The first element of prop.cutoff cannot be 0 and has been removed.")
    }

    return(list(y = y, trt = trt, x.cate = x.cate, x.ps = x.ps, time = time, prop = prop))
  } else if (fun == "crossv") {
    ## Check values of prop.cutoff and prop.multi
    prop.onlyhigh <- prop.bi <- sort(prop.cutoff)
    prop.multi <- sort(prop.multi)
    if (prop.onlyhigh[1] == 0) {
      prop.onlyhigh <- prop.onlyhigh[-1]
      warning("The first element of prop.cutoff cannot be 0 and has been removed.") # only need to show once so no warning in the prop.bi check below
    }
    if (prop.bi[1] == 0) {
      prop.bi <- prop.bi[-1]
    }
    if (prop.bi[length(prop.bi)] == 1) {
      prop.bi <- prop.bi[-length(prop.bi)]
    }
    if (prop.multi[1] != 0) {
      prop.multi <- c(0, prop.multi)
      warning("The first element of prop.multi must be 0 and has been added.")
    }
    if (prop.multi[length(prop.multi)] != 1) {
      prop.multi <- c(prop.multi, 1)
      warning("The last element of prop.multi must be 1 and has been added.")
    }

    return(list(y = y, trt = trt, x.cate = x.cate, x.ps = x.ps, time = time,
                prop.onlyhigh = prop.onlyhigh, prop.bi = prop.bi, prop.multi = prop.multi, cat.trt = cat, initial.predictor.method = initial.predictor.method))
  } else if (fun == "drinf") {
    return(list(y = y, trt = trt, x.cate = x.cate, x.ps = x.ps, time = time, warning = warn.dr))
  }
}


#' Compute the area under the curve using linear or natural spline interpolation
#'
#' This function computes the area under the curve for two vectors where one
#' corresponds to the x values and the other corresponds to the y values.
#' It supports both linear and spline interpolation.
#'
#'
#' @param x A numeric vector of x values.
#' @param y A numeric vector of y values of the same length as x.
#' @param from The value from where to start calculating the area under the curve.
#' Defaults to the smallest x value.
#' @param to The value from where to end the calculation of the area under the curve.
#' Defaults to the greatest x value.
#' @param type The type of interpolation: "linear" or "spline".
#' Defaults to "linear".
#' @param subdivisions An integer indicating how many subdivisions to use for `integrate`
#' (for spline-based approximations).
#' @param ... Additional arguments passed on to `approx` (for linear interpolations).
#' @return A numeric value representing the area under the curve.
#'
#' @importFrom stats approx integrate splinefun
#' @export
auc <- function(x, y, from = min(x, na.rm = TRUE), to = max(x, na.rm = TRUE),
                type = c("linear", "spline"), subdivisions = 100, ...)
{
  type <- match.arg(type)
  stopifnot(length(x) == length(y))
  stopifnot(!is.na(from))
  if (length(unique(x)) < 2)
    return(NA)
  if (type == "linear") {
    # Use approx for linear interpolation
    values <- approx(x, y, xout = sort(unique(c(from, to, x[x > from & x < to]))), ...)
    res <- 0.5 * sum(diff(values$x) * (values$y[-1] + values$y[-length(values$y)]))
  } else {
    # Use spline interpolation
    myfunction <- splinefun(x, y, method = "natural")
    res <- integrate(myfunction, lower = from, upper = to, subdivisions = subdivisions)$value
  }

  return(res)
}

check_positive_integer <- function(x, var_name = "value") {
  if (!is.numeric(x) || x <= 0 || x %% 1 != 0) {
    stop(paste(var_name, "must be a positive integer."))
  }
}

check_boolean <- function(x, var_name = "value") {
  if (!is.logical(x) || length(x) != 1) {
    stop(paste(var_name, "must be a boolean (TRUE or FALSE)."))
  }
}

#' Generate K-fold Indices for Cross-Validation
#'
#' This function generates indices for K-fold cross-validation based on the total sample size `N` and the number of folds `Kfold`.
#' If `reverse = TRUE`, the remainder indices will be assigned in reverse order.
#'
#' @param N Integer. Total sample size (number of observations).
#' @param Kfold Integer. The number of folds to split the data into.
#' @param reverse Logical. Whether to reverse the remainder indices when `N` is not divisible by `Kfold`. Defaults to `FALSE`.
#'
#' @author Thomas Debray
#' @return A vector of length `N` containing the fold assignments (from 1 to `Kfold`).
generate_kfold_indices <- function(N, Kfold, reverse = FALSE) {
  base_index <- rep(seq(Kfold), floor(N / Kfold))  # Base repeated sequence
  remainder <- N %% Kfold                        # Calculate the remainder

  if (remainder > 0) {
    additional_index <- if (reverse) {
      rev(seq(remainder))                         # Reverse the remainder for index0
    } else {
      seq(remainder)                              # Normal order for index1
    }
    base_index <- c(base_index, additional_index)
  }

  return(base_index)
}
