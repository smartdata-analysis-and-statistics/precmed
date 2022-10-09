#' Estimate restricted mean survival time (RMST) based on Cox regression model
#'
#'
#' @param y Observed survival or censoring time; vector of size \code{n}.
#' @param d The event indicator, normally \code{1 = event, 0 = censored}; vector of size \code{n}.
#' @param x.cate Matrix of \code{p.cate} baseline covariates specified in the outcome model; dimension \code{n} by \code{p.cate}.
#' @param xnew Matrix of \code{p.cate} baseline covariates for which we want an estimate of the RMST; dimension \code{m} (observations in the new data set) by \code{p.cate}
#' @param tau0 The truncation time for defining restricted mean time lost.
#'
#' @return The estimated RMST for new subjects with covariates \code{xnew}; vector of size \code{m}.
#'
cox.rmst <- function(y, d, x.cate, xnew, tau0) {
  xnew <- as.matrix(xnew)
  x <- as.matrix(x.cate)

  vx <- apply(x, 2, var)
  xnew <- xnew[, vx > 0, drop = FALSE]
  x <- x[, vx > 0, drop = FALSE]

  mux <- colMeans(x)
  xnew <- t(t(xnew) - mux)
  x <- t(t(x) - mux)

  fit <- coxph(Surv(y, d) ~ x)
  beta <- fit$coef
  detail <- coxph.detail(fit)
  time <- detail$time
  surv0 <- exp(-cumsum(detail$hazard))

  n <- length(xnew[, 1])
  surv0 <- surv0[time <= tau0]
  time <- time[time <= tau0]
  hr <- exp(as.vector(xnew[, is.na(beta) == FALSE, drop = FALSE] %*% beta[is.na(beta) == FALSE, drop = FALSE]))
  surv <- crossprod(t(rep(1, n)), surv0)^hr

  gap <- c(time, tau0) - c(0, time)
  coxest <- colSums(t(cbind(1, surv)) * gap)

  return(coxest)
}

#' Catch errors and warnings when estimating the ATEs in the nested subgroup
#'
#' Storing the errors and warnings that occurred when estimating the ATEs in the nested subgroups.
#' If there are no errors and no warnings, the estimated log.rmtl.ratio and log.hazard.ratio are provided.
#' If there are warnings but no errors, the estimated log.rmtl.ratio and log.hazard.ratio are provided with a warning attribute set.
#' If there are errors, the NA values are returned for log.rmtl.ratio and log.hazard.ratio. A error attribute set is also provided.
#'
#' @param fun The drsurv function...
#' @return A list containing

survCatch <- function(fun) {
  err.msg <- err.call <- NULL
  warn.msg <- warn.call <- NULL
  fit <- withCallingHandlers(
    tryCatch(expr = fun,
             error = function(e) {
               err.msg <<- conditionMessage(e)
               err.call <<- conditionCall(e)
               list(log.rmtl.ratio = NA, log.hazard.ratio = NA)
             }),
    warning = function(w) {
      warn.msg <<- append(warn.msg, conditionMessage(w))
      warn.call <<- append(warn.call, conditionCall(w))
      invokeRestart("muffleWarning")
    })

  if (!is.null(err.msg)) {
    err.call <- toString(as.expression(err.call))
    errors <- paste0("Error in ", err.call, " : ", err.msg)
  } else {
    errors <- NULL
  }

  if (!is.null(warn.msg)) {
    warn.call <- toString(as.expression(warn.call))
    warnings <- paste0("Warning in ", warn.call, " : ", warn.msg)
  } else {
    warnings <- NULL
  }

  out <- c()
  out$log.rmtl.ratio <- fit$log.rmtl.ratio
  out$log.hazard.ratio <- fit$log.hazard.ratio
  out$warnings <- warnings
  out$errors <- errors
  return(out)
}

#' Data preprocessing
#' Apply at the beginning of \code{catefitcount()}, \code{catecvcount()}, \code{catefitsurv()}, and \code{catecvsurv()}, after \code{arg.checks()}
#'
#' @param fun A function for which argument check is needed; "catefit" for \code{catefitcount()} and \code{catefitsurv()},
#' "crossv" for \code{catecvcount()} and \code{catecvsurv()},
#' and "drinf" for \code{drcount.inference()} and \code{drsurv.inference()}. No default.
#' @param cate.model A formula describing the outcome model to be fitted.
#' The outcome must appear on the left-hand side.
#' @param ps.model A formula describing the propensity score model to be fitted.
#' The treatment must appear on the left-hand side. The treatment must be a numeric vector
#' coded as 0/1. If data are from a RCT, specify \code{ps.model} as an intercept-only model.
#' @param ipcw.model A formula describing inverse probability of censoring weighting(IPCW) model to be fitted.
#' If covariates are the same as outcome model, set \code{ipcw.model = NULL}.
#' Otherwise, the left-hand side must be empty and the right-hand side is a covariates model.
#' @param tau0 The truncation time for defining restricted mean time lost. Default is \code{NULL},
#' which corresponds to setting the truncation time as the maximum survival time in the data
#' @param data A data frame containing the variables in the outcome, propensity score, and IPCW models;
#' a data frame with \code{n} rows (1 row per observation).
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
#' @param initial.predictor.method A character vector for the method used to get initial
#' outcome predictions conditional on the covariates. Only applies when \code{score.method}
#' includes \code{'twoReg'} or \code{'contrastReg'}. Allowed values include one of
#' \code{'randomForest'} (survival outcomes only), \code{'boosting'}, \code{'logistic'}
#' (survival outcomes only, fast), \code{'poisson'} (count outcomes only, fast), and
#' \code{'gam'} (count outcomes only). Default is \code{NULL}, which assigns \code{'boosting'}
#' for count outcomes and \code{'randomForest'} for survival outcomes.
#' @param response The type of response variables; \code{count} (default) or \code{survival}.
#'
#' @return A list of elements:
#'            - y: outcome; vector of length \code{n} (observations)
#'            - d : the event indicator; vector of length \code{n}; only if \code{respone = "survival"}
#'            - trt: binary treatment; vector of length \code{n}
#'            - x.ps: matrix of \code{p.ps} baseline covariates specified in the propensity score model (plus intercept); dimension \code{n} by \code{p.ps + 1}
#'            - x.cate: matrix of \code{p.cate} baseline covariates specified in the outcome model; dimension \code{n} by \code{p.cate}
#'            - x.ipcw: matrix of \code{p.ipw} baseline covarites specified in inverse probability of censoring weighting model; dimension \code{n} by \code{p.ipw}
#'            - time: offset; vector of length \code{n}; only if \code{response = "count"}
#'            - if \code{fun = "catefit"}:
#'                - prop: formatted \code{prop.cutoff}
#'                - prop.no1: formatted \code{prop.cutoff} with 1 removed if applicable; otherwise prop.no1 is the same as prop
#'            - if \code{fun = "crossv"}
#'                - prop.onlyhigh: formatted \code{prop.cutoff} with 0 removed if applicable
#'                - prop.bi; formatted \code{prop.cutoff} with 0 and 1 removed if applicable
#'                - prop.multi: formatted \code{prop.multi}, starting with 0 and ending with 1
#'


data.preproc.surv <- function(fun, cate.model, ps.model, ipcw.model = NULL, tau0 = NULL, data, prop.cutoff = NULL, prop.multi = NULL,
                              ps.method, initial.predictor.method = NULL, response = "count") {

  if (!(fun %in% c("catefit", "crossv", "drinf"))) {
    stop("Invalid source function!")
  }
  if (!(response %in% c("count", "survival", "continuous"))) {
    stop("Response not supported!")
  }

  ## cate.model and ps.model as formulas
  cate.model <- as.formula(cate.model)
  ps.model <- as.formula(ps.model)

  ## Extraction step
  ### If count data: extract y, time, trt, x.cate, x.ps in matrix form from cate.model and ps.model
  ### If survival data: extract y, d, trt, x.cate, x.ps, x.ipcw in matrix form from cate.model, ps.model and ipcw.model
  cate.mat <- model.frame(cate.model, data, na.action = 'na.pass')
  ps.mat <- model.frame(ps.model, data, na.action = 'na.pass')
  resp <- model.response(cate.mat)
  if (is.null(resp)) {
    stop("Outcome must be supplied on the left-hand side of the cate.model formula.")
  } else {
    if (response == "count") {
      y <- resp
      time <- model.offset(cate.mat)
      if (is.null(time) == TRUE) { # Eventually some outcomes will not need an offset
        time <- rep(0, nrow(data))
        warning("No offset supplied. Offset set to 0.")
      }
    } else if (response == "survival") {
      y <- resp[, 1]
      d <- resp[, 2]
    }
  }

  # Assign default to tau0 if NULL
  if (response == "survival" & is.null(tau0) == TRUE) {
    tau0 <- max(y[d == 1])
    warning("No value supplied for tau0. Default sets tau0 to the maximum survival time.")
  }

  # trt
  trt <- model.response(ps.mat)
  if (is.null(trt)) stop("Treatment must be supplied on the left-hand side of the ps.model formula.")

  # Covariate matrices
  x.cate <- model.matrix(cate.model, cate.mat)[, -1, drop = FALSE]
  x.ps <- model.matrix(ps.model, ps.mat)
  if (ncol(x.ps) == 1 & ps.method == "lasso") stop("LASSO penalization irrelevant when ps.model specified as a function of an intercept only. Consider setting ps.method='glm'.")

  # Covariate matrices for IPCW
  if (response == "survival") {
    if (is.null(ipcw.model)) {
      x.ipcw <- cbind(x.cate, trt)
    } else {
      ipcw.model <- as.formula(ipcw.model)
      ipcw.mat <- model.frame(ipcw.model, data, na.action = 'na.pass')
      x.ipcw <- model.matrix(ipcw.model, ipcw.mat)[, -1, drop = FALSE]
    }
  }

  ## Check negative count or survival
  if (any(y < 0)) stop(paste0("Negative y values not allowed for ", response, " outcomes."))

  ## Check missing data
  if (response == "count") {
    if (any(is.na(y)) | any(is.na(time)) | any(is.na(trt)) | any(is.na(x.cate)) | any(is.na(x.ps))) stop("Missing data not allowed in cate.model, ps.model or offset.")
  } else if (response == "survival") {
    if (any(is.na(y)) | any(is.na(trt)) | any(is.na(x.cate)) | any(is.na(x.ps)) | any(is.na(x.ipcw))) stop("Missing data not allowed in cate.model or ps.model")
  }

  ## Check treatment binary numeric and coded as 0/1
  cat <- warn.dr <- NULL
  if (as.numeric(length(unique(trt))) > 2) {
    stop("trt must be binary.")
  } else if (!all(unique(trt) %in% c(0, 1))) {
    cat <- cat <- sort(unique(trt))
    trt <- ifelse(trt == cat[1], 0, 1)
    warn.dr <- warning(paste0("Variable trt was recoded to 0/1 with ", cat[1], "->0 and ", cat[2], "->1.\n"))
  } else if (is.factor(trt) == TRUE) {
    trt <- as.numeric(trt == 1)
  }

  ## Check event binary numeric and coded as 0/1 for survival outcomes
  if (response == "survival") {
    if (sum(is.na(d)) == length(d)) {
      stop("The event indicator must be supplied as 0/1.")
    } else if (any(is.na(d))) {
      stop("Missing data in event indicator not allowed.")
    } else if (as.numeric(length(unique(d))) > 2) {
      stop("The event indicator must be binary.")
    }
  }

  ## Assign default to initial.predictor.method if NULL
  if (is.null(initial.predictor.method) == TRUE & fun %in% c("catefit", "crossv")) {
    if (response == "count") {
      initial.predictor.method <- "boosting"
    } else if (response == "survival") {
      initial.predictor.method <- "randomForest"
    }
  }

  ## Output data
  if (response == "count") {
    out.data <- list(y = y, time = time, trt = trt, x.cate = x.cate, x.ps = x.ps, cat.trt = cat, tau0 = tau0, initial.predictor.method = initial.predictor.method)
  } else if (response == "survival") {
    out.data <- list(y = y, d = d, trt = trt, x.cate = x.cate, x.ps = x.ps, x.ipcw = x.ipcw, cat.trt = cat, tau0 = tau0, initial.predictor.method = initial.predictor.method)
  }

  if (fun == "catefit") {
    ## Check values of prop
    prop <- sort(prop.cutoff) # sort the proportions from small to large
    if (prop[1] == 0) {
      prop <- prop[-1] # if first element is 0, remove it because this means we leave out 0% of individuals
      warning("The first element of prop.cutoff cannot be 0 and has been removed.")
    }
    out.data$prop.no1 <- prop
    if (prop[length(prop)] == 1) {
      out.data$prop.no1 <- prop[-length(prop)]
    }
    out.data$prop <- prop
    return(out.data)

  }
  if (fun == "crossv") {
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

    out.data$prop.onlyhigh <- prop.onlyhigh
    out.data$prop.bi <- prop.bi
    out.data$prop.multi <- prop.multi

    return(out.data)

  }

  # fun == drinf
  out.data$warning <- warn.dr
  return(out.data)
}

#' Split the given time-to-event dataset into balanced training and validation sets (within a pre-specified tolerance)
#' Balanced means 1) The ratio of treated and controls is maintained in the training and validation sets
#'                2) The covariate distributions are balanced between the training and validation sets
#'
#' @param y Observed survival or censoring time; vector of size \code{n}.
#' @param d The event indicator, normally \code{1 = event, 0 = censored}; vector of size \code{n}.
#' @param trt Treatment received; vector of size \code{n} with treatment coded as 0/1.
#' @param x.cate Matrix of \code{p.cate} baseline covariates specified in the outcome model; dimension \code{n} by \code{p.cate}.
#' @param x.ps Matrix of \code{p.ps} baseline covariates specified in the propensity score model; dimension \code{n} by \code{p.ps}.
#' @param x.ipcw Matrix of \code{p.ipw} baseline covariate specified in inverse probability of censoring weighting; dimension \code{n} by \code{p.ipw}.
#' @param yf Follow-up time, interpreted as the potential censoring time; vector of size \code{n} if the potential censoring time is known.
#' If unknown, set \code{yf == NULL} and \code{yf} will be taken as \code{y} in the function.
#' @param train.prop A numerical value (in (0, 1)) indicating the proportion of total data used
#' for training. Default is \code{3/4}.
#' @param error.max A numerical value > 0 indicating the tolerance (maximum value of error)
#' for the largest standardized absolute difference in the covariate distributions or in the
#' doubly robust estimated rate ratios between the training and validation sets. This is used
#' to define a balanced training-validation splitting. Default is \code{0.1}.
#' @param max.iter A positive integer value indicating the maximum number of iterations when
#' searching for a balanced training-validation split. Default is \code{5,000}.
#'
#' @return A list of 14 objects, 7training and 7 validation of y, trt, x.cate, x.ps, x.ipcw, time, yf:
#'             y.train          - observed survival or censoring time in the training set; vector of size \code{m} (observations in the training set)
#'             d.train          - event indicator in the training set; vector of size \code{m} coded as 0/1
#'             trt.train        - treatment received in the training set; vector of size \code{m} coded as 0/1
#'             x.cate.train     - baseline covariates for the outcome model in the training set; matrix of dimension \code{m} by \code{p.cate}
#'             x.ps.train       - baseline covariates (plus intercept) for the propensity score model in the training set; matrix of dimension \code{m} by \code{p.ps + 1}
#'             x.ipcw.train      - baseline covariates for inverse probability of censoring in the training set; matrix of dimension \code{m} by \code{p.ipw}
#'             yf.train         - follow-up time in the training set; if known, vector of size \code{m}; if unknown, \code{yf == NULL}
#'             y.valid          - observed survival or censoring time in the validation set; vector of size \code{n-m}
#'             d.valid          - event indicator in the validation set; vector of size \code{n-m} coded as 0/1
#'             trt.valid        - treatment received in the validation set; vector of size \code{n-m} coded as 0/1
#'             x.cate.valid     - baseline covariates for the outcome model in the validation set; matrix of dimension \code{n-m} by \code{p.cate}
#'             x.ps.valid       - baseline covariates (plus intercept) for the propensity score model in the validation set; matrix of dimension \code{n-m} by \code{p.ps + 1}
#'             x.ipcw.valid      - baseline covariates for inverse probability of censoring in the validation set; matrix of dimension \code{n-m} by \code{p.ipw}
#'             yf.valid         - follow-up time in the training set; if known, vector of size \code{n-m}; if unknown, \code{yf == NULL}
#'

balancesurv.split <- function(y, d, trt, x.cate, x.ps, x.ipcw, yf = NULL, train.prop = 3/4, error.max = 0.1, max.iter = 5000) {

  x <- cbind(x.cate, x.ps[, -1], x.ipcw)
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
  error <- errorbest <- Inf
  bestid.valid <- rep(NA, n - m)
  iter <- 0

  while ((error > error.max | is.na(error) == TRUE) & iter <= max.iter) {
    id.valid <- c(sample(id1, n1 - m1, replace = FALSE), sample(id0, n0 - m0, replace = FALSE))

    y.valid <- y[id.valid]
    d.valid <- d[id.valid]
    trt.valid <- trt[id.valid]
    x.valid  <- x[id.valid, , drop = FALSE]
    x.cate.valid <- x.cate[id.valid, , drop = FALSE]
    x.ps.valid <- x.ps[id.valid, , drop = FALSE]
    x.ipcw.valid <- x.ipcw[id.valid, , drop = FALSE]
    yf.valid <- yf[id.valid]

    y.train <- y[-id.valid]
    d.train <- d[-id.valid]
    trt.train <- trt[-id.valid]
    x.train <- x[-id.valid, , drop = FALSE]
    x.cate.train <- x.cate[-id.valid, , drop = FALSE]
    x.ps.train <- x.ps[-id.valid, , drop = FALSE]
    x.ipcw.train <- x.ipcw[-id.valid, , drop = FALSE]
    yf.train <- yf[-id.valid]

    diffx <- colMeans(x.train) - colMeans(x.valid)
    errorx <- max(abs(diffx / sdx))

    est.valid <- coxph(Surv(y.valid, d.valid) ~ trt.valid + x.valid)$coef[1]
    est.train <- coxph(Surv(y.train, d.train) ~ trt.train + x.train)$coef[1]

    error <- max(2 * abs(est.valid - est.train) / (abs(est.valid) + abs(est.train)), errorx)
    if (is.na(error) == FALSE & error < errorbest) {
      bestid.valid <- id.valid
      errorbest <- error
    }
    iter <- iter + 1
  }

  if (iter == max.iter + 1) {
    y.valid <- y[bestid.valid]
    d.valid <- d[bestid.valid]
    trt.valid <- trt[bestid.valid]
    x.cate.valid <- x.cate[bestid.valid, , drop = FALSE]
    x.ps.valid <- x.ps[bestid.valid, , drop = FALSE]
    x.ipcw.valid <- x.ipcw[bestid.valid, , drop = FALSE]

    y.train <- y[-bestid.valid]
    d.train <- d[-bestid.valid]
    trt.train <- trt[-bestid.valid]
    x.cate.train <- x.cate[-bestid.valid, , drop = FALSE]
    x.ps.train <- x.ps[-bestid.valid, , drop = FALSE]
    x.ipcw.train <- x.ipcw[-bestid.valid, , drop = FALSE]
    warning(paste("Maximum iteration reached and the SMD between training and validation set is still greater than error.max (error=", round(errorbest, 4), "). Consider increasing max.iter, decreasing error.max, or increasing sample size.", sep = ""))
  }

  return(list(y.train = y.train, d.train = d.train, trt.train = trt.train, x.cate.train = x.cate.train,
              x.ps.train = x.ps.train, x.ipcw.train = x.ipcw.train, yf.train = yf.train,
              y.valid = y.valid, d.valid = d.valid, trt.valid = trt.valid, x.cate.valid = x.cate.valid,
              x.ps.valid = x.ps.valid, x.ipcw.valid = x.ipcw.valid, yf.valid = yf.valid))
}

