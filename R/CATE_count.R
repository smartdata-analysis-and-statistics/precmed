#' Estimation of the conditional average treatment effect (CATE) score for count data
#'
#' Provides singly robust and doubly robust estimation of CATE score with up to 5 scoring methods
#' among the following: Poisson regression, boosting, two regressions, contrast regression, and
#' negative binomial.
#'
#' @param data A data frame containing the variables in the outcome and propensity score model; a data frame with \code{n} rows
#' (1 row per observation).
#' @param score.method A vector of one or multiple methods to estimate the CATE score.
#' Allowed values are: \code{'boosting'}, \code{'poisson'}, \code{'twoReg'}, \code{'contrastReg'},
#' and \code{'negBin'}.
#' @param cate.model A formula describing the outcome model to be fitted.
#' The outcome must appear on the left-hand side.
#' @param ps.model A formula describing the propensity score model to be fitted.
#' The treatment must appear on the left-hand side. The treatment must be a numeric vector
#' coded as 0/1. If data are from a randomized trial, specify \code{ps.model} as an intercept-only model.
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
#' fit <- catefitcount(data = countExample,
#'                     score.method = "poisson",
#'                     cate.model = y ~ age + female + previous_treatment +
#'                                  previous_cost + previous_number_relapses +
#'                                  offset(log(years)),
#'                     ps.model = trt ~ age + previous_treatment,
#'                     higher.y = FALSE,
#'                     seed = 999, verbose = 1)
#' @export

catefitcount <- function(data,
                         score.method,
                         cate.model,
                         ps.model,
                         ps.method = "glm",
                         initial.predictor.method = "boosting",
                         minPS = 0.01,
                         maxPS = 0.99,
                         higher.y = TRUE,
                         prop.cutoff = seq(0.5, 1, length = 6),
                         xvar.smooth = NULL,
                         tree.depth = 2,
                         n.trees.boosting = 200,
                         B = 3,
                         Kfold = 5,
                         error.maxNR = 1e-3,
                         max.iterNR = 150,
                         tune = c(0.5, 2),
                         seed = NULL,
                         plot.gbmperf = FALSE,
                         verbose = 0,
                         ...) {


  # Set seed once for reproducibility
  set.seed(seed)

  if (verbose >= 1) t.start <- Sys.time()

  #### CHECK ARGUMENTS ####
  arg.checks(
    fun = "catefit", response = "count", data = data, higher.y = higher.y, score.method = score.method, prop.cutoff = prop.cutoff,
    ps.method = ps.method, minPS = minPS, maxPS = maxPS,
    initial.predictor.method = initial.predictor.method,
    tree.depth = tree.depth, n.trees.boosting = n.trees.boosting, B = B, Kfold = Kfold, plot.gbmperf = plot.gbmperf,
    error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune
  )

  #### PRE-PROCESSING ####
  out <- data.preproc(fun = "catefit", cate.model = cate.model, ps.model = ps.model,
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
    warning(paste0("The best boosting iteration was iteration number",
                   n.trees.boosting, " out of ", n.trees.boosting,
                   ". Consider increasing the maximum number of trees and
                  turning on boosting performance plot (plot.gbmperf = TRUE)."))
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

  result$response <- "count"
  result$score.method <- score.method

  class(result) <- "catefit"
  return(result)
}

#' Doubly robust estimators of the coefficients in the two regression
#'
#' @param y Observed outcome; vector of size \code{n}
#' @param x.cate Matrix of \code{p} baseline covariates; dimension \code{n} by \code{p}
#' @param time Log-transformed person-years of follow-up; vector of size \code{n}
#' @param trt Treatment received; vector of size \code{n} units with treatment coded as 0/1
#' @param ps Estimated propensity scores for all observations; vector of size \code{n}
#' @param f.predictor Initial prediction of the outcome (expected number of relapses for one unit of exposure time) conditioned
#' on the covariates \code{x} for one treatment group \code{r}; \code{mu_r(x)}, step 1 in the two regression; vector of size \code{n}
#'
#' @return Doubly robust estimators of the regression coefficients \code{beta_r} in the doubly robust estimating equation
#' where \code{r = 0, 1} is treatment received; vector of size \code{p} + 1 (intercept included)
onearmglmcount.dr <- function(y, x.cate, time, trt, ps, f.predictor) {
  f.predictor <- as.vector(f.predictor)
  y <- y * exp(-time)
  x <- as.matrix(cbind(1, log(f.predictor), x.cate))

  withCallingHandlers({
    fit <- glm(y ~ log(f.predictor) + x.cate, family = "poisson", weights = trt / ps)
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
#' @param y Observed outcome; vector of size \code{n}
#' @param x.cate Matrix of \code{p.cate} baseline covariates; dimension \code{n} by \code{p.cate}
#' @param time Log-transformed person-years of follow-up; vector of size \code{n}
#' @param trt Treatment received; vector of size \code{n} units with treatment coded as 0/1
#' @param ps Estimated propensity scores for all observations; vector of size \code{n}
#' @param f1.predictor Initial predictions of the outcome (expected number of relapses for one unit of exposure time)
#' conditioned on the covariates \code{x} for treatment group trt = 1; \code{mu_1(x)}, step 1 in the two regression; vector of size \code{n}
#' @param f0.predictor  Initial predictions of the outcome (expected number of relapses for one unit of exposure time)
#' conditioned on the covariates \code{x} for treatment group trt = 0; \code{mu_0(x)}, step 1 in the two regression; vector of size \code{n}
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
#'
#' @return coef: Doubly robust estimators of the regression coefficients \code{delta_0}; vector of size \code{p} + 1 (intercept included)
#'         vcov: Variance-covariance matrix of the estimated coefficient \code{delta_0}; matrix of size \code{p} + 1 by \code{p} + 1
#'         converge: Indicator that the Newton Raphson algorithm converged for \code{delta_0}; boolean
twoarmglmcount.dr <- function(y, x.cate, time, trt,
                              ps, f1.predictor, f0.predictor,
                              error.maxNR = 1e-3, max.iterNR = 150, tune = c(0.5, 2.0)) {

  y <- y * exp(-time)
  x.aug <- cbind(1, x.cate)
  p.aug <- length(x.aug[1, ])

  # Estimate coefficients with NR algorithm
  beta <- rep(0, p.aug)
  mae <- Inf
  iter <- 0
  epsilon.min <- Inf
  while (mae > error.maxNR && iter <= max.iterNR && epsilon.min > 0) {
    eta <- as.numeric(exp(x.aug %*% beta))
    error <- (trt * (y - eta * f0.predictor / 2 - f1.predictor / 2) * (1 - ps) - (1 - trt) * (eta * y - eta * f0.predictor / 2 - f1.predictor / 2) * ps) / (eta * ps + (1 - ps))
    score <- colSums(x.aug * error)
    slopewt <- (y + f0.predictor * (trt / ps - 1) / 2 + f1.predictor * ((1 - trt) / (1 - ps) - 1) / 2) * eta * ps * (1 - ps) / (eta * ps + (1 - ps))^2
    slope <- t(x.aug * slopewt) %*% x.aug
    epsilon.min <- eigen(slope, only.values = TRUE)$value[p.aug]
    if (iter == 0) epsilon0 <- epsilon.min + epsilon.min * (epsilon.min < 0) # fixed to all iterations
    beta <- beta + solve(slope + diag(tune[2] * abs(epsilon0), p.aug, p.aug)) %*% score * tune[1] # adding the diagonal matrix to slove potential singualrity issues
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
      score <- colSums(x.aug * error)
      return(sum(score^2))
    }

    initial.value <- lossf(rep(0, p.aug))
    fit <- optim(rep(0, p.aug), fn = lossf) #, control=list(trace=T))
    beta <- fit$par
    converge2 <- 1 * (fit$value < initial.value / 100)
  }

  beta <- as.vector(beta)
  eta <- as.numeric(exp(x.aug %*% beta))
  error <- (trt * (y - eta * f0.predictor / 2 - f1.predictor / 2) * (1 - ps) - (1 - trt) * (eta * y - eta * f0.predictor / 2 - f1.predictor / 2) * ps) / (eta * ps + (1 - ps))
  score <- colSums(x.aug * error)
  slopewt <- (y + f0.predictor * (trt / ps - 1) / 2 + f1.predictor * ((1 - trt) / (1 - ps) - 1) / 2) * eta * ps * (1 - ps) / (eta * ps + (1 - ps))^2
  slope <- t(x.aug * slopewt) %*% x.aug
  sigma <- solve(slope) %*% (t(x.aug * error^2) %*% x.aug) %*% solve(slope)

  return(list(coef = beta, vcov = sigma, converge = (converge1 + converge2 > 0)))
}


#' Estimate the CATE model using specified scoring methods
#'
#' Coefficients of the CATE estimated with boosting, naive Poisson, two regression, contrast regression, negative binomial
#'
#' @param y Observed outcome; vector of size \code{n} (observations)
#' @param x.cate Matrix of \code{p.cate} baseline covariates; dimension \code{n} by \code{p.cate} (covariates in the outcome model)
#' @param x.ps Matrix of \code{p.ps} baseline covariates (plus a leading column of 1 for the intercept);
#' dimension \code{n} by \code{p.ps + 1} (covariates in the propensity score model plus intercept)
#' @param time Log-transformed person-years of follow-up; vector of size \code{n}
#' @param trt Treatment received; vector of size \code{n} units with treatment coded as 0/1
#' @param score.method A vector of one or multiple methods to estimate the CATE score.
#' Allowed values are: \code{'boosting'}, \code{'poisson'}, \code{'twoReg'}, \code{'contrastReg'},
#' \code{'negBin'}. Default specifies all 5 methods.
#' @param ps.method A character value for the method to estimate the propensity score.
#' Allowed values include one of:
#' \code{'glm'} for logistic regression with main effects only (default), or
#' \code{'lasso'} for a logistic regression with main effects and LASSO penalization on
#' two-way interactions (added to the model if interactions are not specified in \code{ps.model}).
#' Relevant only when \code{ps.model} has more than one variable.
#' @param minPS A numerical value (in [0, 1]) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A number above which estimated propensity scores should be trimmed; scalar
#' @param initial.predictor.method A character vector for the method used to get initial
#' outcome predictions conditional on the covariates in \code{cate.model}
#' in \code{score.method = 'twoReg'} and \code{'contrastReg'}. Allowed values include
#' one of \code{'poisson'} (fastest), \code{'boosting'} (default) and \code{'gam'}.
#' @param xvar.smooth A vector of characters indicating the name of the variables used as
#' the smooth terms if \code{initial.predictor.method = 'gam'}. The variables must be selected
#' from the variables listed in \code{cate.model}.
#' Default is \code{NULL}, which uses all variables in \code{cate.model}.
#' @param tree.depth A positive integer specifying the depth of individual trees in boosting
#' (usually 2-3). Used only if \code{score.method = 'boosting'} or
#' if \code{score.method = 'twoReg'} or \code{'contrastReg'} and
#' \code{initial.predictor.method = 'boosting'}. Default is \code{2}.
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
#' @param max.iterNR A positive integer indicating the maximum number of iterations in the
#' Newton Raphson algorithm. Used only if \code{score.method = 'contrastReg'}.
#' Default is \code{150}.
#' @param tune A vector of 2 numerical values > 0 specifying tuning parameters for the
#' Newton Raphson algorithm. \code{tune[1]} is the step size, \code{tune[2]} specifies a
#' quantity to be added to diagonal of the slope matrix to prevent singularity.
#' Used only if \code{score.method = 'contrastReg'}. Default is \code{c(0.5, 2)}.
#' @param ... Additional arguments for \code{gbm()}
#'
#' @return Depending on what score.method is, the outputs is a combination of the following:
#'           result.boosting: Results of boosting fit and best iteration, for trt = 0 and trt = 1 separately
#'           result.poisson: Naive Poisson estimator (beta1 - beta0); vector of length \code{p.cate} + 1
#'           result.twoReg: Two regression estimator (beta1 - beta0); vector of length \code{p.cate} + 1
#'           result.contrastReg: A list of the contrast regression results with 3 elements:
#'               $delta.contrastReg: Contrast regression DR estimator; vector of length \code{p.cate} + 1
#'               $sigma.contrastReg: Variance covariance matrix for delta.contrastReg; matrix of size \code{p.cate} + 1 by \code{p.cate} + 1
#'               $converge.contrastReg: Indicator that the Newton Raphson algorithm converged for \code{delta_0}; boolean
#'           result.negBin: Negative binomial estimator (beta1 - beta0); vector of length \code{p.cate} + 1
#'           best.iter: Largest best iterations for boosting (if used)
#'           fgam: Formula applied in GAM (if used)
intxcount <- function(y, trt, x.cate, x.ps, time,
                      score.method = c("boosting", "poisson", "twoReg", "contrastReg", "negBin"),
                      ps.method = "glm", minPS = 0.01, maxPS = 0.99,
                      initial.predictor.method = "boosting",
                      xvar.smooth = NULL,
                      tree.depth = 2, n.trees.boosting = 200, B = 3, Kfold = 6, plot.gbmperf = TRUE,
                      error.maxNR = 1e-3, max.iterNR = 150, tune = c(0.5, 2.0), ...) {

  result <- vector("list", length(score.method) + 1)
  names(result) <- c(paste0("result.", score.method), "best.iter")

  N1 <- sum(trt)
  N0 <- sum(1 - trt)
  N <- N1 + N0
  p.aug <- ncol(x.cate) + 1

  datatot <- data.frame(y, x.cate, time)
  colnames(datatot) <- c("y", colnames(x.cate), "time")

  ######### cross-fitting  ---------------------------------------------------------------

  index1 <- rep(1:Kfold, floor(N1 / Kfold))
  if (N1 > Kfold * floor(N1 / Kfold)) index1 <- c(index1, 1:(N1 - Kfold * floor(N1 / Kfold)))

  index0 <- rep(1:Kfold, floor(N0 / Kfold))
  if (N0 > Kfold * floor(N0 / Kfold)) index0 <- c(index0, Kfold + 1 - 1:(N0 - Kfold * floor(N0 / Kfold)))

  delta.twoReg.mat <- delta.contrastReg.mat <- matrix(NA, B, p.aug)
  sigma.contrastReg.mat <- matrix(0, p.aug, p.aug)
  converge <- rep(NA, B)
  best.iter <- 0
  fgam <- NULL

  if (any(c("twoReg", "contrastReg") %in% score.method)) {
    for (bb in seq(B)) {
      index1cv <- sample(x = index1, size = N1, replace = FALSE)
      index0cv <- sample(x = index0, size = N0, replace = FALSE)
      index <- rep(NA, N)
      index[trt == 1] <- index1cv
      index[trt == 0] <- index0cv

      f1.predictcv <- f0.predictcv <- pscv <- rep(NA, N)
      for (k in 1:Kfold) {
        datatot_train <- datatot[index != k, ]
        x_ps_train <- x.ps[index != k, , drop = FALSE]
        trt_train <- trt[index != k]

        datatot_valid <- datatot[index == k, ]
        x_ps_valid <- x.ps[index == k, , drop = FALSE]

        data1 <- datatot_train[trt_train == 1, ]
        data0 <- datatot_train[trt_train == 0, ]

        if (initial.predictor.method == "boosting") {
          # if model has a single predictor, GBM must have cv.folds = 0 https://github.com/zoonproject/zoon/issues/130
          cate.cvfold <- ifelse(ncol(x.cate) == 1, 0, 5)
          fit1.boosting <- gbm(y ~ . - time + offset(time), data = data1, distribution = "poisson",
                                 interaction.depth = tree.depth, n.trees = n.trees.boosting, cv.folds = cate.cvfold, ...)
          best1.iter <- max(10, gbm.perf(fit1.boosting, method = "cv", plot.it = plot.gbmperf))


          fit1.boosting <- gbm(y ~ . - time + offset(time), data = data1, distribution = "poisson",
                               interaction.depth = tree.depth, n.trees = n.trees.boosting, cv.folds = cate.cvfold, ...)
          best1.iter <- max(10, gbm.perf(fit1.boosting, method = "cv", plot.it = plot.gbmperf))
          withCallingHandlers({
            f1.predictcv[index == k] <- predict(object = fit1.boosting, newdata = datatot_valid, n.trees = best1.iter, type = "response")
          },
          warning = function(w) {
            if (grepl("does not add the offset", conditionMessage(w)))
              invokeRestart("muffleWarning")  # suppress warning: "predict.gbm does not add the offset to the predicted values."
          })

          fit0.boosting <- gbm(y ~ . - time + offset(time), data = data0, distribution = "poisson",
                               interaction.depth = tree.depth, n.trees = n.trees.boosting, cv.folds = cate.cvfold, ...)
          best0.iter <- max(10, gbm.perf(fit0.boosting, method = "cv", plot.it = plot.gbmperf))
          withCallingHandlers({
            f0.predictcv[index == k] <- predict(object = fit0.boosting, newdata = datatot_valid, n.trees = best0.iter, type = "response")
          },
          warning = function(w) {
            if (grepl("does not add the offset", conditionMessage(w)))
              invokeRestart("muffleWarning")  # suppress warning: "predict.gbm does not add the offset to the predicted values."
          })

          best.iter <- max(best.iter, best1.iter, best0.iter)


        } else if (initial.predictor.method == "poisson") {

          fit1.pois <- glm(y ~ . - time + offset(time), data = data1, family = "poisson")
          f1.predictcv[index == k] <- predict(object = fit1.pois, newdata = datatot_valid, type = "response")

          fit0.pois <- glm(y ~ . - time + offset(time), data = data0, family = "poisson")
          f0.predictcv[index == k] <- predict(object = fit0.pois, newdata = datatot_valid, type = "response")

        } else if (initial.predictor.method == "gam") {

          xvars <- colnames(x.cate)
          if (is.null(xvar.smooth)) {
            fgam <- paste0("y ~ ", paste0("s(", xvars, ")", collapse = "+"))
          } else {
            xvar.smooth2 <- xvars[str_detect(xvars, paste(paste0(xvar.smooth, "$"), collapse = "|"))] # conver to the preprocessed names
            xvar.linear <- setdiff(xvars, xvar.smooth2) # the remaining xvars in x.cate but not in xvar.smooth are linear predictors
            fgam <- paste0("y ~ ", paste0(xvar.linear, collapse = "+"), "+", paste0("s(", xvar.smooth2, ")", collapse = "+"))
          }
          fit1.gam <- mgcv::gam(as.formula(fgam), offset = time, data = data1, family = "poisson")
          f1.predictcv[index == k] <- predict(object = fit1.gam, newdata = datatot_valid, type = "response")

          fit0.gam <- mgcv::gam(as.formula(fgam), offset = time, data = data0, family = "poisson")
          f0.predictcv[index == k] <- predict(object = fit0.gam, newdata = datatot_valid, type = "response")
        }

        if (ps.method == "glm") {
          pscv[index == k] <- glm.simplereg.ps(x.ps = x_ps_train, trt = trt_train, xnew = x_ps_valid, minPS = minPS, maxPS = maxPS)
        } else {
          pscv[index == k] <- glm.ps(x.ps = x_ps_train, trt = trt_train, xnew = x_ps_valid, minPS = minPS, maxPS = maxPS)
        }

      }#end of Kfold loops

      if ("twoReg" %in% score.method) {
        ## bb-th cross fitting two regression estimator
        beta1.final <- onearmglmcount.dr(y = y, x.cate = x.cate, time = time, trt = trt, ps = pscv, f.predictor = f1.predictcv)
        beta0.final <- onearmglmcount.dr(y = y, x.cate = x.cate, time = time, trt = 1 - trt, ps = 1 - pscv, f.predictor = f0.predictcv)
        delta.twoReg.mat[bb, ] <- as.vector(beta1.final - beta0.final)
      }#end of if ("twoReg" %in% score.method)

      if ("contrastReg" %in% score.method) {
        ## bb-th cross fitting contrast regression estimator
        fit_two <- twoarmglmcount.dr(y = y, x.cate = x.cate, time = time, trt = trt, ps = pscv,
                                     f1.predictor = f1.predictcv, f0.predictor = f0.predictcv,
                                     error.maxNR = error.maxNR, max.iterNR = max.iterNR, tune = tune)
        delta.contrastReg.mat[bb, ] <- fit_two$coef
        converge[bb] <- fit_two$converge
        if (converge[bb]) sigma.contrastReg.mat <- sigma.contrastReg.mat + fit_two$vcov
      }#end of if ("contrastReg" %in% score.method)
    }#end of B loops
  }#end of if c("twoReg", "contrastReg") %in% score.method

  if ("boosting" %in% score.method) {
    ## Boosting method based on the entire data (score 1)
    data1 <- datatot[trt == 1, ]
    # if model has a single predictor, GBM must have cv.folds = 0 https://github.com/zoonproject/zoon/issues/130
    cate.cvfold <- ifelse(ncol(x.cate) == 1, 0, 5)
    cate.gbmmethod <- ifelse(ncol(x.cate) == 1, "OOB", "cv")
    if (cate.gbmmethod == "OOB") warning("If the model has a single predictor, GBM must use OOB")

    fit1.boosting <- gbm(y ~ . - time + offset(time), data = data1, distribution = "poisson",
                           interaction.depth = tree.depth, n.trees = n.trees.boosting, cv.folds = cate.cvfold, ...)
    best1.iter <- max(10, gbm.perf(fit1.boosting, method = cate.gbmmethod, plot.it = plot.gbmperf))

    data0 <- datatot[trt == 0, ]
    fit0.boosting <- gbm(y ~ . - time + offset(time), data = data0, distribution = "poisson",
                           interaction.depth = tree.depth, n.trees = n.trees.boosting, cv.folds = cate.cvfold, ...)
    best0.iter <- max(10, gbm.perf(fit0.boosting, method = cate.gbmmethod, plot.it = plot.gbmperf))

    result$result.boosting <- list(fit0.boosting = fit0.boosting, best0.iter = best0.iter, fit1.boosting = fit1.boosting, best1.iter = best1.iter)
    best.iter <- max(best.iter, best1.iter, best0.iter)
  }

  if ("poisson" %in% score.method) {
    ## Naive Poisson regression method  (score 2)
    beta1.ini <- glm(y ~ x.cate + offset(time), family = "poisson", subset = (trt == 1))$coef
    beta0.ini <- glm(y ~ x.cate + offset(time), family = "poisson", subset = (trt == 0))$coef
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
    if (converge.contrastReg) {
      delta.contrastReg <- colMeans(delta.contrastReg.mat[converge == TRUE, , drop = FALSE])
      sigma.contrastReg <- sigma.contrastReg.mat/sum(converge)
    } else {
      delta.contrastReg <- colMeans(delta.contrastReg.mat)
      sigma.contrastReg <- sigma.contrastReg.mat
    }
    names(delta.contrastReg) <- colnames(sigma.contrastReg) <- rownames(sigma.contrastReg) <- c("(Intercept)", colnames(x.cate))
    result$result.contrastReg <- list(delta.contrastReg = delta.contrastReg,
                                      sigma.contrastReg = sigma.contrastReg,
                                      converge.contrastReg = converge.contrastReg)
  }


  if ("negBin" %in% score.method) {
    ## Naive negative binomial regression method  (score 5)
    beta1.ini <- glm.nb(y ~ x.cate + offset(time), subset = (trt == 1), maxit = 500)$coef
    beta0.ini <- glm.nb(y ~ x.cate + offset(time), subset = (trt == 0), maxit = 500)$coef
    delta.negBin <- beta1.ini - beta0.ini
    names(delta.negBin) <- c("(Intercept)", colnames(x.cate))
    result$result.negBin <- delta.negBin
  }

  result$best.iter <- best.iter
  result$fgam <- fgam

  return(result)
}



#' Calculate the log CATE score given the baseline covariates and follow-up time for specified scoring method methods
#'
#' Based on intxcount results of the CATE coefficients estimated with boosting, naive Poisson, two regression, contrast regression, negative binomial
#'
#' @param fit List of objects generated from intxcount: outputs of boosting, naive Poisson, two regression, contrast regression, negative binomial
#' @param x.cate Matrix of \code{p.cate} baseline covariates; dimension \code{n} (observations) by \code{p.cate} (covariates in the outcome model)
#' @param time Log-transformed person-years of follow-up; vector of size \code{n}
#' @param score.method A vector of one or multiple methods to estimate the CATE score.
#' Allowed values are: \code{'boosting'}, \code{'poisson'}, \code{'twoReg'}, \code{'contrastReg'},
#' \code{'negBin'}. Default specifies all 5 methods.
#'
#' @return score.boosting: Estimated log CATE score for all \code{n} observations with the boosting method; vector of size \code{n}
#'         score.poisson: Estimated log CATE score for all \code{n} observations with the naive Poisson method; vector of size \code{n}
#'         score.twoReg: Estimated log CATE score for all \code{n} observations with the two regression method; vector of size \code{n}
#'         score.contrastReg: Estimated log CATE score for all \code{n} observations with the contrast regression method; vector of size \code{n}
#'         score.negBin: Estimated log CATE score for all \code{n} observations with the naive Poisson method; vector of size \code{n}
#'         score = NA if the corresponding method is not called

scorecount <- function(fit, x.cate, time,
                       score.method = c("boosting", "poisson", "twoReg", "contrastReg", "negBin")) {

  result <- vector("list", length(score.method))
  names(result) <- paste0("score.", score.method)
  x.aug <- cbind(1, x.cate)

  if ("boosting" %in% score.method) {
    fit0.boosting <- fit$result.boosting$fit0.boosting
    best0.iter <- fit$result.boosting$best0.iter
    fit1.boosting <- fit$result.boosting$fit1.boosting
    best1.iter <- fit$result.boosting$best1.iter

    datanew <- data.frame(x = x.cate, time = median(time))
    colnames(datanew) <- c(colnames(x.cate), "time")
    withCallingHandlers({
      predict0 <- predict(object = fit0.boosting, newdata = datanew, n.trees = best0.iter)
      predict1 <- predict(object = fit1.boosting, newdata = datanew, n.trees = best1.iter)
    },
    warning = function(w) {
      if (grepl("does not add the offset", conditionMessage(w)))
        invokeRestart("muffleWarning") # suppress warning: "predict.gbm does not add the offset to the predicted values."
    })
    result$score.boosting <- predict1 - predict0
  }

  if ("poisson" %in% score.method) {
    delta.poisson <- fit$result.poisson
    result$score.poisson <- as.numeric(as.matrix(x.aug[,is.na(delta.poisson) == FALSE]) %*% delta.poisson[is.na(delta.poisson) == FALSE])
  }

  if ("twoReg" %in% score.method) {
    delta.twoReg <- fit$result.twoReg
    result$score.twoReg <- as.numeric(as.matrix(x.aug[,is.na(delta.twoReg) == FALSE]) %*% delta.twoReg[is.na(delta.twoReg) == FALSE])
  }

  if ("contrastReg" %in% score.method) {
    delta.contrastReg <- fit$result.contrastReg$delta.contrastReg
    result$score.contrastReg <- as.numeric(as.matrix(x.aug[,is.na(delta.contrastReg) == FALSE]) %*% delta.contrastReg[is.na(delta.contrastReg) == FALSE])
  }

  if ("negBin" %in% score.method) {
    delta.negBin <- fit$result.negBin
    result$score.negBin <- as.numeric(as.matrix(x.aug[,is.na(delta.negBin) == FALSE]) %*% delta.negBin[is.na(delta.negBin) == FALSE])
  }

  return(result)
}
