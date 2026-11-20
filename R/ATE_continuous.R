#' Doubly robust estimator of and inference for the average treatment effect for
#' continuous data
#'
#' Continuous-outcome average treatment effect estimation is not yet implemented.
#'
#' @param data A data frame containing the variables in the outcome and
#' treatment models; a data frame with \code{n} rows (1 row per
#' observation).
#' @param outcome_model A formula specifying the outcome model, with the
#' continuous outcome on the left and baseline covariates on the right. Do not
#' include treatment. Continuous-outcome fitting is not yet implemented.
#' @param ps_model A formula describing the treatment model (propensity-score model) to be fitted.
#' The treatment must appear on the left-hand side. The treatment must be a
#' numeric vector coded as 0/1. If data are from a RCT, specify \code{ps_model}
#' as an intercept-only model.
#' @param ps_method A character value for the method to estimate the propensity
#' score. Allowed values include one of: \code{'glm'} for logistic regression
#' with main effects only (default), or \code{'lasso'} for a logistic regression
#' with main effects and LASSO penalization on two-way interactions (added to
#' the model if interactions are not specified in \code{ps_model}). Relevant
#' only when \code{ps_model} has more than one variable.
#' @param minPS A numerical value between 0 and 1 below which estimated
#' propensity scores should be truncated. Default is \code{0.01}.
#' @param maxPS A numerical value between 0 and 1 above which estimated
#' propensity scores should be truncated. Must be strictly greater than
#' \code{minPS}. Default is \code{0.99}.
#' @param n.boot A numeric value indicating the number of bootstrap samples
#' used. Default is \code{500}.
#' @param plot.boot A logical value indicating whether histograms of the
#' bootstrapped treatment effect estimates should be produced at every
#' \code{n.boot/10}-th iteration and whether the final histogram should be
#' outputted. Default is \code{FALSE}.
#' @param seed An optional integer specifying an initial randomization seed for
#' reproducibility. Default is \code{NULL}, corresponding to no seed.
#'
#' @param verbose An integer value indicating whether intermediate progress
#' messages and histograms should be printed. \code{1} indicates messages are
#' printed and \code{0} otherwise. Default is \code{0}.
#' @param ps.control Named options for LASSO propensity-score interaction screening.
#' Supported entry: \code{solver = "C"} (default) or \code{solver = "R"}.
#' Only affects screening when \code{ps_method = "lasso"}; the LASSO fit and
#' interaction selection rule are unchanged. Kept separate from boosting's
#' \code{...} and other numerical controls.
#' @param ... Additional arguments.
#' Deprecated argument names from precmed 1.1.0 remain supported for
#' backward compatibility and emit warnings.
#' @return Continuous-outcome fitting is not implemented; this function stops
#' with an error and does not return an object.
#'
#' @details Continuous-outcome estimation is not currently available.
#'
#' @examples
#' # This module is not implemented yet!
#'
#' @importFrom ggplot2 ggplot geom_histogram geom_vline
#' @importFrom dplyr mutate

atefitmean <- function(data,
                       outcome_model,
                       ps_model,
                       ps_method = "glm",
                       minPS = 0.01,
                       maxPS = 0.99,
                       n.boot = 500,
                       plot.boot = FALSE,
                       seed = NULL,
                       verbose = 0, ps.control = list(), ...) {
  .args <- as.list(match.call())[-1L]
  aliases <- intersect(names(.args), c("cate.model", "ps.model", "ps.method"))
  .args <- resolve_ate_aliases(.args, reject_unknown = TRUE,
                               allowed = names(formals(sys.function())))
  if ("cate.model" %in% aliases) outcome_model <- eval(.args$outcome_model, parent.frame())
  if ("ps.model" %in% aliases) ps_model <- eval(.args$ps_model, parent.frame())
  if ("ps.method" %in% aliases) ps_method <- eval(.args$ps_method, parent.frame())
  if (!missing(ps_model)) ps.model <- ps_model
  ps.method <- ps_method
  if (missing(outcome_model)) stop("outcome_model must be supplied.")
  cate.model <- outcome_model

  resolve_ps_control(ps.control)


  stop("This functionality is not implemented yet")
  interactions <- TRUE

  # Set seed once for reproducibility
  set.seed(seed)

  # CHECK ARGUMENTS
  # TODO: Phoebe I think we removed the arg.check on interactions but here it is?
  arg.checks(fun = "drinf", response = "continuous", data = data,
             ps.method = ps.method,
             minPS = minPS,
             maxPS = maxPS,
             interactions = interactions,
             n.boot = n.boot,
             plot.boot = plot.boot)

  # PRE-PROCESSING
  preproc <- data.preproc.mean(fun = "drinf", cate.model = cate.model,
                               init.model = NULL, #TODO: Needs to be specified
                               ps.model = ps.model,
                               data = data, ps.method = ps.method)
  y <- preproc$y
  trt <- preproc$trt
  x.ps <- preproc$x.ps
  x.cate <- preproc$x.cate

  # Point estimate
  est <- drmean(y = y, x.cate = x.cate, x.ps = x.ps, trt = trt, ps.method = ps.method, minPS = minPS, maxPS = maxPS, interactions = interactions, ps.control = ps.control)
  meandiff <- est$mean.diff

  # Apply bootstrap
  n <- length(y)
  if (is.na(meandiff) == TRUE) {
    stop("Impossible to use bootstrap, Mean Difference is NA.")
  } else {
    trt.boot <- rep(0, n.boot)
    for (i in 1:n.boot) {
      idsub.boot <- sample(n, size = n, replace = TRUE)
      trt.boot[i] <- drmean(y = y[idsub.boot], x.cate = x.cate[idsub.boot, , drop = FALSE], x.ps = x.ps[idsub.boot, , drop = FALSE], trt = trt[idsub.boot], ps.method = ps.method, minPS = minPS, maxPS = maxPS, interactions = interactions, ps.control = ps.control)$mean.diff

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

#' Doubly robust estimator of the average treatment effect for continuous data
#'
#' Doubly robust estimator of the average treatment effect between two
#' treatments, which is the mean difference of treatment 1 over treatment 0 for
#' continuous outcomes.
#'
#' @param y A numeric vector of size \code{n} with each element representing
#' the observed continuous outcome for each subject.
#' @param trt A numeric vector (in \{0, 1\}) of size \code{n} with each element
#' representing the treatment received for each subject.
#' @param x.cate A numeric matrix of dimension \code{n} by \code{p.cate} with
#' each column representing each baseline covariate specified in the outcome
#' model for all subjects.
#' @param x.ps A numeric matrix of dimension \code{n} by \code{p.ps + 1} with
#' a leading column of 1 as the intercept and each remaining column representing
#' each baseline covariate specified in the treatment model for all
#' subjects
#' @param ps.method A character value for the method to estimate the propensity
#' score. Allowed values include one of:
#' \code{'glm'} for logistic regression with main effects only (default), or
#' \code{'lasso'} for a logistic regression with main effects and LASSO
#' penalization on two-way interactions (added to the model if interactions are
#' not specified in \code{ps.model}). Relevant only when \code{ps.model} has
#' more than one variable.
#' @param minPS A numerical value between 0 and 1 below which estimated propensity
#' scores should be truncated. Default is \code{0.01}.
#' @param maxPS A numerical value between 0 and 1 above which estimated propensity
#' scores should be truncated. Must be strictly greater than \code{minPS}.
#' Default is \code{0.99}.
#' @param interactions A logical value indicating whether the outcome model
#' should assume interactions between \code{x} and \code{trt}. If \code{TRUE},
#' interactions will be assumed only if at least 10 patients received each
#' treatment option. Default is \code{TRUE}.
#'
#' @param ps.control Named options for LASSO propensity-score interaction screening.
#' Supported entry: \code{solver = "C"} (default) or \code{solver = "R"}.
#' Only affects screening when \code{ps.method = "lasso"}; the LASSO fit and
#' interaction selection rule are unchanged. Kept separate from boosting's
#' \code{...} and other numerical controls.
#' @return Return a list of 4 elements:
#' \itemize{
#'   \item{\code{mean.diff}: } A numeric value of the estimated mean difference.
#'   \item{\code{mean.diff0}: } A numeric value of the estimated mean difference
#'   in treatment group 0.
#'   \item{\code{mean.diff1}: } A numeric value of the estimated mean difference
#'   in treatment group 1.
#' }
#'

drmean <- function(y, trt, x.cate, x.ps,
                   ps.method = "glm", minPS = 0.01, maxPS = 0.99,
                   interactions = TRUE, ps.control = list()) {
  ps.args <- resolve_ps_control(ps.control)


  n <- length(y)

  withCallingHandlers({
    if (ps.method == "lasso") {
      ps <- glm.ps(trt = trt, x.ps = x.ps, xnew = NULL, minPS = minPS, maxPS = maxPS, solver = ps.args$solver)
    } else {
      ps <- glm.simplereg.ps(trt = trt, x.ps = x.ps, xnew = NULL, minPS = minPS, maxPS = maxPS)
    }
  },
  warning = function(w) {
    if (grepl("numerically 0|1 occurred", conditionMessage(w)) | grepl("glm.fit: algorithm did not converge", conditionMessage(w)))
      invokeRestart("muffleWarning")  # suppress warning: "glm.fit: fitted probabilities/rates numerically 0 occurred." or "glm.fit: algorithm did not converge."
  })

  if (min(sum(trt), n - sum(trt)) > 10 & interactions == TRUE) {
    withCallingHandlers({
      # Calculate rate by treatment arm by fitting two separate models
      fit1 <- glm(y ~ x.cate , family = "gaussian", subset = (trt == 1))
      fit0 <- glm(y ~ x.cate , family = "gaussian", subset = (trt == 0))
    },
    warning = function(w) {
      if (grepl("numerically 0|1 occurred|algorithm did not converge", conditionMessage(w)))
        invokeRestart("muffleWarning")  # suppress warning: "glm.fit: fitted rates numerically 0 occurred." or "glm.fit: algorithm did not converge."
    })
    beta1 <- fit1$coef
    y1hat <- cbind(1, x.cate)[, is.na(beta1) == FALSE, drop = FALSE] %*% beta1[is.na(beta1) == FALSE]
    beta0 <- fit0$coef
    y0hat <- cbind(1, x.cate)[, is.na(beta0) == FALSE, drop = FALSE] %*% beta0[is.na(beta0) == FALSE]
  } else {
    withCallingHandlers({
      # Calculate rate by treatment arm by fitting a single model
      fit <- glm(y ~ x.cate + trt , family = "gaussian")
    },
    warning = function(w) {
      if (grepl("numerically 0|1 occurred|algorithm did not converge", conditionMessage(w)))
        invokeRestart("muffleWarning")  # suppress warning: "glm.fit: fitted rates numerically 0 occurred." or "glm.fit: algorithm did not converge."
    })
    beta <- fit$coef
    y1hat <- cbind(1, x.cate, 1)[, is.na(beta) == FALSE, drop = FALSE] %*% beta[is.na(beta) == FALSE]
    y0hat <- cbind(1, x.cate, 0)[, is.na(beta) == FALSE, drop = FALSE] %*% beta[is.na(beta) == FALSE]
  }
  mu1dr <- sum(y1hat + (y - y1hat) * trt / ps) / n
  mu0dr <- sum(y0hat + (y - y0hat) * (1 - trt) / (1 - ps)) / n

  ###
  if (is.na(mu1dr - mu0dr) == FALSE) {
    md <- mu1dr - mu0dr
  } else {
    md <- NA
  }

  return(list(mean.diff = md, mean.diff0 = mu0dr, mean.diff1 = mu1dr))
}

# Matrix-only equivalent of drmean(..., ps.method = "glm", interactions = TRUE).
# outcome.design includes the intercept and is constructed once per split search.
# Author: Thomas Debray <tdebray@fromdatatowisdom.com>
drmean.matrix <- function(y, trt, outcome.design, x.ps, families,
                          minPS = 0.01, maxPS = 0.99) {
  withCallingHandlers({
    beta.ps <- stats::glm.fit(x.ps, trt, family = families$propensity,
                              intercept = FALSE)$coefficients
    keep <- !is.na(beta.ps)
    ps <- as.vector(1 / (1 + exp(-x.ps[, keep, drop = FALSE] %*% beta.ps[keep])))
    ps <- pmax(pmin(ps, maxPS), minPS)

    predict.mean <- function(beta, design) {
      keep <- !is.na(beta)
      design[, keep, drop = FALSE] %*% beta[keep]
    }
    if (min(sum(trt), length(y) - sum(trt)) > 10) {
      id1 <- trt == 1
      id0 <- trt == 0
      beta1 <- stats::glm.fit(outcome.design[id1, , drop = FALSE], y[id1],
                               family = families$outcome)$coefficients
      beta0 <- stats::glm.fit(outcome.design[id0, , drop = FALSE], y[id0],
                               family = families$outcome)$coefficients
      y1hat <- predict.mean(beta1, outcome.design)
      y0hat <- predict.mean(beta0, outcome.design)
    } else {
      beta <- stats::glm.fit(cbind(outcome.design, trt), y,
                              family = families$outcome)$coefficients
      y1hat <- predict.mean(beta, cbind(outcome.design, 1))
      y0hat <- predict.mean(beta, cbind(outcome.design, 0))
    }
    mu1dr <- sum(y1hat + (y - y1hat) * trt / ps) / length(y)
    mu0dr <- sum(y0hat + (y - y0hat) * (1 - trt) / (1 - ps)) / length(y)
    md <- if (!is.na(mu1dr - mu0dr)) mu1dr - mu0dr else NA
    list(mean.diff = md, mean.diff0 = mu0dr, mean.diff1 = mu1dr)
  }, warning = function(w) {
    if (grepl("numerically 0|1 occurred|glm.fit: algorithm did not converge|glm.fit: fitted rates numerically 0 occurred",
              conditionMessage(w))) invokeRestart("muffleWarning")
  })
}



#' Estimate the ATE of the mean difference in multiple bi-level subgroups
#' defined by the proportions
#'
#' If only care about the higher subgroup (above cutoff), only need
#' trt.est.high so set \code{onlyhigh} to be TRUE. Scores are adjusted to the
#' opposite sign if \code{higher.y} == FALSE; scores stay the same if
#' \code{higher.y} == TRUE. This is because \code{estcount.bilevel.subgroups}()
#' always takes the subgroup of the top highest adjusted scores,and higher
#' adjusted scores should always represent high responders in treatment group 1.
#'
#' @param y Observed outcome; vector of size \code{n} (observations)
#' @param trt Treatment received; vector of size \code{n} units with treatment coded as 0/1
#' @param x.cate Matrix of \code{p.cate} baseline covariates; dimension \code{n} by \code{p.cate} (covariates in the outcome model)
#' @param x.ps Matrix of \code{p.ps} baseline covariates (plus a leading column of 1 for the intercept);
#' dimension \code{n} by \code{p.ps + 1} (covariates in the treatment model plus intercept)
#' @param score Estimated CATE scores for all \code{n} observations from one of the six methods
#' (boosting, linear regression, two regressions, contrast regression, random forest, generalized additive model); vector of size \code{n}
#' @param prop Proportions corresponding to percentiles in the estimated CATE scores that define subgroups to calculate ATE for;
#' vector of floats in `(0, 1]` (if onlyhigh=T) or in `(0, 1)` (if onlyhigh=F):
#'              Each element of \code{prop} represents the high/low cutoff in each bi-level subgroup and the length of \code{prop}
#'              is number of bi-level subgroups
#' @param onlyhigh Indicator of returning only the ATEs in the higher-than-cutoff category of the bi-level subgroups; boolean
#' @param higher.y A logical value indicating whether higher (\code{TRUE}) or lower (\code{FALSE})
#' values of the outcome are more desirable. Default is \code{TRUE}.
#' @param ps.method A character value for the method to estimate the propensity score. Allowed values include one of:
#' \code{'glm'} for logistic regression with main effects only (default), or
#' \code{'lasso'} for a logistic regression with main effects and LASSO penalization on
#' two-way interactions (added to the model if interactions are not specified in \code{ps.model}).
#' Relevant only when \code{ps.model} has more than one variable.
#' @param minPS A numerical value (in `[0, 1]`) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in `(0, 1]`) above which estimated propensity scores should be
#' truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.
#'
#' @param ps.control Named options for LASSO propensity-score interaction screening.
#' Supported entry: \code{solver = "C"} (default) or \code{solver = "R"}.
#' Only affects screening when \code{ps.method = "lasso"}; the LASSO fit and
#' interaction selection rule are unchanged. Kept separate from boosting's
#' \code{...} and other numerical controls.
#' @return ate.est.high: estimated ATEs in the multiple bi-level subgroups that are in the higher-than-cutoff category;
#' vector of size equal to the length of prop; always returned
#'         ate.est.low: estimated ATEs in the multiple bi-level subgroups that are in the lower-than-cutoff category;
#'         vector of size equal to the length of prop; returned only when \code{onlyhigh} == TRUE
estmean.bilevel.subgroups <- function(y, x.cate, x.ps, trt, score, higher.y,
                                      prop, onlyhigh,
                                      ps.method = "glm", minPS = 0.01, maxPS = 0.99, ps.control = list()) {

  if (higher.y == FALSE) score <- -score
  cut <- quantile(score, 1 - prop)
  n.subgroup <- length(prop)
  ate.est.high <- rep(0, n.subgroup)
  if (onlyhigh == FALSE) ate.est.low <- rep(0, n.subgroup)
  for (b in 1:n.subgroup) {
    idsub1 <- (score >= cut[b])
    ate.est.high[b] <- drmean(y = y[idsub1], x.cate = x.cate[idsub1, , drop = FALSE], x.ps = x.ps[idsub1, , drop = FALSE],
                              trt = trt[idsub1], ps.method = ps.method, minPS = minPS, maxPS = maxPS,
                              interactions = FALSE, ps.control = ps.control)$mean.diff
    if (onlyhigh == FALSE) {
      idsub0 <- (score < cut[b])
      ate.est.low[b] <- drmean(y = y[idsub0], x.cate = x.cate[idsub0, , drop = FALSE], x.ps = x.ps[idsub0, , drop = FALSE],
                               trt = trt[idsub0], ps.method = ps.method, minPS = minPS, maxPS = maxPS,
                               interactions = FALSE, ps.control = ps.control)$mean.diff
    }
  }
  if (onlyhigh == TRUE) {
    return((ate.est.high))
  } else {
    return(list(ate.est.high = (ate.est.high), ate.est.low = (ate.est.low)))
  }
}


#' Estimate the ATE of the mean difference in one multilevel subgroup defined by the proportions
#'
#' Scores are adjusted to the opposite sign if \code{higher.y} == FALSE; scores stay the same if \code{higher.y} == TRUE;
#'  this is because subgroups defined in estmean.multilevel.subgroup() start from the lowest to the highest adjusted scores,
#'  and higher adjusted scores should always represent high responders of trt=1
#'
#' @param y Observed outcome; vector of size \code{n} (observations)
#' @param trt Treatment received; vector of size \code{n} units with treatment coded as 0/1
#' @param x.cate Matrix of \code{p.cate} baseline covariates; dimension \code{n} by \code{p.cate} (covariates in the outcome model)
#' @param x.ps Matrix of \code{p.ps} baseline covariates (plus a leading column of 1 for the intercept);
#' dimension \code{n} by \code{p.ps + 1} (covariates in the treatment model plus intercept)
#' @param score Estimated CATE scores for all \code{n} observations from one of the six methods
#' (boosting, linear regression, two regressions, contrast regression, random forest, generalized additive model); vector of size \code{n}
#' @param prop Proportions corresponding to percentiles in the estimated CATE scores that define subgroups to calculate ATE for;
#' vector of floats in `[0, 1]` always starting with 0 and ending with 1:
#'              Each element of \code{prop} represents inclusive cutoffs in the multilevel subgroup and the length of \code{prop}
#'              is number of categories in the multilevel subgroup
#' @param higher.y A logical value indicating whether higher (\code{TRUE}) or lower (\code{FALSE})
#' values of the outcome are more desirable. Default is \code{TRUE}.
#' @param ps.method A character value for the method to estimate the propensity score. Allowed values include one of:
#' \code{'glm'} for logistic regression with main effects only (default), or
#' \code{'lasso'} for a logistic regression with main effects and LASSO penalization on
#' two-way interactions (added to the model if interactions are not specified in \code{ps.model}).
#' Relevant only when \code{ps.model} has more than one variable.
#' @param minPS A numerical value (in `[0, 1]`) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in `(0, 1]`) above which estimated propensity scores should be
#' truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.
#'
#' @param ps.control Named options for LASSO propensity-score interaction screening.
#' Supported entry: \code{solver = "C"} (default) or \code{solver = "R"}.
#' Only affects screening when \code{ps.method = "lasso"}; the LASSO fit and
#' interaction selection rule are unchanged. Kept separate from boosting's
#' \code{...} and other numerical controls.
#' @return estimated ATEs of all categories in the one multilevel subgroup; vector of size equal to the length of categories in the multilevel subgroup
estmean.multilevel.subgroup <- function(y, x.cate, x.ps, trt, score, higher.y, prop,
                                        ps.method = "glm", minPS = 0.01, maxPS = 0.99, ps.control = list()){

  if (higher.y == FALSE) score <- -score
  cut <- quantile(score, prop)
  n.subgroup <- length(prop) - 1
  ate.est <- rep(NA, n.subgroup)
  for (b in 1:n.subgroup) {
    idsub <- (score >= cut[b] & score <= cut[b + 1])
    ate.est[b] <- drmean(y = y[idsub], x.cate = x.cate[idsub, , drop = FALSE],
                         x.ps = x.ps[idsub, , drop = FALSE],
                         trt = trt[idsub],
                         ps.method = ps.method, minPS = minPS, maxPS = maxPS, interactions = FALSE, ps.control = ps.control)$mean.diff
  }
  return((ate.est))
}
