#' Doubly robust estimator of and inference for the average treatment effect for count, survival and continuous data
#'
#' Doubly robust estimator of the average treatment effect between two treatments, which is the rate ratio
#' for count outcomes, the restricted mean time lost ratio for survival outcomes and the mean difference for continuous outcome. Bootstrap is used for
#' inference.
#'
#' @param response A string describing the type of outcome in the data. Allowed values include
#' "count" (see \code{\link{catecvcount}()}), "survival" (see
#' \code{\link{catecvsurv}()}) and "continuous" (see \code{\link{catecvmean}()}).
#' @param cate.model A formula describing the outcome model to be fitted.
#' The outcome must appear on the left-hand side. For survival outcomes,
#' a \code{Surv} object must be used to describe the outcome.
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
#' be printed. \code{1} indicates messages are printed and \code{0} otherwise. Default is \code{0}.
#'
#' @return For count response, see description of outputs in \code{\link{atefitcount}()}.
#' For survival response, see description of outputs in \code{\link{atefitsurv}()}.
#'
#' @details For count response, see details in \code{\link{atefitcount}()}.
#' For survival response, see details in \code{\link{atefitsurv}()}.
#'
#' @examples
#'
#' # Count outcome
#' output <- atefit(response = "count",
#'                  data = countExample,
#'                  cate.model = y ~ age + female + previous_treatment +
#'                               previous_cost + previous_number_relapses +
#'                               offset(log(years)),
#'                  ps.model = trt ~ age + previous_treatment,
#'                  verbose = 1,
#'                  n.boot = 50, seed = 999)
#' output
#' plot(output)
#'
#' \dontrun{
#'
#' # Survival outcome
#' tau0 <- with(survivalExample,
#'                  min(quantile(y[trt == "drug1"], 0.95), quantile(y[trt == "drug0"], 0.95)))
#'
#' output2 <- atefit(response = "survival",
#'                   data = survivalExample,
#'                   cate.model = survival::Surv(y, d) ~ age + female +
#'                         previous_cost + previous_number_relapses,
#'                         ps.model = trt ~ age + previous_treatment,
#'                   tau0 = tau0,
#'                   seed = 999)
#' output2
#' plot(output2)
#'
#'# Continuous outcome
#' output3 <- atefit(response = "continuous",
#'                        cate.model = y ~ age +
#'                                         previous_treatment +
#'                                         previous_cost +
#'                                         previous_status_measure,
#'                        previous_cost + previous_number_relapses + offset(log(years)),
#'                        ps.model = trt ~ previous_treatment,
#'                        data = meanExample,
#'                        seed = 999)
#' print(output3)
#' output$plot
#'}
#' @export
#'
#' @importFrom ggplot2 ggplot geom_histogram geom_vline
#' @importFrom dplyr mutate
#' @importFrom tidyr gather

atefit <- function(response,
                   data,
                   cate.model,
                   ps.model,
                   ps.method = "glm",
                   ipcw.model = NULL,
                   ipcw.method = "breslow",
                   minPS = 0.01,
                   maxPS = 0.99,
                   verbose = 0,
                   followup.time = NULL,
                   tau0 = NULL,
                   surv.min = 0.025,
                   interactions = TRUE,
                   n.boot = 500,
                   seed = NULL) {

  stopifnot("`response` must be either `count`, `survival`, or `continuous`." = any(response== c("count", "survival", "continuous")))
  .args <- as.list(match.call())[-1]
  .args$response <- NULL
  switch(response,
         count = {
           atefitout <- do.call(atefitcount, .args)
         },
         survival = {
           atefitout <- do.call(atefitsurv, .args)
         },
         continuous = {
           atefitout <- do.call(atefitmean, .args)
         }
  )
  return(atefitout)

}

#' Histogram of bootstrap estimates
#'
#' @param x An object of class \code{"atefit"}.
#' @param bins Number of bins
#' @param alpha Opacity
#' @param title The text for the title
#' @param theme Defaults to \code{theme_classic()}. Other options include \code{theme_grey()}, \code{theme_bw()}, \code{theme_light()}, \code{theme_dark()}, and \code{theme_void()}
#' @param ... Other parameters
#'
#' @details Create a histogram displaying the distribution of the bootstrap estimates.
#' The red vertical reference line represents the final estimate.
#' @author Thomas Debray
#'
#' @importFrom dplyr %>% add_row mutate
#' @importFrom ggplot2 ggplot geom_histogram geom_vline theme_classic
#' @importFrom tidyr gather
#' @export
#'
plot.atefit <- function(x, bins, alpha = 0.7, title, theme = theme_classic(), ...) {

  if (missing(bins)) {
    bins <- max(1, min(50, floor(x$n.boot/10))) # need at least one bin
  }
  if (missing(title)) {
    title <- paste0(x$n.boot, " bootstrap iterations")
  }

  if (x$response == "count") {
    x$trt.boot %>% as.data.frame() %>%
      ggplot(aes(x = .data$.)) +
      geom_histogram(bins = bins, alpha = alpha) +
      theme +
      geom_vline(xintercept = x$log.rate.ratio$estimate, linetype = 1, color = "red") +
      labs(x = "Bootstrap values of the log rate ratio", y = "Frequency",
           title = title)
  } else if (x$response == "survival") {

    est.df <- data.frame(key = character(), value = numeric())
    est.df <- est.df %>% add_row(key = "rmst1", value = x$rmst1$estimate)
    est.df <- est.df %>% add_row(key = "rmst0", value = x$rmst0$estimate)
    est.df <- est.df %>% add_row(key = "log.rmtl.ratio", value = x$log.rmtl.ratio$estimate)
    est.df <- est.df %>% add_row(key = "log.hazard.ratio", value = x$log.hazard.ratio$estimate)

    x$trt.boot %>% as.data.frame() %>%
      gather(key = "key", value = "value") %>% na.omit() %>%
      mutate(key = factor(.data$key, levels = c("rmst1", "rmst0", "log.rmtl.ratio", "log.hazard.ratio"))) %>%
      ggplot(aes(x = .data$value)) + geom_histogram(bins = bins, alpha = alpha) +
      facet_wrap(. ~ .data$key, scales = "free") +
      geom_vline(data = est.df, aes(xintercept = .data$value), linetype = 1, color = "red") +
      theme +
      labs(x = "Bootstrap values", y = "Frequency", title = title)
  }

}

#' Print function for atefit
#'
#' @param x An object of class \code{"atefit"}.
#' @param ... Other parameters
#'
#' @details Display the estimated treatment effects for survival outcomes (log
#' restricted mean time lost ratio and log hazard ratio) and count outcomes
#' (the log rate ratio).
#'
#' @author Thomas Debray
#'
#' @export
#'
print.atefit <- function(x, ...) {

  if (x$response == "count") {
    cat(paste0("Average treatment effect:\n\n"))
    print(x$log.rate.ratio)
    cat(paste0("\nEstimated event rates:\n\n"))
    print(rbind(x$rate1,x$rate0))
  }
  if (x$response == "survival") {
    cat(paste0("Average treatment effect:\n\n"))
    print(rbind(x$log.rmtl.ratio,x$log.hazard.ratio))
    cat(paste0("\nEstimated RMST:\n\n"))
    print(rbind(x$rmst1,x$rmst0))
  }

  if (!is.null(x$warning)) {
    cat("\n")
    warning(x$warning)
  }

}
