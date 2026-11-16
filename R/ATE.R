#' Estimate average treatment effects for binary, count and survival outcomes
#'
#' Doubly robust average treatment effect estimation for two treatments. The
#' effect is a risk difference for binary outcomes, a rate ratio for count
#' outcomes, or a restricted mean time lost ratio for survival outcomes.
#' Bootstrap resampling is used for inference. Continuous outcomes are not
#' yet implemented.
#'
#' @param response A string describing the type of outcome in the data. Allowed values include
#' "binary", "count" (see \code{\link{catecvcount}()}), "survival" (see
#' \code{\link{catecvsurv}()}) and "continuous" (see \code{\link{catecvmean}()}).
#' @param cate.model A formula describing the outcome model to be fitted.
#' The outcome must appear on the left-hand side. For survival outcomes,
#' a \code{Surv} object must be used to describe the outcome.
#' @param ps.model A formula describing the treatment model (propensity-score model) to be fitted. The treatment must
#' appear on the left-hand side. For binary outcomes, logical treatment values
#' \code{FALSE} and \code{TRUE} are mapped to 0 and 1, respectively.
#' Numeric 0/1 values or two distinct treatment labels are also accepted.
#' If data are from a randomized controlled trial, specify \code{ps.model = trt ~ 1} as an intercept-only model.
#' @param data A data frame containing the variables in the outcome, treatment,
#' and censoring models (if specified); a data frame with \code{n} rows (1 row per observation).
#' @param ipcw.model A formula describing the censoring model used for inverse
#' probability of censoring weighting (IPCW). The left-hand side must be empty. Only applies for survival outcomes.
#' Default is \code{NULL}, which uses a censoring model with the same covariates
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
#' @param minPS A numerical value (in `[0, 1]`) below which estimated propensity scores should be
#' truncated. Default is \code{0.01}.
#' @param maxPS A numerical value (in `(0, 1]`) above which estimated propensity scores should be
#' truncated. Must be strictly greater than \code{minPS}. Default is \code{0.99}.
#' @param interactions A logical value indicating whether the outcome model should assume interactions
#' between \code{x} and \code{trt}. For binary outcomes this controls all
#' treatment-covariate interactions in the logistic outcome model. For count outcomes, if \code{TRUE}, interactions will
#' be assumed only if at least 10 patients received each treatment option. Default is \code{TRUE}.
#' @param n.boot A numeric value indicating the number of bootstrap samples used. Default is \code{500}.
#' @param seed An optional integer specifying an initial randomization seed for reproducibility.
#' Default is \code{NULL}, corresponding to no seed.
#' @param verbose An integer value indicating whether intermediate progress messages and histograms should
#' be printed. \code{1} indicates messages are printed and \code{0} otherwise. Default is \code{0}.
#'
#' Binary outcomes use \code{response = "binary"} and estimate a doubly robust
#' risk difference (treatment 1 minus treatment 0) with a logistic outcome model.
#' The result includes \code{risk.difference}, \code{risk0}, \code{risk1}, and
#' bootstrap estimates. If any bootstrap fit fails, inference is unavailable.
#' @return For binary response, see \code{\link{atefitbinary}()}.
#' For count response, see description of outputs in \code{\link{atefitcount}()}.
#' For survival response, see description of outputs in \code{\link{atefitsurv}()}.
#'
#' @details For binary response, see \code{\link{atefitbinary}()}.
#' For count response, see details in \code{\link{atefitcount}()}.
#' For survival response, see details in \code{\link{atefitsurv}()}.
#'
#' @references
#' Hoogland, J., IntHout, J., Belias, M., et al. (2021). A tutorial on
#' individualized treatment effect prediction from randomized trials with a
#' binary endpoint. \emph{Statistics in Medicine}, 40(26), 5961--5981.
#' \doi{10.1002/sim.9154}
#'
#' @details Failed or non-finite bootstrap estimates are recorded in
#' \code{bootstrap.diagnostics}, with replicate number, affected measure and
#' reason. A warning summarizes failures; inference for each affected measure
#' is unavailable rather than based on a subset of replicates. See also
#' \code{diagnostics()} to inspect these records.
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
#'                  n.boot = 50,
#'                  seed = 999)
#' output
#' plot(output)
#'
#' \donttest{
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
#'}
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#' @author Xiaotong (Phoebe) Jiang \email{xiaotong.phoebe.jiang@gmail.com}
#' @export
#'
#' @importFrom ggplot2 ggplot geom_histogram geom_vline
#' @importFrom dplyr mutate
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
                   followup.time = NULL,
                   tau0 = NULL,
                   surv.min = 0.025,
                   interactions = TRUE,
                   n.boot = 500,
                   seed = NULL,
                   verbose = 0, rmst.control = list(), ps.control = list()) {
  check_response(response)
  resolve_ps_control(ps.control)


  if (length(rmst.control) && response != "survival")
    stop("rmst.control is currently supported only for survival outcomes.")
  .args <- as.list(match.call())[-1]
  .args$response <- NULL
  if (response != "survival") .args$rmst.control <- NULL
  switch(response,
         binary = {
           atefitout <- do.call(atefitbinary, .args, envir = parent.frame())
         },
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

# Author: Thomas Debray
precmed_plot_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "#E5E9ED", linewidth = 0.4),
      strip.background = ggplot2::element_rect(fill = "#F2F5F7", color = NA),
      strip.text = ggplot2::element_text(face = "bold", color = "#263746", hjust = 0,
                                         margin = ggplot2::margin(6, 8, 6, 8)),
      strip.text.y.right = ggplot2::element_text(angle = 0),
      axis.title = ggplot2::element_text(color = "#263746"),
      axis.text = ggplot2::element_text(color = "#52616B"),
      plot.caption = ggplot2::element_text(color = "#52616B", hjust = 0, margin = ggplot2::margin(t = 10))
    )
}

#' Histogram of bootstrap estimates
#'
#' @param x An object of class \code{"atefit"}.
#' @param bins Number of bins
#' @param alpha Histogram bar opacity, between 0 and 1. Defaults to 1.
#' @param title The text for the title
#' @param theme A ggplot2 theme. By default, a customized minimal theme is used;
#' supply, for example, \code{ggplot2::theme_classic()} to use another theme.
#' @param null.value Reference value on the plotted effect scale. Defaults to
#' 0 (no risk difference or a ratio of 1); use \code{NULL} or \code{NA} to omit the null line.
#' @param metrics Character vector of metrics to plot, in display order. For
#' binary fits, choose from \code{"risk0"}, \code{"risk1"}, and
#' \code{"risk.difference"}. For
#' count fits, choose from \code{"rate0"}, \code{"rate1"}, and
#' \code{"log.rate.ratio"}; for survival fits, choose from \code{"rmst0"},
#' \code{"rmst1"}, \code{"log.rmtl.ratio"}, and
#' \code{"log.cumulative.hazard.ratio"}. Older survival fits may instead
#' contain \code{"log.hazard.ratio"}. The default \code{NULL} plots all
#' metrics available in the object.
#' @param ... Other parameters
#'
#' @return A plot of the class \code{ggplot}, displaying the estimated ATE across
#' the bootstrap samples
#'
#' @details Create a histogram displaying the distribution of the bootstrap estimates.
#' The solid red vertical line represents the full-sample estimate. A dashed
#' black line marks \code{null.value} in risk-difference and log-ratio panels;
#' it is not shown in treatment-specific risk, rate or RMST panels.
#' Treatment-specific facets use the original treatment labels when available.
#' Older binary or count objects without treatment-specific bootstrap estimates
#' show only the risk difference or log rate ratio, respectively.
#' @author Thomas Debray
#'
#' @importFrom dplyr %>% add_row mutate
#' @importFrom ggplot2 ggplot geom_histogram geom_vline waiver
#' @importFrom tidyr gather
#' @export
#'
plot.atefit <- function(x, bins, alpha = 1, title = waiver(),
                        theme = precmed_plot_theme(),
                        null.value = 0, metrics = NULL, ...) {

  if (!is.null(null.value) && (length(null.value) != 1L ||
      !(is.numeric(null.value) || (is.logical(null.value) && is.na(null.value))) ||
      (!is.na(null.value) && !is.finite(null.value))))
    stop("null.value must be a finite numeric value, NA, or NULL.")
  show.null <- !is.null(null.value) && !is.na(null.value)

  if (missing(bins)) {
    bins <- max(1, min(50, floor(x$n.boot/10))) # need at least one bin
  }

  if (identical(x$response, "binary")) {
    has.risks <- !is.null(x$risk.boot) && all(c("risk0", "risk1") %in% names(x$risk.boot)) &&
      !is.null(x$risk0$estimate) && !is.null(x$risk1$estimate)
    available <- c(if (has.risks) c("risk0", "risk1"), "risk.difference")
    if (is.null(metrics)) metrics <- available
    if (!is.character(metrics) || !length(metrics) || anyNA(metrics))
      stop("metrics must be a non-empty vector of distinct available metric names: ", paste(available, collapse = ", "), ".")
    if (anyDuplicated(metrics) || any(!metrics %in% available))
      stop("metrics must be a non-empty vector of distinct available metric names: ", paste(available, collapse = ", "), ".")
    if (has.risks) {
      treatment.labels <- x$treatment.labels
      if (is.null(treatment.labels) || length(treatment.labels) != 2L ||
          anyNA(treatment.labels) || any(!nzchar(treatment.labels)))
        treatment.labels <- c("treatment 0", "treatment 1")
      facet.labels <- c(risk0 = paste0("Risk (", treatment.labels[1], ")"),
                        risk1 = paste0("Risk (", treatment.labels[2], ")"),
                        risk.difference = "Risk difference")
      boot.df <- data.frame(risk0 = x$risk.boot$risk0, risk1 = x$risk.boot$risk1,
                            risk.difference = x$trt.boot) %>%
        gather(key = "key", value = "value") %>%
        dplyr::filter(.data$key %in% metrics) %>%
        mutate(key = factor(.data$key, levels = metrics))
      estimates <- c(risk0 = x$risk0$estimate, risk1 = x$risk1$estimate,
                     risk.difference = x$risk.difference$estimate)
      est.df <- data.frame(key = factor(metrics, levels = metrics), value = unname(estimates[metrics]))
      p <- ggplot(boot.df, aes(x = .data$value)) +
        geom_histogram(bins = bins, alpha = alpha, fill = "#7293AD", color = "white", linewidth = 0.3) +
        facet_wrap(. ~ .data$key, nrow = 1, scales = "free",
                   labeller = ggplot2::as_labeller(facet.labels)) +
        geom_vline(data = est.df, aes(xintercept = .data$value), color = "#B43B35", linewidth = 0.8) +
        theme + labs(x = "Bootstrap estimate", y = "Frequency", title = title,
                     caption = if (show.null && "risk.difference" %in% metrics) "Red: full-sample estimate   |   Dashed: null value" else "Red: full-sample estimate")
      if (show.null && "risk.difference" %in% metrics) {
        null.df <- data.frame(key = factor("risk.difference", levels = metrics), value = null.value)
        p <- p + geom_vline(data = null.df, aes(xintercept = .data$value), linetype = 2, color = "#303030", linewidth = 0.7)
      }
      return(p)
    }
    p <- ggplot(data.frame(value = x$trt.boot), aes(x = .data$value)) +
      geom_histogram(bins = bins, alpha = alpha, fill = "#7293AD", color = "white") +
      geom_vline(xintercept = x$risk.difference$estimate, color = "#B43B35") +
      theme + labs(x = "Bootstrap risk difference", y = "Frequency", title = title)
    if (show.null) p <- p + geom_vline(xintercept = null.value, linetype = 2)
    return(p)
  }
  available <- if (x$response == "count") {
    c(if (!is.null(x$rate.boot)) c("rate0", "rate1"), "log.rate.ratio")
  } else if (x$response == "survival") {
    intersect(intersect(c("rmst1", "rmst0", "log.rmtl.ratio", "log.hazard.ratio",
                          "log.cumulative.hazard.ratio"), colnames(x$trt.boot)), names(x))
  } else character()
  if (is.null(metrics)) metrics <- available
  if (!is.character(metrics) || !length(metrics) || anyNA(metrics) ||
      anyDuplicated(metrics) || any(!metrics %in% available))
    stop("metrics must be a non-empty vector of distinct available metric names: ",
         paste(available, collapse = ", "), ".")

  if (x$response == "count") {
    if (!is.null(x$rate.boot)) {
      treatment.labels <- x$treatment.labels
      if (is.null(treatment.labels) || length(treatment.labels) != 2L ||
          anyNA(treatment.labels) || any(!nzchar(treatment.labels)))
        treatment.labels <- c("treatment 0", "treatment 1")
      keys <- metrics
      facet.labels <- c(rate0 = paste0("Rate (", treatment.labels[1], ")"),
                        rate1 = paste0("Rate (", treatment.labels[2], ")"),
                        log.rate.ratio = "Log rate ratio")
      boot.df <- data.frame(x$rate.boot, log.rate.ratio = x$trt.boot) %>%
        gather(key = "key", value = "value") %>%
        dplyr::filter(.data$key %in% keys) %>%
        mutate(key = factor(.data$key, levels = keys))
      estimates <- c(rate0 = x$rate0$estimate, rate1 = x$rate1$estimate,
                     log.rate.ratio = x$log.rate.ratio$estimate)
      est.df <- data.frame(key = factor(keys, levels = keys), value = unname(estimates[keys]))
      p <- ggplot(boot.df, aes(x = .data$value)) +
        geom_histogram(bins = bins, alpha = alpha, fill = "#7293AD", color = "white", linewidth = 0.3) +
        facet_wrap(. ~ .data$key, nrow = 1, scales = "free",
                   labeller = ggplot2::as_labeller(facet.labels)) +
        geom_vline(data = est.df, aes(xintercept = .data$value), linetype = 1, color = "#B43B35", linewidth = 0.8) +
        theme +
        labs(x = "Bootstrap estimate", y = "Frequency", title = title,
             caption = if (show.null && "log.rate.ratio" %in% keys) "Red: full-sample estimate   |   Dashed: null value" else "Red: full-sample estimate")
      if (show.null && "log.rate.ratio" %in% keys) {
        null.df <- data.frame(key = factor("log.rate.ratio", levels = keys), value = null.value)
        p <- p + geom_vline(data = null.df, aes(xintercept = .data$value), linetype = 2, color = "#303030", linewidth = 0.7)
      }
      return(p)
    }
    p <- x$trt.boot %>% as.data.frame() %>%
      ggplot(aes(x = .data$.)) +
      geom_histogram(bins = bins, alpha = alpha, fill = "#7293AD", color = "white", linewidth = 0.3) +
      theme +
      geom_vline(xintercept = x$log.rate.ratio$estimate, linetype = 1, color = "#B43B35", linewidth = 0.8) +
      labs(x = "Bootstrap values of the log rate ratio", y = "Frequency",
           title = title, caption = if (show.null) "Red: full-sample estimate   |   Dashed: null value" else "Red: full-sample estimate")
    if (show.null) p <- p + geom_vline(xintercept = null.value, linetype = 2, color = "#303030", linewidth = 0.7)
    p
  } else if (x$response == "survival") {

    treatment.labels <- x$treatment.labels
    if (is.null(treatment.labels) || length(treatment.labels) != 2L ||
        anyNA(treatment.labels) || any(!nzchar(treatment.labels)))
      treatment.labels <- c("treatment 0", "treatment 1")
    facet.labels <- c(rmst1 = paste0("RMST (", treatment.labels[2], ")"),
                      rmst0 = paste0("RMST (", treatment.labels[1], ")"),
                      log.rmtl.ratio = "Log RMTL ratio",
                      log.cumulative.hazard.ratio = "Log cumulative hazard ratio",
                      log.hazard.ratio = "Log hazard ratio")
    est.df <- data.frame(key = character(), value = numeric())
    est.df <- est.df %>% add_row(key = "rmst1", value = x$rmst1$estimate)
    est.df <- est.df %>% add_row(key = "rmst0", value = x$rmst0$estimate)
    est.df <- est.df %>% add_row(key = "log.rmtl.ratio", value = x$log.rmtl.ratio$estimate)
    if (!is.null(x$log.hazard.ratio))
      est.df <- est.df %>% add_row(key = "log.hazard.ratio", value = x$log.hazard.ratio$estimate)

    if (!is.null(x$log.cumulative.hazard.ratio))
      est.df <- est.df %>% add_row(key = "log.cumulative.hazard.ratio", value = x$log.cumulative.hazard.ratio$estimate)

    est.df <- est.df[match(metrics, est.df$key), , drop = FALSE]
    est.df$key <- factor(est.df$key, levels = metrics)

    p <- x$trt.boot %>% as.data.frame() %>%
      gather(key = "key", value = "value") %>% na.omit() %>%
      dplyr::filter(.data$key %in% metrics) %>%
      mutate(key = factor(.data$key, levels = metrics)) %>%
      ggplot(aes(x = .data$value)) + geom_histogram(bins = bins, alpha = alpha, fill = "#7293AD", color = "white", linewidth = 0.3) +
      facet_wrap(. ~ .data$key, scales = "free", labeller = ggplot2::as_labeller(facet.labels)) +
      geom_vline(data = est.df, aes(xintercept = .data$value), linetype = 1, color = "#B43B35", linewidth = 0.8) +
      theme +
      labs(x = "Bootstrap estimate", y = "Frequency", title = title,
           caption = if (show.null && any(grepl("^log\\.", metrics))) "Red: full-sample estimate   |   Dashed: null value" else "Red: full-sample estimate")
    if (show.null && any(grepl("^log\\.", metrics))) {
      null.df <- data.frame(key = factor(metrics[grepl("^log\\.", metrics)], levels = metrics), value = null.value)
      p <- p + geom_vline(data = null.df, aes(xintercept = .data$value),
                          linetype = 2, color = "#303030", linewidth = 0.7)
    }
    p
  }

}

#' Print function for atefit
#'
#' @param x An object of class \code{"atefit"}.
#' @param ... Other parameters
#'
#' @details Display the estimated treatment effects for survival outcomes (log
#' restricted mean time lost ratio and log cumulative hazard ratio) and count outcomes
#' (the log rate ratio).
#'
#' @return No return value
#'
#' @author Thomas Debray
#'
#' @export
#'
print.atefit <- function(x, ...) {

  if (x$response == "binary") {
    cat("Average treatment effect (risk difference):\n\n")
    print(x$risk.difference)
    cat("\nEstimated risks:\n\n")
    print(rbind(x$risk1, x$risk0))
  }
  if (x$response == "count") {
    cat(paste0("Average treatment effect:\n\n"))
    print(x$log.rate.ratio)
    cat(paste0("\nEstimated event rates:\n\n"))
    print(rbind(x$rate1,x$rate0))
  }
  if (x$response == "survival") {
    cat(paste0("Average treatment effect:\n\n"))
    effects <- x$log.rmtl.ratio
    if (!is.null(x$log.hazard.ratio)) effects <- rbind(effects, x$log.hazard.ratio)
    if (!is.null(x$log.cumulative.hazard.ratio)) effects <- rbind(effects, x$log.cumulative.hazard.ratio)
    print(effects)
    cat(paste0("\nEstimated RMST:\n\n"))
    print(rbind(x$rmst1,x$rmst0))
  }

  if (!is.null(x$warning)) {
    cat("\n")
    message(x$warning)
  }
}
