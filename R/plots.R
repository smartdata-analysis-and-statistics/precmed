# ------------------------------------------------------------------
#
# Project: Precision Medicine MS (precmed) - Comprehensive R package
#
# Purpose: Plot functions for precmed package
#
# Platform: Windows
# R Version: 4.1.0
#



#' Two side-by-side line plots of validation curves from the \code{"precmed"} object
#'
#' Provides validation curves in two side-by-side plots, visualizing the estimated ATEs in a series
#' of nested subgroups in the training set and validation set separately, where each line represents
#' one scoring method specified in \code{\link{catecv}()} or \code{\link{catecvmean}()}. This should be run
#' only after results of \code{\link{catecv}()} or \code{\link{catecvmean}()} have been obtained.
#'
#' @param x An object of class \code{"precmed"}.
#' @param cv.i A positive integer indicating the index of the CV iteration results to be plotted.
#' Allowed values are: a positive integer \eqn{<=} \code{cv.n} in \code{\link{catecv}()} or
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
#' @details \code{\link{plot}()} takes in outputs from \code{\link{catecv}()} and generates two plots
#' of validation curves side-by-side, one for the training set and one for validation set.
#' Separate validation curves are produced for each scoring method specified via \code{score.method}
#' in \code{\link{catecv}()} or \code{\link{catecvmean}()}.
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
#' @seealso \code{\link{abc}()} and \code{\link{boxplot}()} for \code{"precmed"} objects.
#'
#' @examples
#' \dontrun{
#' # Count outcome
#' cv_count <- catecv(response = "count",
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
#' tau0 <- with(survivalExample,
#'              min(quantile(y[trt == "drug1"], 0.95), quantile(y[trt == "drug0"], 0.95)))
#' cv_surv <- catecv(response = "survival",
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
#'
#' # default setting, plot RMTL ratios in both training and validation sets
#' plot(x = cv_surv)
#'
#' # plot hazard ratio
#' plot(x = cv_surv, plot.hr = TRUE)
#'
#' # Continuous outcome
#' cv_mean <- catecvmean(cate.model = y ~ age +
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
#' plot(x = cv_mean)
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

plot.precmed <- function(x,
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
  if (plot.hr == TRUE & x$response != "survival") stop("Hazard ratio plot is only available for precmed object with x$response equals to 'survival'.")
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
      plot.ratio <- "Rate ratio"
    } else if (x$response == "survival"){
      if (plot.hr == TRUE){
        plot.ratio <- "Hazard ratio"
      } else {
        plot.ratio <- "RMTL ratio"
      }
    } else if (x$response == "continuous"){
      plot.ratio <- "Mean difference"
    }

    if (is.null(x$formulas$trt_labels) == TRUE & x$response == "continuous"){
      ylab <- paste0(plot.ratio, " ", trt, "=1 - ", trt, "=0\nin each subgroup")

    } else if (is.null(x$formulas$trt_labels) == TRUE & x$response != "continuous") {
      ylab <- paste0(plot.ratio, " ", trt, "=1 - ", trt, "=0\nin each subgroup")
    } else {
      ylab <- paste0(plot.ratio, " ", trt, "=", x$formulas$trt_labels[2], " - ", trt, "=", x$formulas$trt_labels[1], "\nin each subgroup")
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


#' A set of box plots of estimated ATEs from the \code{"precmed"} object
#'
#' Provides box plots which depict distributions of estimated ATEs for each multi-category subgroup in
#' the validation set across all cross-validation iterations. The subgroups are mutually exclusive and
#' are categorized by the CATE score percentiles (\code{prop.multi} specified in \code{\link{catecv}()} or
#' \code{\link{catecvmean}()}). Box plots of mutually exclusive subgroups are constructed separately by scoring
#' method specified in \code{\link{catecv}()}. This should be run only after results of \code{\link{catecv}()} or
#' \code{\link{catecvmean}()}) have been obtained.
#'
#' @param x An object of class \code{"precmed"}.
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
#' @details \code{\link{boxplot}()} takes in outputs from \code{\link{catecv}()} and generates
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
#' @seealso \code{\link{plot}} and \code{\link{abc}()} for \code{"precmed"} objects.
#'
#' @examples
#' \dontrun{
#' # Count outcome
#' cv_count <- catecv(response = "count",
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
#' tau0 <- with(survivalExample,
#'              min(quantile(y[trt == "drug1"], 0.95), quantile(y[trt == "drug0"], 0.95)))
#' cv_surv <- catecv(response = "survival",
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
#' boxplot(x = cv_surv, ylab = "RMTL ratio of drug1 vs drug0 in each subgroup")
#'
#'# Continuous outcome
#' cv_mean <- catecvmean(cate.model = y ~ age +
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
#' boxplot(x = cv_mean)
#'}
#'
#' @export
#'
#' @importFrom graphics boxplot
#' @importFrom dplyr mutate_at vars contains
#' @importFrom ggplot2 ggplot aes facet_wrap geom_boxplot geom_hline geom_line labs scale_color_manual scale_linetype_manual scale_y_continuous theme theme_classic
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data


boxplot.precmed <- function(x, ylab = NULL,
                            plot.hr = FALSE,
                            grayscale = FALSE, ...) {

  # Must be count or survival data for this function
  stopifnot(x$response == "count" | x$response == "survival" | x$response == "continuous")

  # Argument checks on plot.hr
  if (!(plot.hr %in% c(TRUE, FALSE))) stop("plot.hr has to be boolean.")
  if (plot.hr == TRUE & x$response != "survival") stop("Hazard ratio plot is only available for precmed objects with x$response is \"survival\".")

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
      plot.ratio = "Rate ratio"
    } else if (x$response == "survival"){
      if (plot.hr == TRUE){
        plot.ratio = "Hazard ratio"
      } else {
        plot.ratio = "RMTL ratio"
      }
    } else if (x$response == "continuous"){
      plot.ratio = "Mean difference"
    }

    if (is.null(x$formulas$trt_labels) == TRUE & x$response == "continuous"){
      ylab <- paste0(plot.ratio, " ", trt, "=1 - ", trt, "=0\nin each subgroup")

    } else if (is.null(x$formulas$trt_labels) == TRUE & x$response != "continuous") {
      ylab <- paste0(plot.ratio, " ", trt, "=1 - ", trt, "=0\nin each subgroup")
    } else {
      ylab <- paste0(plot.ratio, " ", trt, "=", x$formulas$trt_labels[2], " - ", trt, "=", x$formulas$trt_labels[1], "\nin each subgroup")
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
