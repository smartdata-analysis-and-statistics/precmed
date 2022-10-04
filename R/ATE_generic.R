#' Histogram of bootstrap estimates
#'
#' @param x An object of class \code{"atefit"}.
#' @param bins Number of bins.
#' @param alpha Opacity
#' @param ... Other parameters
#'
#' Create a histogram displaying the distribution of the bootstrap estimates.
#' The red vertical reference line represents the final estimate.
#'
#' @importFrom dplyr %>% add_row mutate
#' @importFrom ggplot2 ggplot geom_histogram geom_vline
#' @importFrom tidyr gather
#' @export
#'
plot.atefit <- function(x, bins = 50, alpha = 0.7, ...) {

  if (x$response == "count") {
    n.boot <- length(x$trt.boot)

    x$trt.boot %>% as.data.frame() %>%
      ggplot(aes(x = .data$.)) +
      geom_histogram(bins = bins, alpha = alpha) +
      theme_classic() +
      geom_vline(xintercept = x$log.rate.ratio$estimate, linetype = 1, color = "red") +
      labs(x = "Bootstrap values of the log rate ratio", y = "Frequency",
           title = paste0(n.boot, " bootstrap iterations"))
  } else if (x$response == "survival") {
    n.boot <- nrow(x$trt.boot)

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
      theme_classic() +
      labs(x = "Bootstrap values", y = "Frequency", title = paste0(n.boot, " bootstrap iterations"))
  }

}

#' Print function for atefit
#'
#' @param x An object of class \code{"atefit"}.
#' @param ... Other parameters
#'
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
