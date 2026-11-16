#' Simulate a binary relapse outcome from relapse counts
#'
#' Generate synthetic relapse counts over a common follow-up period with
#' \code{simulate_count()}, then define an event when the relapse count exceeds
#' a chosen threshold. The data-generating model and its
#' coefficients are illustrative, not clinical estimates.
#'
#' @param n Positive integer sample size. Default is 4,000.
#' @param seed Integer seed, or \code{NULL} to use the current random-number
#'   state. Supplying a seed sets R's random-number state.
#' @param followup Positive finite common follow-up in years. Default is 1.
#' @param threshold Non-negative integer relapse count above which an event is
#'   recorded. The default, 0, defines an event as at least one relapse.
#' @param ... Additional arguments to \code{simulate_count()}, such as
#'   \code{conditional_rr}, \code{ps_weights}, \code{lp_weights}, or
#'   \code{dispersion}. Do not supply \code{followup_range}; follow-up is
#'   fixed by \code{followup}.
#' @return A data frame with the same covariates, treatment and \code{years} as
#'   \code{simulate_count()}, a 0/1 outcome \code{y}, and the underlying count
#'   \code{relapse_count}. A value of 1 means the count exceeds the threshold.
#'   Treatment labels are \code{drug0} and \code{drug1}.
#' @details The binary outcome is fully observed at the common horizon; there
#'   is no censoring. A treatment effect on the relapse rate generally implies
#'   a different effect on the event probability, depending on follow-up,
#'   the threshold, covariates and count distribution. The model does not
#'   reconstruct any patient population or individual patient records.
#' @examples
#' d <- simulate_binary(n = 500, seed = 1, followup = 1, threshold = 0)
#' table(d$trt, d$y)
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#' @export
simulate_binary <- function(n = 4000L, seed = 20261114L, followup = 1,
                            threshold = 0L, ...) {
  if (!is.numeric(followup) || length(followup) != 1L ||
      !is.finite(followup) || followup <= 0)
    stop("followup must be a positive finite number of years.")
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold < 0 || threshold != floor(threshold))
    stop("threshold must be a non-negative integer.")
  out <- simulate_count(n = n, seed = seed,
                        followup_range = c(followup, followup), ...)
  out$relapse_count <- out$y
  out$y <- as.integer(out$relapse_count > threshold)
  out
}
