#' Simulate an illustrative count-outcome dataset
#'
#' Generate fully synthetic data for demonstrating count treatment-effect
#' methods. Treatment is assigned by a logistic treatment model. Relapse counts
#' follow a Poisson or negative binomial distribution with a log-rate model,
#' treatment interactions and varying follow-up. Coefficients are illustrative,
#' not clinical estimates or a reconstruction of the MarketScan data.
#'
#' @param n Positive integer sample size. Default is 4,000.
#' @param seed Integer seed, or \code{NULL} to use the current random-number
#'   state. Supplying a seed sets R's random-number state.
#' @param ps_weights Named coefficients for treatment assignment. Recognized
#'   names are \code{intercept}, \code{age}, \code{drugB}, and \code{drugC};
#'   named subsets override the defaults. Age is measured in decades from 48
#'   years, and drugA is the reference previous treatment.
#' @param lp_weights Named coefficients for the log relapse rate. Recognized
#'   names are \code{intercept}, \code{age}, \code{female},
#'   \code{previous_cost}, \code{previous_number_relapses},
#'   \code{treatment_age}, and \code{treatment_relapses}; named subsets
#'   override the defaults. Age is measured in decades from 48 years and cost
#'   is centered at 13,824 and divided by 20,172 US dollars. The intercept is
#'   the log annual rate under drug0 at the reference covariate values.
#' @param conditional_rr Positive rate ratio for drug1 versus drug0 at age 48
#'   and zero previous relapses. Elsewhere its log is \code{log(conditional_rr)
#'   + treatment_age * (age_in_years - 48)/10 + treatment_relapses *
#'   previous_number_relapses}.
#' @param followup_range Two finite positive values specifying the minimum and
#'   maximum follow-up in years. Follow-up is sampled uniformly and independently
#'   of covariates, treatment and the outcome process; equal values give constant
#'   follow-up.
#' @param dispersion Positive negative-binomial size parameter, or \code{Inf}
#'   (the default) for Poisson counts. For finite size, the count variance is
#'   \code{mean + mean^2 / dispersion}.
#' @return A data frame with centered \code{age}, \code{female},
#'   \code{previous_treatment}, scaled \code{previous_cost},
#'   \code{previous_number_symptoms}, \code{previous_number_relapses},
#'   \code{trt}, count outcome \code{y}, and follow-up in years \code{years}.
#'   Age and cost retain centering and scaling attributes. Treatment labels are
#'   \code{drug0} and \code{drug1}.
#' @details The model assumes a constant individual event rate over follow-up;
#'   conditional on covariates, treatment and follow-up, expected counts equal
#'   the annual rate multiplied by follow-up. Extreme coefficients may cause
#'   poor overlap or invalid rates. The generator is for illustration, not a
#'   calibrated representation of any patient population.
#' @examples
#' toy <- simulate_count(n = 500, seed = 1)
#' table(toy$trt)
#' toy_overdispersed <- simulate_count(n = 500, seed = 1,
#'   conditional_rr = 0.8, dispersion = 2)
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#' @export
simulate_count <- function(
    n = 4000L, seed = 20261114L,
    ps_weights = c(intercept = 0.11, age = 1.15, drugB = 0.25, drugC = -0.20),
    lp_weights = c(intercept = log(0.5), age = 0.25, female = 0.12,
                   previous_cost = 0.10, previous_number_relapses = 0.16,
                   treatment_age = -0.20, treatment_relapses = -0.12),
    conditional_rr = 0.9, followup_range = c(0.5, 3), dispersion = Inf) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 ||
      n != floor(n) || n > .Machine$integer.max)
    stop("n must be a positive integer.")
  merge_weights <- function(x, defaults, label) {
    if (!is.numeric(x) || is.null(names(x)) || anyDuplicated(names(x)) ||
        any(!nzchar(names(x))) || any(!names(x) %in% names(defaults)) || any(!is.finite(x)))
      stop(label, " must contain finite, uniquely named recognized coefficients.")
    defaults[names(x)] <- x
    defaults
  }
  ps_weights <- merge_weights(ps_weights,
    c(intercept = 0.11, age = 1.15, drugB = 0.25, drugC = -0.20), "ps_weights")
  lp_weights <- merge_weights(lp_weights,
    c(intercept = log(0.5), age = 0.25, female = 0.12, previous_cost = 0.10,
      previous_number_relapses = 0.16, treatment_age = -0.20, treatment_relapses = -0.12), "lp_weights")
  if (!is.numeric(conditional_rr) || length(conditional_rr) != 1L ||
      !is.finite(conditional_rr) || conditional_rr <= 0)
    stop("conditional_rr must be a finite positive number.")
  if (!is.numeric(followup_range) || length(followup_range) != 2L ||
      any(!is.finite(followup_range)) || any(followup_range <= 0) || followup_range[1] > followup_range[2])
    stop("followup_range must contain two positive finite values in increasing order.")
  if (!is.numeric(dispersion) || length(dispersion) != 1L || is.na(dispersion) || dispersion <= 0)
    stop("dispersion must be positive or Inf.")
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed != floor(seed) || abs(seed) > .Machine$integer.max)
      stop("seed must be an integer or NULL.")
    set.seed(seed)
  }
  age <- pmin(65, pmax(18, stats::rnorm(n, 48, 9))) - 48
  female <- stats::rbinom(n, 1, 0.72)
  previous <- sample(c("drugA", "drugB", "drugC"), n, TRUE, prob = c(0.45, 0.35, 0.20))
  cost <- (stats::rgamma(n, shape = 2, scale = 6912) - 13824) / 20172
  symptoms <- sample(c("0", "1", ">=2"), n, TRUE, prob = c(0.45, 0.35, 0.20))
  relapses <- pmin(stats::rpois(n, 0.8), 5L)
  p <- stats::plogis(ps_weights[["intercept"]] + ps_weights[["age"]] * age / 10 +
    ps_weights[["drugB"]] * (previous == "drugB") + ps_weights[["drugC"]] * (previous == "drugC"))
  treatment <- stats::rbinom(n, 1, p)
  years <- stats::runif(n, followup_range[1], followup_range[2])
  lp <- lp_weights[["intercept"]] + lp_weights[["age"]] * age / 10 +
    lp_weights[["female"]] * female + lp_weights[["previous_cost"]] * cost +
    lp_weights[["previous_number_relapses"]] * relapses + treatment *
    (log(conditional_rr) + lp_weights[["treatment_age"]] * age / 10 +
       lp_weights[["treatment_relapses"]] * relapses)
  means <- exp(lp) * years
  if (any(!is.finite(means)) || any(means <= 0))
    stop("Outcome coefficients produced invalid expected counts; use less extreme values.")
  y <- if (is.infinite(dispersion)) stats::rpois(n, means) else stats::rnbinom(n, mu = means, size = dispersion)
  if (any(!is.finite(y))) stop("Outcome coefficients produced invalid counts; use less extreme values.")
  out <- data.frame(age = age, female = female,
    previous_treatment = factor(previous, levels = c("drugA", "drugB", "drugC")),
    previous_cost = cost, previous_number_symptoms = factor(symptoms, levels = c("0", "1", ">=2")),
    previous_number_relapses = relapses,
    trt = factor(ifelse(treatment == 1, "drug1", "drug0"), levels = c("drug0", "drug1")),
    y = y, years = years)
  attr(out$age, "scaled:center") <- 48
  attr(out$previous_cost, "scaled:center") <- 13824
  attr(out$previous_cost, "scaled:scale") <- 20172
  out
}
