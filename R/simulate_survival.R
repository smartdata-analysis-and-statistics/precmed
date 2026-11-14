#' Simulate an illustrative survival dataset
#'
#' Generate fully synthetic data for demonstrating survival treatment-effect
#' methods. Baseline covariates are generated independently of patient records.
#' Treatment is assigned using a logistic propensity model; the numbers in
#' each treatment group arise from those probabilities and vary by sample.
#' Event times follow a Weibull proportional hazards model with treatment
#' interactions with age and previous relapses. Censoring is independent of covariates and
#' event times. These data are illustrative and are not clinical estimates.
#'
#' @param n Positive integer giving the total sample size. Default is 4,000.
#'   Treatment group sizes are determined by the propensity model.
#' @param seed An integer seed for reproducibility, or \code{NULL} to use the
#'   current random-number state. Supplying a seed sets R's random-number state.
#' @param ps_weights Named numeric coefficients for the logistic propensity model.
#'   Recognized names are \code{intercept}, \code{age},
#'   \code{drugB}, and \code{drugC}; named subsets override the defaults.
#'   Age is measured in decades from 48 years;
#'   previous-treatment drugA is the reference category. The default intercept
#'   yields approximately 54% assigned drug1, similar to the aggregate treatment
#'   split in the package's original survival example.
#' @param lp_weights Named numeric coefficients for the log event hazard.
#'   Recognized names are \code{intercept}, \code{age}, \code{female},
#'   \code{previous_cost}, \code{previous_number_relapses}, and
#'   \code{treatment_age}, and \code{treatment_relapses}; named subsets override
#'   the defaults. Age is measured in decades from 48 years. The interaction
#'   coefficients describe changes in the log HR per decade and per previous
#'   relapse, respectively. Previous cost is centered at 14,362 and divided by
#'   24,266 dollars.
#' @param conditional_hr Positive conditional instantaneous hazard ratio for
#'   drug1 versus drug0 at age 48 with zero previous relapses. For other values,
#'   the log HR is \code{log(conditional_hr) + treatment_age *
#'   (age_in_years - 48)/10 + treatment_relapses * previous_number_relapses},
#'   where \code{age_in_years} is the uncentered age.
#'
#' @return A data frame with \code{n} rows and nine columns:
#'   centered \code{age}, \code{female}, \code{previous_treatment}, scaled
#'   \code{previous_cost}, \code{previous_number_symptoms},
#'   \code{previous_number_relapses}, \code{trt}, observed time \code{y}, and
#'   event indicator \code{d}. The age and cost columns retain centering and
#'   scaling attributes. Times are in days; 1 indicates relapse and 0 censoring.
#' @details Event times use a Weibull baseline with shape 1.15 and scale 650
#'   days. Censoring times are exponential with mean 450 days and administrative
#'   censoring at 900 days. Extreme propensity coefficients can produce poor
#'   treatment overlap or a sample with no patients in one treatment group.
#' @examples
#' toy <- simulate_survival(n = 200, seed = 1)
#' table(toy$trt)
#' toy_unequal <- simulate_survival(n = 200, seed = 1,
#'   ps_weights = c(intercept = -1), conditional_hr = 0.8)
#' table(toy_unequal$trt)
#' @export
simulate_survival <- function(
    n = 4000L, seed = 20261114L,
    ps_weights = c(intercept = 0.11, age = 0.25, drugB = 0.25, drugC = -0.20),
    lp_weights = c(intercept = 0, age = 0.20, female = 0.12, previous_cost = 0.10,
                   previous_number_relapses = 0.16, treatment_age = -0.28, treatment_relapses = -0.16),
    conditional_hr = exp(-0.35)) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 ||
      n != floor(n) || n > .Machine$integer.max)
    stop("n must be a positive integer.")
  merge_weights <- function(x, defaults, label) {
    if (!is.numeric(x) || is.null(names(x)) || anyDuplicated(names(x)) ||
        any(!names(x) %in% names(defaults)) || any(!is.finite(x)))
      stop(label, " must contain finite, uniquely named recognized coefficients.")
    defaults[names(x)] <- x
    defaults
  }
  ps_weights <- merge_weights(ps_weights,
    c(intercept = 0.11, age = 0.25, drugB = 0.25, drugC = -0.20), "ps_weights")
  lp_weights <- merge_weights(lp_weights,
    c(intercept = 0, age = 0.20, female = 0.12, previous_cost = 0.10,
      previous_number_relapses = 0.16, treatment_age = -0.28, treatment_relapses = -0.16), "lp_weights")
  if (!is.numeric(conditional_hr) || length(conditional_hr) != 1L ||
      !is.finite(conditional_hr) || conditional_hr <= 0)
    stop("conditional_hr must be a finite positive number.")
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed != floor(seed) || abs(seed) > .Machine$integer.max)
      stop("seed must be an integer or NULL.")
    set.seed(seed)
  }

  age_years <- pmin(65, pmax(18, stats::rnorm(n, 48, 9)))
  female <- stats::rbinom(n, 1, 0.72)
  prev <- sample(c("drugA", "drugB", "drugC"), n, TRUE, prob = c(0.45, 0.35, 0.20))
  cost <- stats::rgamma(n, shape = 2, scale = 7181)
  symptoms <- sample(c("0", "1", ">=2"), n, TRUE, prob = c(0.45, 0.35, 0.20))
  relapses <- pmin(stats::rpois(n, 0.8), 5L)
  p <- stats::plogis(ps_weights[["intercept"]] + ps_weights[["age"]] * (age_years - 48) / 10 +
                      ps_weights[["drugB"]] * (prev == "drugB") + ps_weights[["drugC"]] * (prev == "drugC"))
  treatment <- stats::rbinom(n, 1, p)
  dat <- data.frame(age_years, female, prev, cost, symptoms, relapses, treatment)
  age <- dat$age_years - 48
  cost <- (dat$cost - 14362) / 24266
  lp <- lp_weights[["intercept"]] + lp_weights[["age"]] * age / 10 + lp_weights[["female"]] * dat$female +
    lp_weights[["previous_cost"]] * cost + lp_weights[["previous_number_relapses"]] * dat$relapses +
    dat$treatment * (log(conditional_hr) + lp_weights[["treatment_age"]] * age / 10 +
                     lp_weights[["treatment_relapses"]] * dat$relapses)
  event_time <- 650 * (-log(stats::runif(nrow(dat))) / exp(lp))^(1 / 1.15)
  if (any(!is.finite(event_time)) || any(event_time <= 0))
    stop("Outcome coefficients produced invalid event times; use less extreme values.")
  censor_time <- pmin(stats::rexp(nrow(dat), rate = 1 / 450), 900)
  out <- data.frame(age = age, female = dat$female,
                    previous_treatment = factor(dat$prev, levels = c("drugA", "drugB", "drugC")),
                    previous_cost = cost,
                    previous_number_symptoms = factor(dat$symptoms, levels = c("0", "1", ">=2")),
                    previous_number_relapses = dat$relapses,
                    trt = factor(ifelse(dat$treatment == 1, "drug1", "drug0"), levels = c("drug0", "drug1")),
                    y = pmin(event_time, censor_time), d = as.integer(event_time <= censor_time))
  attr(out$age, "scaled:center") <- 48
  attr(out$previous_cost, "scaled:center") <- 14362
  attr(out$previous_cost, "scaled:scale") <- 24266
  out
}
