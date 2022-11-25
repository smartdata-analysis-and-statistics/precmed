#' Simulated data with a continuous outcome
#'
#' A dataset containing a continuous outcome and 6 baseline covariates
#'
#' @docType data
#'
#' @usage data(meanExample)
#'
#' @format A dataframe with 4000 rows (patients) and 9 variables:
#' \describe{
#'   \item{age}{age at baseline, centered to 48 years old, in years}
#'   \item{female}{sex, 0 for male, 1 for female}
#'   \item{previous_treatment}{previous treatment, "drugA", "drugB", or "drugC"}
#'   \item{previous_cost}{previous medical cost, in US dollars}
#'   \item{previous_number_symptoms}{previous number of symptoms, "0", "1", or ">=2"}
#'   \item{previous_number_relapses}{previous number of relapses}
#'   \item{trt}{current treatment, "drug0" or "drug1"}
#'   \item{y}{count outcome, current number of relapses}#'
#' }
#'
#' @keywords datasets
#'
#' @examples
#' data(meanExample)
#' str(meanExample)
"meanExample"

#' Simulated data with count outcome
#'
#' A dataset containing a count outcome, a length of follow-up and 6 baseline covariates
#'
#' @docType data
#'
#' @usage data(countExample)
#'
#' @format A dataframe with 4000 rows (patients) and 9 variables:
#' \describe{
#'   \item{age}{age at baseline, centered to 48 years old, in years}
#'   \item{female}{sex, 0 for male, 1 for female}
#'   \item{previous_treatment}{previous treatment, "drugA", "drugB", or "drugC"}
#'   \item{previous_cost}{previous medical cost, in US dollars}
#'   \item{previous_number_symptoms}{previous number of symptoms, "0", "1", or ">=2"}
#'   \item{previous_number_relapses}{previous number of relapses}
#'   \item{trt}{current treatment, "drug0" or "drug1"}
#'   \item{y}{count outcome, current number of relapses}
#'   \item{years}{length of follow-up, in years}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' data(countExample)
#' str(countExample)
#' rate <- countExample$y / countExample$years
"countExample"

#' Simulated data with survival outcome
#'
#' A dataset containing a time-to-event outcome, an event indicator, treatment group, and 6 baseline covariates
#'
#' @docType data
#'
#' @usage data(survivalExample)
#'
#' @format A dataframe with 4000 rows (patients) and 9 variables:
#' \describe{
#'   \item{age}{age at baseline, centered to 48 years old, in years}
#'   \item{female}{sex, 0 for male, 1 for female}
#'   \item{previous_treatment}{previous treatment, "drugA", "drugB", or "drugC"}
#'   \item{previous_cost}{previous medical cost, in US dollars}
#'   \item{previous_number_symptoms}{previous number of symptoms, "0", "1", or ">=2"}
#'   \item{previous_number_relapses}{previous number of relapses}
#'   \item{trt}{current treatment, "drug0" or "drug1"}
#'   \item{y}{time to first relapse or censoring}
#'   \item{d}{event indicator, 1: relapse, 0: censored}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' data(survivalExample)
#' str(survivalExample)
"survivalExample"

