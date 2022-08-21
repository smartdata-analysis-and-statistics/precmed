## code to prepare `countExample` dataset


library(tidyverse, warn.conflicts = FALSE)
library(MASS)

simcountdata <- function(n, seed = 999,
                    beta = c(-0.5, -0.25, 0, 0.25, 0.5),
                    beta.x = c(-1.54, -0.01, 0.06, 0.25, 0.5, 0.13, 0.0000003),
                    percentiles = seq(0, 1, by = 0.2)){

  #' Generate simulated count data with settings based on real-world data
  #' Assume randomized treatment and independent covariates
  #'
  #' @param n sample size; integer
  #' @param seed randomization seed; integer
  #' @param beta coefficients characterizing treatment effect heterogeneity; vector of length 5
  #'             beta[1]*trt*I(high responder to DMF) +
  #'             beta[2]*trt*I(moderate responder to DMF) +
  #'             beta[3]*trt*I(neutral) +
  #'             beta[4]*trt*I(moderate responder to TERI) +
  #'             beta[5]*trt*I(high responder to TERI)
  #'             In the absence of treatment effect heterogeneity, set all beta[1:5] with the same values
  #' @param beta.x coefficients for main effects of other covariates in the rate; vector of length 7
  #'               beta.x[1] (intercept)
  #'               beta.x[2]*ageatindex_centered
  #'               beta.x[3]*female
  #'               beta.x[4]*prerelapse_num
  #'               beta.x[5]*prevDMTefficacy== mediumhigh (reference low efficacy)
  #'               beta.x[6]*prevDMTefficacy== none (reference low efficacy)
  #'               beta.x[7]*premedicalcost
  #' @param percentiles percentiles (cumulative, monotone increasing) to define each responder subgroup; vector of floats of length 6
  #' @return data - a data frame of \code{n} rows and 8 columns of simulated data; data.frame

  library(truncnorm)
  library(magrittr)
  library(reshape2)
  library(fastDummies)

  set.seed(seed)

  if (percentiles[1] != 0 | percentiles[length(percentiles)] != 1 | length(percentiles) != 6){message("Wrong values of percentiles!")}

  # Create an empty shell
  ds <- data.frame(matrix(NA, nrow = n, ncol = 10))
  colnames(ds) <- c("treatment", "ageatindex_centered", "female", "prerelapse_num",
                    "prevDMTefficacy", "premedicalcost", "numSymptoms",
                    "postrelapse_num", "finalpostdayscount", "group")

  # Define X, A, and time
  ds %<>%
    mutate(trt =                 rbinom(n = n, size = 1, prob = 0.75),
           treatment =           ifelse(trt == 1, "DMF", "TERI"),
           female =              rbinom(n = n, size = 1, prob = 0.75),
           ageatindex_centered = round(rtruncnorm(n, a = 18, b = 64, mean = 48, sd = 12), 0) - 48, # rounded to integers
           prerelapse_num =      rpois(n = n, lambda = 0.44),
           prevDMTefficacy =     sample(x = c("None", "Low efficacy", "Medium and high efficacy"),
                                        size = n, replace = TRUE, prob = c(0.45, 0.44, 0.11)),
           premedicalcost =      pmin(round(exp(rnorm(n = n, mean = 8.9, sd = 1.14)), 2), 600000), # rounded to 2 decimal points
           numSymptoms =         sample(x = c("0", "1", ">=2"),
                                        size = n, replace = TRUE, prob = c(0.67, 0.24, 0.09)), # nuisance variable; do not include in the score
           finalpostdayscount =  ceiling(rgamma(n = n, shape = 0.9, scale = 500)), # rounded up to integers
           finalpostdayscount =  ifelse(finalpostdayscount > 2096, 2096, finalpostdayscount), # truncate at the max follow up day, 2096
           finalpostdayscount =  ifelse((finalpostdayscount > 2090) & (runif(1, 0, 1) < .5), 29, finalpostdayscount), # mimic the 1 month peak;  move roughly half of the large values to 29
           group =              "Simulated")

  # Define Y
  xmat.score <- model.matrix(~ ageatindex_centered + female + prerelapse_num + prevDMTefficacy + premedicalcost, ds) %>% as.matrix()
  gamma <- matrix(c(-0.33, # Intercept
                    -0.001, # Age
                    0.05, # female
                    -0.002, # prerelapse_num
                    0.33, 0.02, # Medium/high and none DMT efficacy
                    -0.0000005), nrow = 7) # premedicalcost
  ds <- ds %>% mutate(score = exp(xmat.score %*% gamma),
                      Iscore = cut(score, quantile(score, percentiles), include.lowest = T, labels = seq(1, 5)))

  xmat.rate <- model.matrix(~ ageatindex_centered + female + prerelapse_num + prevDMTefficacy + premedicalcost +
                              Iscore + trt + trt*Iscore, ds) %>% as.matrix()

  betas <- matrix(c(beta.x[1], # Intercept
                    beta.x[2], # Age
                    beta.x[3], # female
                    beta.x[4], # prerelapse_num
                    beta.x[5], beta.x[6], # Medium/high and none DMT efficacy
                    beta.x[7], # premedicalcost
                    # Previous beta.x: -1.54 Intercept, -0.01 Age, 0.06 female, 0.47 prerelapse_num, 0.68, 0.13 Medium/high and none DMT efficacy, 0.000003 premedicalcost
                    0, 0, 0, 0, # Iscore categories 2:5 (reference group is Iscore1, high responders to DMF)
                    beta[1], beta[2] - beta[1], beta[3] - beta[1], beta[4] - beta[1], beta[5] - beta[1]))
  rate <-  exp(xmat.rate %*% betas)
  ds <- ds %>% mutate(postrelapse_num = rnegbin(n = n, mu = rate*finalpostdayscount/365.25, theta = 3))

  return(list(data = ds, betas = betas, percentiles = percentiles))
}


# Generate data
data <- simcountdata(n = 4000, seed = 34,
                      beta = c(log(0.4), log(0.5), log(1), log(1.1), log(1.2)), # moderate level of heterogeneity
                      beta.x = c(-1.54, -0.01, 0.06, 0.25, 0.5, 0.13, 0.0000003))$data

# Rename some variables
countExample <- data %>% rename(previous_treatment = prevDMTefficacy,
                                age = ageatindex_centered,
                                y = postrelapse_num,
                                years = finalpostdayscount,
                                previous_number_relapses = prerelapse_num,
                                previous_number_symptoms = numSymptoms,
                                previous_cost = premedicalcost) %>%
  mutate(previous_treatment = factor(previous_treatment, labels = c("drugA", "drugB", "drugC")),
         previous_number_symptoms = factor(previous_number_symptoms, labels = c("0", "1", ">=2")),
         years = years / 365.25,
         trt = factor(trt, levels = c(0, 1), labels = c("drug0", "drug1")),
         previous_cost = scale(previous_cost),
         age = age + 48,
         age = scale(age, scale = FALSE)) %>%
  dplyr::select(age, female, previous_treatment, previous_cost, previous_number_symptoms, previous_number_relapses, trt, y, years)

remove(list = c("simcountdata", "data"))

# save.image("GitHub/PMMS/PrecMed/data/countExample.RData")
usethis::use_data(countExample, overwrite = TRUE)
