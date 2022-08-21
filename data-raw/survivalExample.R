## code to prepare `survivalExample` dataset


library(tidyverse, warn.conflicts = FALSE)
library(MASS)

simsurvdata <- function(n, rateC = 1, seed = NULL){

  #' Generate simulated data as an example for survival data
  #' Assume randomized treatment and independent covariates
  #'
  #' @param n sample size; integer
  #' @param seed randomization seed; integer
  #' @param rateC censoring rate; float
  #' @return A list of two elements:
  #'          data - a data frame of \code{n} rows and `p+4` columns of simulated data; data.frame
  #'            Columns:
  #'                y: observed time to event
  #'                d: event indicator, 1 : event, 0 : censoring
  #'                trt: treatment group
  #'                x.X1, ..., x.Xp: covariate matrix (\code{n} X \code{p})
  #'                yf: maximum follow-up time
  #'          censor.prop - table of censoring proportions of the generated data (overall and by treatment)

  library(truncnorm)
  library(magrittr)
  set.seed(seed)

  if(!is.null(seed)) set.seed(seed)

  # Create an empty shell
  ds <- data.frame(matrix(NA, nrow = n, ncol = 7))
  colnames(ds) <- c("treatment", "ageatindex_centered", "female", "prerelapse_num",
                    "prevDMTefficacy", "premedicalcost", "numSymptoms")

  # Define X, A, and time
  ds %<>%
    mutate(female =              rbinom(n = n, size = 1, prob = 0.75),
           ageatindex_centered = round(rtruncnorm(n, a = 18, b = 64, mean = 48, sd = 12), 0) - 48, # rounded to integers
           prerelapse_num =      rpois(n = n, lambda = 0.44),
           prevDMTefficacy =     sample(x = c("None", "Low efficacy", "Medium and high efficacy"),
                                        size = n, replace = TRUE, prob = c(0.45, 0.44, 0.11)),
           premedicalcost =      pmin(round(exp(rnorm(n = n, mean = 8.9, sd = 1.14)), 2), 600000), # rounded to 2 decimal points
           numSymptoms =         sample(x = c("0", "1", ">=2"),
                                        size = n, replace = TRUE, prob = c(0.67, 0.24, 0.09)), # nuisance variable; do not include in the score
           trt =                 rbinom(n, 1, 1/(1 + exp(-(1 + ageatindex_centered + prerelapse_num + premedicalcost/10000)))), # assign treatment with PS
           treatment =           ifelse(trt == 1, "DMF", "TERI"),
           ttime =               rexp(n) * exp(0.3 * (ageatindex_centered + (prerelapse_num > 0)) * (2 * trt - 1) + prerelapse_num + premedicalcost/10000000) / 3, # event time
           ctime =               pmin(rexp(n) * exp(1 - rateC + premedicalcost/10000), runif(n) * 4) * 1.5, # censoring time
           timetodiscontinuation =  pmin(ttime, ctime) * 350,
           discontinuationindicator =  1 * (ttime < ctime)
  )

  # Show proportion of censoring
  tbl <- table(ds$discontinuationindicator, ds$trt)
  tbl <- cbind(rowSums(tbl), tbl[,2:1])
  tbl <- round(tbl[1,]/colSums(tbl),3)
  names(tbl) <- c('Overall', 'trt = 1', 'trt = 0')

  cat('# Proportion of censoring','\n')
  print(tbl)

  return(list(data = ds, censor.prop = tbl))
}


# Generate data
data <- simsurvdata(n = 4000, rateC = 3.5, seed = 999)$data

# Rename some variables
survivalExample <- data %>% rename(previous_treatment = prevDMTefficacy,
                                   age = ageatindex_centered,
                                   y = timetodiscontinuation,
                                   d = discontinuationindicator,
                                   previous_number_relapses = prerelapse_num,
                                   previous_number_symptoms = numSymptoms,
                                   previous_cost = premedicalcost) %>%
  mutate(previous_treatment = factor(previous_treatment, labels = c("drugA", "drugB", "drugC")),
         previous_number_symptoms = factor(previous_number_symptoms, labels = c("0", "1", ">=2")),
         trt = factor(trt, levels = c(0, 1), labels = c("drug0", "drug1")),
         previous_cost = scale(previous_cost),
         age = age + 48,
         age = scale(age, scale = FALSE)) %>%
  dplyr::select(age, female, previous_treatment, previous_cost, previous_number_symptoms, previous_number_relapses, trt, y, d)

remove(list = c("simsurvdata", "data"))

# save.image("GitHub/PMMS/PrecMed/data/survivalExample.RData")
# setwd("~/Innovation/PMMS/PrecMed/data")
# saveRDS(survivalExample, "survivalExample.rda")

usethis::use_data(survivalExample, overwrite = TRUE)

