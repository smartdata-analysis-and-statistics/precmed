## code to prepare `meanExample` dataset

##TODO: Organize for continuous case.
## Questions: Should I use the same distribution as pre_statusmeasure?
## Questions: Should I add beta as count data example?
library(tidyverse, warn.conflicts = FALSE)
library(MASS)

simmeandata <- function(n, seed = 999,
                        beta =  c(0, 1, 0, 1, 0, 0, 0, 0.0000001),
                        beta.x = c(0, 0, 0, 0, 1, 0, 0, 0),
                        beta.ps = c(0, 1/3, 0, 0, 0, 0, 0, 0)){
  
  
  library(truncnorm)
  library(magrittr)
  library(reshape2)
  library(fastDummies)
  
  set.seed(seed)
  
  # Create an empty shell
  ds <- data.frame(matrix(NA, nrow = n, ncol = 8))
  colnames(ds) <- c("treatment", "ageatindex_centered", "female", "pre_statusmeasure",
                    "prevDMTefficacy", "premedicalcost", "numSymptoms",
                    "group")
  
  # Define X, A, and time
  ds %<>%
    mutate(
      # Covariates
      female =              rbinom(n = n, size = 1, prob = 0.75),
      ageatindex_centered = round(rtruncnorm(n, a = 18, b = 64, mean = 48, sd = 12), 0) - 48, # rounded to integers
      
      pre_statusmeasure =   rnorm(n = n, mean = 0, sd = 4),
      prevDMTefficacy =     sample(x = c("None", "Low_efficacy", "Medium_and_high_efficacy"),
                                   size = n, replace = TRUE, prob = c(0.45, 0.44, 0.11)),
      premedicalcost =      pmin(round(exp(rnorm(n = n, mean = 8.9, sd = 1.14)), 2), 600000), # rounded to 2 decimal points
      numSymptoms =         sample(x = c("0", "1", ">=2"),
                                   size = n, replace = TRUE, prob = c(0.67, 0.24, 0.09)), # nuisance variable; do not include in the score
      group =              "Simulated")
  
  # Treatments
  xmat =               model.matrix(~ ageatindex_centered + female + pre_statusmeasure +  I(pre_statusmeasure > 0) + prevDMTefficacy + premedicalcost, ds) %>% as.matrix()
  betas.ps =           matrix(c(beta.ps[1], # Intercept
                                beta.ps[2], # Age
                                beta.ps[3], # female
                                beta.ps[4], # pre_statusmesure
                                beta.ps[5], # I(pre_statusmeasure >0)
                                beta.ps[6], beta.ps[7], # Medium/high and none DMT efficacy
                                beta.ps[8] # premedicalcost
  ))
  prop =               exp(xmat %*% betas.ps) / (1 + exp(xmat %*% betas.ps))
  trt =                rbinom(n = n, size = 1, prob = prop)
  ds <- ds %>% mutate(treatment = ifelse(trt == 1, "drug1", "drug0"))
  
  # Define Y
  betas =              matrix(c(beta[1], # Intercept
                                beta[2], # Age
                                beta[3], # female
                                beta[4], # pre_statusmesure
                                beta[5], # I(pre_statusmeasure >0)
                                beta[6], beta[7], # Medium/high and none DMT efficacy
                                beta[8] # premedicalcost
  ))
  
  betas.x =            matrix(c(beta.x[1], # Intercept
                                beta.x[2], # Age
                                beta.x[3], # female
                                beta.x[4], # pre_statusmesure
                                beta.x[5], # I(pre_statusmeasure >0)
                                beta.x[6], beta.x[7], # Medium/high and none DMT efficacy
                                beta.x[8] # premedicalcost
  ))
  
  ds <- ds %>% mutate(statusmeasure =  xmat %*% betas + xmat %*% betas.x * (2 * trt - 1) + rnorm(n, 0, 3))
  
  return(list(data = ds))
}


# Generate data
data <- simmeandata(n = 4000, seed = 34)$data

# Rename some variables
meanExample <- data %>% rename(previous_treatment = prevDMTefficacy,
                               age = ageatindex_centered,
                               y = statusmeasure,
                               trt = treatment,
                               previous_status_measure = pre_statusmeasure,
                               previous_number_symptoms = numSymptoms,
                               previous_cost = premedicalcost) %>%
  mutate(previous_treatment = factor(previous_treatment, labels = c("drugA", "drugB", "drugC")),
         previous_number_symptoms = factor(previous_number_symptoms, labels = c("0", "1", ">=2")),
         previous_cost = scale(previous_cost),
         age = age + 48,
         age = scale(age, scale = FALSE)) %>%
  dplyr::select(age, female, previous_treatment, previous_cost, previous_number_symptoms, previous_status_measure, trt, y)

remove(list = c("simmeandata", "data"))

usethis::use_data(meanExample, overwrite = TRUE)
