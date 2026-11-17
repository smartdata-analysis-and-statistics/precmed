# Understanding doubly robust methods

How can we estimate a treatment effect when treated and untreated
patients differ at baseline? This introduction compares inverse
probability weighting, outcome-model standardization, and doubly robust
estimation. It then shows how to use `precmed` for a binary outcome. For
the formulas and outcome-specific assumptions, see [Methods for
estimating average treatment
effects](https://smartdata-analysis-and-statistics.github.io/precmed/articles/ATE-theory.md).

## The question and the models

Suppose the outcome is whether a patient has a relapse during one year.
We want the **population-average risk difference**: the proportion who
would relapse if everyone received treatment 1 minus the proportion who
would relapse if everyone received treatment 0. A negative difference
means fewer relapses under treatment 1. This marginal effect is distinct
from a treatment effect conditional on a particular patient’s
covariates.

In observational data, the observed difference in relapse proportions
may reflect differences in baseline characteristics as well as
treatment. Two approaches can help adjust for measured covariates:

- **Inverse probability weighting (IPW)** uses a **treatment model** to
  estimate each patient’s propensity score: the probability of receiving
  treatment 1 given their covariates. Patients receive more weight if
  they received a treatment that was less likely given their covariates.
  The weighted relapse outcomes are then compared between treatment
  groups. Propensity scores can also be used for matching,
  stratification, or regression adjustment; see Chapter 7 for an
  introduction to these methods (Jiang et al. 2026).
- **Outcome-model standardization**, also called **g-computation**, uses
  an **outcome model** to predict relapse given treatment and
  covariates. For each patient in the target population, the model
  predicts the probability of relapse under treatment 0 and again under
  treatment 1, using the patient’s observed covariates. Averaging the
  predictions separately for each treatment gives two estimated
  population-average risks; their difference estimates the risk
  difference. For example, average predicted risks of 20% under
  treatment 0 and 15% under treatment 1 give a risk difference of -5
  percentage points.

Both approaches need causal identification assumptions, including
well-defined treatments and consistency, no unmeasured treatment–outcome
confounding given the adjustment covariates, and positivity: each
treatment must be possible at relevant covariate values. They also need
suitable modeling and regularity assumptions. Outcome-model
standardization generally relies on a correctly specified outcome model,
whereas IPW relies on a correctly specified treatment model. These are
sometimes called **singly robust** approaches with respect to model
misspecification.

## What does doubly robust mean?

A doubly robust estimator combines the outcome and treatment models in a
particular way. For example, augmented inverse probability weighting
(AIPW) starts with outcome-model predictions and adds a
propensity-weighted correction based on the observed prediction errors.
Under the identification and relevant regularity assumptions, this
estimator is consistent if **either** the outcome model **or** the
treatment model is correctly specified; both may be correct (Kurz 2022).
Merely fitting two models does not make an estimator doubly robust.

| Outcome model | Treatment model | AIPW consistency under the other assumptions |
|:--------------|:----------------|:---------------------------------------------|
| Correct       | Correct         | Yes                                          |
| Correct       | Incorrect       | Yes                                          |
| Incorrect     | Correct         | Yes                                          |
| Incorrect     | Incorrect       | Not guaranteed                               |

This is an asymptotic property, not a promise of negligible bias in a
finite sample or valid confidence intervals without further conditions.
Double robustness does not fix unmeasured confounding, poor overlap,
measurement errors, or misspecification of both models. Extreme
propensity scores can make estimates unstable. Truncating propensity
scores may improve stability, but if it changes the true probabilities
it can remove the guarantee based solely on a correct treatment model.
Examine overlap and model plausibility rather than treating double
robustness as a substitute for study design.

## An example with a binary outcome

The following synthetic example compares three estimates of the same
population-average risk difference, treatment 1 minus treatment 0, for
relapse within one year. Treatment assignment depends on age and
previous treatment. The outcome is fully observed. The simulation is
illustrative, not a clinical estimate; its treatment model is correctly
specified by construction, while the generated outcome need not follow
the logistic outcome model exactly. A single simulated dataset cannot
establish which estimator is unbiased.

``` r

library(precmed)
dat <- simulate_binary(n = 4000, seed = 2026, followup = 1, threshold = 0)
table(dat$trt)
#> 
#> drug0 drug1 
#>  1858  2142
```

First, fit a logistic outcome model and average predictions for the same
patients under each treatment. Then fit a logistic treatment model and
calculate IPW means. These are illustrative base R calculations, not
separate `precmed` estimators. We include treatment interactions in the
outcome model and use the same covariates in both approaches below.

``` r

outcome_fit <- glm(y ~ trt * (age + previous_treatment),
                   data = dat, family = binomial())
treatment_fit <- glm(trt ~ age + previous_treatment,
                     data = dat, family = binomial())

under_0 <- under_1 <- dat
under_0$trt <- factor("drug0", levels = levels(dat$trt))
under_1$trt <- factor("drug1", levels = levels(dat$trt))
standardized <- mean(predict(outcome_fit, under_1, type = "response")) -
  mean(predict(outcome_fit, under_0, type = "response"))

# Match the default propensity-score limits used by atefit().
ps <- pmin(pmax(predict(treatment_fit, type = "response"), 0.01), 0.99)
treated <- dat$trt == "drug1"
ipw <- mean(treated * dat$y / ps) - mean((!treated) * dat$y / (1 - ps))
```

Now fit the doubly robust AIPW estimator with
[`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md).
`cate.model` specifies the outcome and covariates;
[`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
adds treatment and, by default, treatment–covariate interactions.
`ps.model` specifies the treatment model. The 50 bootstrap replicates
keep this example quick; use more for a substantive analysis and examine
bootstrap diagnostics.

``` r

fit <- atefit(response = "binary", data = dat,
              cate.model = y ~ age + previous_treatment,
              ps.model = trt ~ age + previous_treatment,
              n.boot = 50, seed = 999)
c(outcome_standardization = standardized, ipw = ipw,
  doubly_robust = unname(fit$risk.difference$estimate))
#> outcome_standardization                     ipw           doubly_robust 
#>             -0.10298071             -0.09978176             -0.10762720
confint(fit)
#>                      2.5 %      97.5 %
#> risk.difference -0.1421899 -0.07306452
```

The estimates need not agree in a finite sample, even when the necessary
assumptions hold. The displayed confidence interval applies to the
[`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
estimate, not the two illustrative calculations. A confidence interval
describes sampling uncertainty under the procedure; it does not account
for unmeasured confounding or establish that either model is correct.

## Choosing and using a method

Doubly robust estimation can be useful when there is uncertainty about
the correct specification of the outcome and treatment models. A singly
robust approach, such as outcome-model standardization or inverse
probability weighting, generally relies on correct specification of its
respective model. A doubly robust estimator such as AIPW is consistent
if **either the outcome model or the treatment model is correctly
specified**, provided the required causal identification and regularity
assumptions also hold. Note that doubly robust methods do not protect
against unmeasured confounding, poor treatment overlap, or
misspecification of both models.

Doubly robust estimation is not automatically more accurate in every
dataset. Specify both models carefully using subject-matter knowledge,
inspect treatment overlap and extreme weights, assess model fit and
sensitivity to modeling choices, and report the target population and
effect scale. Randomization can support exchangeability, but overlap and
other assumptions still matter. For a simple randomized trial, an
intercept-only treatment model may be appropriate; account for the
actual assignment design when choosing it.

In `precmed`,
[`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
estimates a marginal risk difference for binary outcomes, a marginal
rate ratio for count outcomes, and, for survival outcomes at a chosen
horizon, a marginal restricted mean time lost ratio, cumulative hazard
ratio, and event risk difference (treatment 1 minus treatment 0). The
[`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
methods instead construct covariate-conditional treatment-effect scores;
some scoring methods are doubly robust, while others use different
approaches. Score fitting and subsequent validation are distinct tasks,
and an ATE estimate or its bootstrap interval is not validation of a
CATE score. See [Methods for developing and validating CATE
scores](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Theoretical-details.md)
for method-specific details and assumptions.

For survival outcomes, `precmed` also uses a **censoring model** to
supply probabilities for inverse probability of censoring weighting
(IPCW). Double robustness with respect to the outcome and treatment
models does not protect against a misspecified censoring model. The
estimand and support at the chosen survival horizon matter too. For
example, an out-of-range restricted mean survival estimate may trigger
an outcome-model fallback that does not retain the doubly robust
guarantee. See [Methods for estimating average treatment
effects](https://smartdata-analysis-and-statistics.github.io/precmed/articles/ATE-theory.md)
for these qualifications and the survival estimands.

## Other `precmed` vignettes in this series

- [Estimating and validating treatment effects for binary
  outcomes](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Binary-examples.md)
- [Estimating and validating treatment effects for count
  outcomes](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Count-examples.md)
- [Estimating and validating treatment effects for survival
  outcomes](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Survival-examples.md)
- [Customizing CATE models and
  validation](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Additional-examples.md)
- [Methods for estimating average treatment
  effects](https://smartdata-analysis-and-statistics.github.io/precmed/articles/ATE-theory.md)
- [Methods for developing and validating CATE
  scores](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Theoretical-details.md)

## References

Jiang, Tammy, Gabrielle Simoneau, Changyu Shen, Fabio Pellegrini, and
Carl de Moor. 2026. “Confounding Adjustment Using Propensity Score
Methods.” Chap. 7 in *Comparative Effectiveness and Personalized
Medicine Research Using Real-World Data*, edited by Thomas P. A. Debray,
Tri-Long Nguyen, and Robert W. Platt. Chapman & Hall/CRC.
<https://doi.org/10.1201/9781003300809-7>.

Kurz, Christoph F. 2022. “Augmented Inverse Probability Weighting and
the Double Robustness Property.” *Medical Decision Making* 42 (2):
156–67. <https://doi.org/10.1177/0272989X211027181>.
