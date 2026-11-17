# Quick start

This example uses **precmed** to explore whether treatment effects
differ across patients, and to assess a prediction model using
cross-validation. It uses the included simulated multiple sclerosis data
and a count outcome: the number of relapses during follow-up.

## Is precmed suitable for me?

| Question | Guidance |
|----|----|
| What data can I use? | Individual-level data from a randomized trial or an observational study. |
| What treatments can I compare? | Two treatments at a time. |
| What outcomes are supported? | Binary, count and time-to-event outcomes for ATE and CATE fitting in the development version. Continuous outcomes are not yet implemented. |
| What can I estimate? | An overall average treatment effect (ATE) and conditional average treatment effects (CATEs) for patients with given observed characteristics. |
| Can I assess model performance? | [`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md) cross-validation, subgroup curves and ABC are illustrated here for count outcomes and in the survival tutorial. |
| Can I interpret the estimates causally? | Only when the relevant assumptions are justified. Observational analyses require, among other assumptions, adequate control of confounding and sufficient treatment overlap. |

## Install and load

Install the released version from CRAN, then load the package.
Binary-outcome fitting is currently available in the development
version, which can be installed from GitHub:

``` r

install.packages("precmed")
```

``` r

devtools::install_github("smartdata-analysis-and-statistics/precmed")
```

``` r

library(precmed)
data(countExample)
```

The data contain 4,000 patients. `y` is the number of relapses, `years`
is follow-up time, and `trt` identifies the two treatments (`drug0` and
`drug1`). The other variables used below are baseline patient
characteristics. Age is centered and previous cost is centered and
scaled in this dataset.

## Estimate the overall treatment effect

Start with
[`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
to estimate the average treatment effect across the study population.
For count outcomes, it estimates a relapse rate ratio comparing `drug1`
with `drug0`. The outcome model includes a log follow-up offset, and the
treatment model (propensity-score model) describes treatment assignment
using baseline characteristics.

``` r

fit_ate <- atefit(
  response = "count",
  data = countExample,
  outcome_model = y ~ age + female + previous_treatment +
    previous_cost + previous_number_relapses + offset(log(years)),
  ps.model = trt ~ age + previous_treatment,
  outcome_method = "poisson",
  n.boot = 50,
  seed = 999,
  verbose = 0
)
```

The fitted object contains the log rate ratio and bootstrap-based
uncertainty. Exponentiate the estimate and confidence limits to obtain
the rate ratio and its 95% confidence interval; a rate ratio below 1
indicates a lower estimated relapse rate with `drug1`. The confidence
interval uses the bootstrap standard error. The estimator is doubly
robust: under the required causal and other regularity assumptions, it
is consistent if either the outcome model or the treatment model is
correctly specified. The 50 bootstrap replicates used here keep the
example short; use more for a substantive analysis.

``` r

rate_ratio <- exp(coef(fit_ate))
rate_ratio_ci <- exp(confint(fit_ate))
```

| Rate ratio | 95% CI lower | 95% CI upper |
|-----------:|-------------:|-------------:|
|       1.07 |         0.77 |         1.46 |

The overall effect can mask differences across patients. Next, we
explore treatment effect heterogeneity.

## Fit and validate a model

[`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)
fits a conditional average treatment effect (CATE) scoring model and
evaluates it with cross-validation. Here, the outcome model includes a
log follow-up offset to model relapse rates. The treatment model
describes treatment assignment using age and previous treatment. Because
fewer relapses are better, we set `higher.y = FALSE`.

``` r

fit_cv <- catecv(
  response = "count",
  data = countExample,
  score.method = "poisson",
  cate.model = y ~ age + female + previous_treatment +
    previous_cost + previous_number_relapses + offset(log(years)),
  ps.model = trt ~ age + previous_treatment,
  higher.y = FALSE,
  cv.n = 5,
  seed = 999
)
```

## Inspect the validation curves

``` r

plot(fit_cv)
```

![Training and validation curves for estimated treatment effects across
nested patient subgroups.](quick-start_files/figure-html/plot-1.png)

The curves summarize estimated treatment effects in nested subgroups
ranked by the model’s predicted treatment effect. Compare the training
and validation panels: a pattern that persists in validation provides
more support for the model’s ability to distinguish subgroups than one
seen only in training. These results are illustrative, and internal
cross-validation alone does not establish performance in a new
population.

## Next steps

- Read [Understanding doubly robust
  methods](https://smartdata-analysis-and-statistics.github.io/precmed/articles/doubly-robust.md)
  for an introduction to the estimator and its assumptions.
- Follow the [binary ATE
  example](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Binary-examples.md)
  and [binary CATE
  example](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Binary-CATE.md)
  for fixed-horizon response and subgroup validation.
- Follow the [count ATE
  example](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Count-examples.md)
  and [count CATE
  example](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Count-CATE.md)
  to estimate overall effects and assess scoring methods.
- See the [survival ATE
  example](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Survival-examples.md)
  and [survival CATE
  example](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Survival-CATE.md)
  for time-to-event data.
- Explore [model and validation
  options](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Additional-examples.md)
  and the methods for [average treatment
  effects](https://smartdata-analysis-and-statistics.github.io/precmed/articles/ATE-theory.md)
  and [CATE
  scores](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Theoretical-details.md).

For observational data, causal interpretation depends on assumptions
including adequate control of confounding and sufficient treatment
overlap. A CATE describes an average effect conditional on observed
characteristics; it does not reveal an individual patient’s actual
benefit.

For further background, see [Chapter 5: Understanding and defining
causal effects](https://doi.org/10.1201/9781003300809-5) for marginal
and conditional effects, and [Chapter 17: Introduction to methods for
personalizing medicine](https://doi.org/10.1201/9781003300809-17) for
treatment-effect prediction.
