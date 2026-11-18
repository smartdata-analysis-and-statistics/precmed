# precmed: Estimation and validation of heterogeneous treatment effects

Tools for estimating population-average and covariate-conditional
treatment effects in randomized and observational studies with binary,
count, or time-to-event outcomes. Causal interpretation requires
appropriate identification and modeling assumptions.

## Details

Use
[`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
to estimate population-average treatment effects. Use
[`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
to develop conditional treatment-effect scores and
[`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)
to evaluate scores using cross-validation. Estimating a
population-average effect and validating a conditional treatment-effect
score are distinct tasks. Use
[`diagnostics()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/diagnostics.md)
to inspect available diagnostics and
[`subgroup_estimates()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/subgroup_estimates.md)
to extract subgroup results.

ATE estimation uses augmented inverse probability weighting (AIPW), with
censoring adjustment for survival outcomes. CATE score methods include
T-learners that contrast predictions from separate treatment-specific
outcome models, and doubly robust two-regression and contrast-regression
procedures. The ratio-based doubly robust methods build on Yadlowsky et
al.; the T-learner framework is described by Kunzel et al. See the
outcome-specific help for available methods and their assumptions.

## Outcomes and effect scales

For binary outcomes,
[`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
estimates the population-average risk difference (treatment 1 minus
treatment 0). For count outcomes it estimates a population-average rate
ratio. For survival outcomes at a chosen horizon, it estimates a
restricted mean time lost ratio, cumulative hazard ratio, and event risk
difference (treatment 1 minus treatment 0). The CATE scores describe
effects conditional on covariates; their scales depend on the outcome.
See
[`atefit`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
and
[`catefit`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
for details.

## The question and the models

In observational studies, treatment groups may differ in baseline
characteristics that also affect outcomes. A treatment model estimates
propensity scores, the probabilities of treatment assignment given
covariates. An outcome model predicts the outcome given treatment and
covariates. These models can be used to adjust for measured covariates
when estimating treatment effects. For survival outcomes, a censoring
model supplies probabilities used for inverse probability of censoring
weighting.

## What does doubly robust mean?

[`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
uses doubly robust estimators of population-average treatment effects
for binary, count, and survival outcomes. These estimators combine
outcome and treatment models so that, under the required causal
identification and regularity assumptions, they are consistent if either
the outcome model or the treatment model is correctly specified. Both
models may be correct. Double robustness does not protect against
unmeasured confounding, poor treatment overlap, or misspecification of
both models. The property is asymptotic and does not guarantee
negligible bias in a finite sample. For survival outcomes, it also does
not protect against a misspecified censoring model. If an out-of-range
restricted mean survival estimate triggers an outcome-model fallback,
the substituted estimate does not retain the doubly robust guarantee.
See
[`atefit`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
for details.

## Online guides and vignettes

The package website provides a [quick
start](https://smartdata-analysis-and-statistics.github.io/precmed/articles/quick-start.html),
an introduction to [doubly robust
methods](https://smartdata-analysis-and-statistics.github.io/precmed/articles/doubly-robust.html),
and [all
vignettes](https://smartdata-analysis-and-statistics.github.io/precmed/articles/index.html),
including examples and methodological details for different outcomes.
These resources are available online; the vignettes are not included in
the CRAN package build.

## References

Debray, T. P. A., Nguyen, T.-L., and Platt, R. W. (eds.) (2026).
*Comparative Effectiveness and Personalized Medicine Research Using
Real-World Data*. Chapman & Hall/CRC.
[doi:10.1201/9781003300809](https://doi.org/10.1201/9781003300809) .

Kunzel, S. R., Sekhon, J. S., Bickel, P. J., and Yu, B. (2019).
Metalearners for estimating heterogeneous treatment effects using
machine learning. *Proceedings of the National Academy of Sciences*,
116, 4156–4165.
[doi:10.1073/pnas.1804597116](https://doi.org/10.1073/pnas.1804597116) .

Yadlowsky, S., Pellegrini, F., Lionetto, F., Braune, S., and Tian, L.
(2020). Estimation and validation of ratio-based conditional average
treatment effects using observational data.
[doi:10.1080/01621459.2020.1772080](https://doi.org/10.1080/01621459.2020.1772080)
.

## See also

Useful links:

- <https://github.com/smartdata-analysis-and-statistics/precmed>

- <https://smartdata-analysis-and-statistics.github.io/precmed/>

- Report bugs at
  <https://github.com/smartdata-analysis-and-statistics/precmed/issues>

## Author

**Maintainer**: Thomas Debray <tdebray@fromdatatowisdom.com>
([ORCID](https://orcid.org/0000-0002-1790-2719))

Authors:

- Thomas Debray <tdebray@fromdatatowisdom.com>
  ([ORCID](https://orcid.org/0000-0002-1790-2719))

- Lu Tian <lutian@stanford.edu>
  ([ORCID](https://orcid.org/0000-0002-5893-0169))

- Xiaotong Jiang <xiaotong.phoebe.jiang@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-3698-4526))

- Gabrielle Simoneau <gabrielle.simoneau@mail.mcgill.ca>
  ([ORCID](https://orcid.org/0000-0001-9310-6274))

Other contributors:

- Biogen MA Inc. (Copyright holder of the original implementation
  through version 1.1.0) \[copyright holder\]

- Smart Data Analysis and Statistics B.V. (Copyright holder of its
  contributions after version 1.1.0) \[copyright holder\]

- Stan Wijn ([ORCID](https://orcid.org/0000-0003-3782-6677))
  \[contributor\]

- Joana Caldas \[contributor\]

## Examples

``` r
# A small illustrative example; use more bootstrap replicates in practice.
dat <- simulate_binary(n = 400, seed = 2026)
fit <- atefit(response = "binary", data = dat,
              outcome_model = y ~ age + previous_treatment,
              ps_model = trt ~ age + previous_treatment,
              n.boot = 20, seed = 2026)
coef(fit)
#> risk.difference 
#>      -0.1420822 
diagnostics(fit)
#> [1] method    split     iteration group     subgroup  type      measure  
#> [8] message  
#> <0 rows> (or 0-length row.names)
```
