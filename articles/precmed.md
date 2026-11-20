# Getting started with precmed

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/precmed)](https://cran.r-project.org/package=precmed)
[![metacran
downloads](https://cranlogs.r-pkg.org/badges/last-month/precmed)](https://cran.r-project.org/package=precmed)

**precmed** was developed to help researchers with the implementation of
precision medicine in R. A key objective of precision medicine is to
determine the optimal treatment separately for each patient instead of
applying a common treatment to all patients. Personalizing treatment
decisions becomes particularly relevant when treatment response differs
across patients, or when patients have different preferences about
benefits and harms. This package offers statistical methods to develop
and validate prediction models for estimating individualized treatment
effects. These treatment effects are also known as the conditional
average treatment effects (CATEs) and describe how different subgroups
of patients respond to the same treatment. Presently, **precmed**
focuses on the personalization of two competitive treatments using
randomized data from a clinical trial ([Zhao et al.
2013](#ref-zhao2013effectively)) or using real-world data (RWD) from a
non-randomized study ([Yadlowsky et al.
2021](#ref-yadlowsky2020estimation)).

The treatment model estimates propensity scores, and the censoring model
provides probabilities used for inverse probability of censoring
weighting (IPCW). Use `ps_model` for the treatment formula in
[`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md),
[`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md),
and
[`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md).
For survival outcomes, use `censoring_model` for the censoring formula
and `censoring_method` to choose how it is fitted.

For background on marginal and conditional treatment effects, see
Chapter 5 ([Nguyen et al. 2026](#ref-nguyen2026causaleffects)). Chapter
17 introduces CATE prediction and model development ([Hoogland et al.
2026](#ref-hoogland2026introduction)).

## Install precmed

The `precmed` package can be installed from CRAN as follows:

``` r

install.packages("precmed")
```

The latest version can be installed from GitHub as follows:

``` r

install.packages("devtools")
devtools::install_github(repo = "smartdata-analysis-and-statistics/precmed")
```

## Package capabilities

The **precmed** package contains functions to

- Estimate the average treatment effect (ATE) for binary, count and
  survival outcomes using
  [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
- Fit CATE scores for binary, count and survival outcomes using
  [`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
- Develop and cross-validate CATE scores using
  [`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md),
  with count and survival examples
- Plot the proportion of subjects with an estimated treatment effect no
  less than $`c`$ over a range of values for $`c`$([Zhao et al.
  2013](#ref-zhao2013effectively)).
- Compute the area between the average treatment difference curve of
  competing models for the CATE using
  [`abc()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/abc.md)
  ([Zhao et al. 2013](#ref-zhao2013effectively))

For survival outcomes,
[`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
reports the population-average event risk difference at the chosen
follow-up horizon (treatment 1 minus treatment 0), alongside the
restricted mean time lost ratio and cumulative hazard ratio. Survival
CATE scores are on the log restricted mean time lost ratio scale;
subgroup estimates and validation plots can also report risk differences
at the horizon.

## Recommended workflow

Chapter 19 discusses validation strategies and the assessment of
discrimination, calibration, and decision accuracy ([Thao et al.
2026](#ref-thao2026validation)). Chapter 20 provides practical guidance
on validation curves, ABC statistics, and responder diagnostics,
including examples using `precmed` ([Jiang et al.
2026](#ref-jiang2026visualization)).

For count and survival outcomes, the following workflow estimates the
ATE and develops and validates a CATE model:

1.  Use
    [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
    to estimate the overall ATE between two treatment groups and its
    uncertainty using bootstrap resampling.
2.  Compare up to five modelling approaches (e.g., Poisson regression,
    boosting) for estimating the CATE using cross-validation through
    [`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md).
3.  Select the best modelling approach using 3 metrics:
    - Compare the steepness of the validation curves in the validation
      samples across methods using
      [`plot()`](https://rdrr.io/r/graphics/plot.default.html). Two
      side-by-side plots are generated, visualizing the estimated
      average treatment effects in a series of nested subgroups. On the
      left side the curve is shown for the training set, and on the
      right side the curve is shown for the validation set. Each line in
      the plots represents one scoring method (e.g., boosting,
      randomForest) included in the fit.
    - The area between curves (ABC) using
      [`abc()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/abc.md)
      quantifies a model’s ability to capture treatment effect
      heterogeneity. Higher ABC values are preferable as they indicate
      that more treatment effect heterogeneity is captured by the
      scoring method.
    - Compare the distribution of the estimated ATE across different
      levels of the CATE score percentiles using
      [`boxplot()`](https://rdrr.io/r/graphics/boxplot.html).
4.  Fit the chosen approach using
    [`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
    and assess it in an independent dataset where possible.

For binary outcomes,
[`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
and
[`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
are available in the development version. The [binary CATE
example](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Binary-CATE.md)
illustrates held-out subgroup assessment. These binary features are not
in CRAN release 1.1.0.

In the vignettes, we will adopt a different workflow to gradually expose
the user from simple to more complex methods.

### Specify model inputs

When applying
[`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md),
the user has to (at least) input:

- `response`: `binary`, `count`, or `survival`
- `data`: a data frame with individual patient data
- `outcome_model`: a formula describing the outcome and score covariates
  (e.g., outcome ~ age + gender + previous_treatment)
- `outcome_method`: outcome-model fitting methods
- `benefit_method`: score-construction methods paired with
  `outcome_method`
- `ps_model`: a formula describing the treatment model to adjust for
  confounding (e.g., treatment ~ age + previous_treatment)

Use the same model arguments in
[`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)
to validate a score.

## Overview of the vignettes

- [Understanding doubly robust
  methods](https://smartdata-analysis-and-statistics.github.io/precmed/articles/doubly-robust.md)
- [Methods for estimating average treatment
  effects](https://smartdata-analysis-and-statistics.github.io/precmed/articles/ATE-theory.md)
- [Binary ATE
  example](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Binary-examples.md)
- [Count ATE
  example](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Count-examples.md)
- [Survival ATE
  example](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Survival-examples.md)
- [Methods for developing and validating CATE
  scores](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Theoretical-details.md)
- [Binary CATE
  example](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Binary-CATE.md)
- [Count CATE
  example](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Count-CATE.md)
- [Survival CATE
  example](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Survival-CATE.md)
- [Customizing CATE models and
  validation](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Additional-examples.md)

## Abbreviations

| Abbreviation | Full Name                                |
|:-------------|:-----------------------------------------|
| ABC          | Area between curves                      |
| ATE          | Average treatment effect                 |
| CATE         | Conditional average treatment effect     |
| contrastReg  | Contrast regression                      |
| CV           | Cross-validation                         |
| GAM          | Generalized additive model               |
| GBM          | Gradient boosting machine                |
| IPCW         | Inverse probability of censoring weights |
| PM           | Precision medicine                       |
| PS           | Propensity score                         |
| RCT          | Randomized controlled trial              |
| RR           | Rate ratio                               |
| twoReg       | Two regressions                          |

## References

Hoogland, Jeroen, Thomas P. A. Debray, Guillaume Desachy, and Tri-Long
Nguyen. 2026. “Introduction to Methods for Personalizing Medicine.”
Chap. 17 in *Comparative Effectiveness and Personalized Medicine
Research Using Real-World Data*, edited by Thomas P. A. Debray, Tri-Long
Nguyen, and Robert W. Platt. Chapman & Hall/CRC.
<https://doi.org/10.1201/9781003300809-17>.

Jiang, Xiaotong, Gabrielle Simoneau, Fabio Pellegrini, Changyu Shen,
Carl de Moor, and Lu Tian. 2026. “Visualization and Interpretation of
Individualized Treatment Rule Results.” Chap. 20 in *Comparative
Effectiveness and Personalized Medicine Research Using Real-World Data*,
edited by Thomas P. A. Debray, Tri-Long Nguyen, and Robert W. Platt.
Chapman & Hall/CRC. <https://doi.org/10.1201/9781003300809-20>.

Nguyen, Tri-Long, Jessica L. Rohmann, Tat Thang Vo, Robert W. Platt, and
Thomas P. A. Debray. 2026. “Understanding and Defining Causal Effects.”
Chap. 5 in *Comparative Effectiveness and Personalized Medicine Research
Using Real-World Data*, edited by Thomas P. A. Debray, Tri-Long Nguyen,
and Robert W. Platt. Chapman & Hall/CRC.
<https://doi.org/10.1201/9781003300809-5>.

Thao, Le Thi Phuong, Orestis Efthimiou, Tibor V. Varga, Jeroen Hoogland,
and Thomas P. A. Debray. 2026. “Validation of Prediction Models for
Patient Outcomes and Individualized Treatment Effect.” Chap. 19 in
*Comparative Effectiveness and Personalized Medicine Research Using
Real-World Data*, edited by Thomas P. A. Debray, Tri-Long Nguyen, and
Robert W. Platt. Chapman & Hall/CRC.
<https://doi.org/10.1201/9781003300809-19>.

Yadlowsky, Steve, Fabio Pellegrini, Federica Lionetto, Stefan Braune,
and Lu Tian. 2021. “Estimation and Validation of Ratio-Based Conditional
Average Treatment Effects Using Observational Data.” *Journal of the
American Statistical Association* 116 (533): 335–52.
<https://doi.org/10.1080/01621459.2020.1772080>.

Zhao, Lihui, Lu Tian, Tianxi Cai, Brian Claggett, and Lee-Jen Wei. 2013.
“Effectively Selecting a Target Population for a Future Comparative
Study.” *Journal of the American Statistical Association* 108 (502):
527–39. <https://doi.org/10.1080/01621459.2013.770705>.
