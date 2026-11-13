# A general introduction

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/precmed)](https://cran.r-project.org/package=precmed)
[![metacran
downloads](https://cranlogs.r-pkg.org/badges/last-month/precmed)](https://cran.r-project.org/package=precmed)

# precmed: Precision Medicine in R

A doubly robust precision medicine approach to estimate and validate
conditional average treatment effects.

## General introduction

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
2020](#ref-yadlowsky2020estimation)).

### How to install?

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

- Estimate the average treatment effect (ATE) for count and survival
  outcome data using
  [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
- Estimate the conditional average treatment effect (CATE) using
  [`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
- Development and Cross-validation of CATE using
  [`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md)
- Plot the proportion of subjects with an estimated treatment effect no
  less than $`c`$ over a range of values for $`c`$([Zhao et al.
  2013](#ref-zhao2013effectively)).
- Compute the area between the average treatment difference curve of
  competing models for the CATE using
  [`abc()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/abc.md)
  ([Zhao et al. 2013](#ref-zhao2013effectively))

## Recommended workflow

We recommend the following workflow to develop a model for estimating
the CATE in order to identify treatment effect heterogeneity:

1.  Compare up to five modelling approaches (e.g., Poisson regression,
    boosting) for estimating the CATE using cross-validation through
    [`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md).
2.  Select the best modelling approach using 3 metrics:
    - Compare the steepness of the validation curves in the validation
      samples across methods using
      [`plot()`](https://rdrr.io/r/graphics/plot.default.html). Two
      side-by-side plots are generated, visualizing the estimated
      average treatment effects in a series of nested subgroups. On the
      left side the curve is shown for the training set, and on the
      right side the curve is shown for the validation set. Each line in
      the plots represents one scoring method (e.g., boosting,
      randomForest) specified under the argument `score.method`.
    - The area between curves (ABC) using
      [`abc()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/abc.md)
      quantifies a model’s ability to capture treatment effect
      heterogeneity. Higher ABC values are preferable as they indicate
      that more treatment effect heterogeneity is captured by the
      scoring method.
    - Compare the distribution of the estimated ATE across different
      levels of the CATE score percentiles using
      [`boxplot()`](https://rdrr.io/r/graphics/boxplot.html).
3.  Apply the best modelling approach in the original data or in a new
    external dataset using
    [`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md).
4.  Optional. Use
    [`atefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.md)
    to estimate ATE between 2 treatment groups with a doubly robust
    estimator and estimate the variability of the ATE with a bootstrap
    approach.

In the vignettes, we will adopt a different workflow to gradually expose
the user from simple to more complex methods.

### User input

When applying
[`catefit()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.md)
or
[`catecv()`](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.md),
the user has to (at least) input:

- `response`: type of outcome/response (either `count` or `survival`)  
- `data`: a data frame with individual patient data  
- `score.method`: methods to estimate the CATE (e.g., `boosting`,
  `poisson`, `twoReg`, `contrastReg`)  
- `cate.model`: a formula describing the outcome model (e.g., outcome ~
  age + gender + previous_treatment)  
- `ps.model`: a formula describing the propensity score model to adjust
  for confounding (e.g., treatment ~ age + previous_treatment)

## Overview of the vignettes

[1. General
introduction](https://smartdata-analysis-and-statistics.github.io/precmed/articles/precmed.md)  
[2. Examples for count
outcome](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Count-examples.md)  
[3. Examples for survival
outcome](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Survival-examples.md)  
[4. Additional
examples](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Additional-examples.md)  
[5. Theoretical
details](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Theoretical-details.md)

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

Yadlowsky, Steve, Fabio Pellegrini, Federica Lionetto, Stefan Braune,
and Lu Tian. 2020. “Estimation and Validation of Ratio-Based Conditional
Average Treatment Effects Using Observational Data.” *Journal of the
American Statistical Association*, 1–18.

Zhao, Lihui, Lu Tian, Tianxi Cai, Brian Claggett, and Lee-Jen Wei. 2013.
“Effectively Selecting a Target Population for a Future Comparative
Study.” *Journal of the American Statistical Association* 108 (502):
527–39.
