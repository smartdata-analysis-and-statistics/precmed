
# precmed <img src="man/figures/logo.png" align="right" height="120"/>

<!-- badges: start -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/precmed)](https://cran.r-project.org/package=precmed)
[![metacran
downloads](https://cranlogs.r-pkg.org/badges/last-month/precmed)](https://cran.r-project.org/package=precmed)
[![DOI](https://zenodo.org/badge/523236629.svg)](https://zenodo.org/badge/latestdoi/523236629)
[![R-CMD-check](https://github.com/amices/ggmice/workflows/R-CMD-check/badge.svg)](https://github.com/smartdata-analysis-and-statistics/precmed/actions)
<!-- badges: end -->

## Overview

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
randomized data from a clinical trial (Zhao et al. 2013) or using
real-world data (RWD) from a non-randomized study (Yadlowsky et al.
2020).

## Installation

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

The main functions in the **precmed** package are:

| Function                                                                                                | Description                                                                                                                            |
|---------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------|
| [catefit](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.html)()         | Estimation of the conditional average treatment effect (CATE)                                                                          |
| [atefit](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.html)()           | Doubly robust estimator for the average treatment effect (ATE)                                                                         |
| [catecv](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.html)()           | Development and cross-validation of the CATE                                                                                           |
| [abc](https://smartdata-analysis-and-statistics.github.io/precmed/reference/abc.html)()                 | Compute the area between the average treatment difference curve of competing models for the CATE (Zhao et al. 2013)                    |
| [plot](https://smartdata-analysis-and-statistics.github.io/precmed/reference/plot.precmed.html)()       | Two side-by-side line plots of validation curves from the `precmed` object                                                             |
| [boxplot](https://smartdata-analysis-and-statistics.github.io/precmed/reference/boxplot.precmed.html)() | Plot the proportion of subjects with an estimated treatment effect no less than $c$ over a range of values for $c$ (Zhao et al. 2013). |

For more info:
<https://smartdata-analysis-and-statistics.github.io/precmed/>

## Recommended workflow

We recommend the following workflow to develop a model for estimating
the CATE in order to identify treatment effect heterogeneity:

1.  Compare up to five modelling approaches (e.g., Poisson regression,
    boosting) for estimating the CATE using cross-validation through
    [catecv](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.html).
2.  Select the best modelling approach using 3 metrics:
    -   Compare the steepness of the validation curves in the validation
        samples across methods using `plot()`. Two side-by-side plots
        are generated, visualizing the estimated average treatment
        effects in a series of nested subgroups. On the left side the
        curve is shown for the training set, and on the right side the
        curve is shown for the validation set. Each line in the plots
        represents one scoring method (e.g., boosting, randomForest)
        specified under the argument `score.method`.
    -   The area between curves (ABC) using
        [abc](https://smartdata-analysis-and-statistics.github.io/precmed/reference/abc.html)
        quantifies a model’s ability to capture treatment effect
        heterogeneity. Higher ABC values are preferable as they indicate
        that more treatment effect heterogeneity is captured by the
        scoring method.
    -   Compare the distribution of the estimated ATE across different
        levels of the CATE score percentiles using `boxplot()`.
3.  Apply the best modelling approach in the original data or in a new
    external dataset using `catefit()`.
4.  Optional. Use `atefit()` to estimate ATE between 2 treatment groups
    with a doubly robust estimator and estimate the variability of the
    ATE with a bootstrap approach.

In the vignettes, we will adopt a different workflow to gradually expose
the user from simple to more complex methods.

## User input

When applying `catefit()` or `catecv()`, the user has to (at least)
input:

-   `response`: type of outcome/response (either `count` or
    `survival`)  
-   `data`: a data frame with individual patient data  
-   `score.method`: methods to estimate the CATE (e.g., `boosting`,
    `poisson`, `twoReg`, `contrastReg`)  
-   `cate.model`: a formula describing the outcome model (e.g., outcome
    \~ age + gender + previous_treatment)  
-   `ps.model`: a formula describing the propensity score model to
    adjust for confounding (e.g., treatment \~ age + previous_treatment)

## Vignettes

1.  [Examples with count outcome of the entire
    workflow](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Count-examples.html)  
2.  [Examples with survival outcome of the entire
    workflow](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Survival-examples.html)  
3.  [Additional examples for the `precmed`
    package](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Additional-examples.html)  
4.  [Theoretical
    details](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Theoretical-details.html)

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-yadlowsky2020estimation" class="csl-entry">

Yadlowsky, Steve, Fabio Pellegrini, Federica Lionetto, Stefan Braune,
and Lu Tian. 2020. “Estimation and Validation of Ratio-Based Conditional
Average Treatment Effects Using Observational Data.” *Journal of the
American Statistical Association*, 1–18.

</div>

<div id="ref-zhao2013effectively" class="csl-entry">

Zhao, Lihui, Lu Tian, Tianxi Cai, Brian Claggett, and Lee-Jen Wei. 2013.
“Effectively Selecting a Target Population for a Future Comparative
Study.” *Journal of the American Statistical Association* 108 (502):
527–39. <https://doi.org/10.1080/01621459.2013.770705>.

</div>

</div>
