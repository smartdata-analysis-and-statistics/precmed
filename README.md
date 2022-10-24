
# precmed <img src="man/figures/logo.png" align="right" height="120"/>

<!-- badges: start -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/precmed)](https://cran.r-project.org/package=precmed)
[![metacran
downloads](https://cranlogs.r-pkg.org/badges/last-month/precmed)](https://cran.r-project.org/package=precmed)
<!-- badges: end -->

## Overview

**precmed** was developed to help researchers with precision medicine.
In precision medicine, we want to determine the optimal treatment for
each patient instead of applying one treatment to all patients because
not all patients respond to treatments in the same way. This is also
called *treatment effect heterogeneity* in which different subgroups of
patients respond differently to the same treatment. **precmed** can help
to estimate and internally validate these different treatment effects,
also called the conditional average treatment effects (CATEs) which
describe how different subgroups respond to the same treatment. Users of
the package can estimate these effects (CATEs) between two treatments
(e.g., treatment A or treatment B) from either randomized controlled
trials (RCTs) or real world data (RWD).

## Installation

``` r
# The `precmed` package can be installed from CRAN:
install.packages("precmed")

# Alternatively, the lastest version can be installed from github:
install.packages("devtools")
devtools::install_github(repo = "smartdata-analysis-and-statistics/precmed")
```

## Capabilities of `precmed` package

The `precmed` package contains functions to

-   Estimate average treatment effects for count, survival and contiuous
    data
-   Estimate conditional average treatment effect (CATE) score for
    count, survival and continuous data
-   Cross-validation of CATE score
-   Compute area between curves
-   Plot validation curves

## Main functions

The main functions in the `precmed` package are:

| Function name | Description                                                                                                          |
|---------------|----------------------------------------------------------------------------------------------------------------------|
| catecv()      | Cross-validation of the conditional average treatment effect (CATE) score for count, survival or continuous outcomes |
| catefit()     | Estimation of the conditional average treatment effect (CATE) score for count, survival and continuous data          |
| atefit()      | Doubly robust estimator of and inference for the average treatment effect for count, survival and continuous data    |
| abc()         | Compute the area between curves from the `precmed` object                                                            |
| plot()        | Two side-by-side line plots of validation curves from the `precmed` object                                           |
| boxplot()     | A set of box plots of estimated ATEs from the `precmed` object                                                       |

## Background materials

1.  [Zhao et al. (2013), Effectively selecting a target population for a
    future comparative
    study](https://doi.org/10.1080/01621459.2013.770705)

2.  [Yadlowsky et al. (2020), Estimation and validation of ratio-based
    conditional average treatment effects using observational
    data](https://doi.org/10.1080/01621459.2020.1772080)

For more info:
<https://smartdata-analysis-and-statistics.github.io/precmed/>

## Vignettes

1.  [General introduction to
    precmed](https://smartdata-analysis-and-statistics.github.io/precmed/docs/articles/general-introduction.html)  
2.  [Examples with count outcome of the entire
    workflow](https://smartdata-analysis-and-statistics.github.io/precmed/docs/articles/count-examples.html)  
3.  [Examples with survival outcome of the entire
    workflow](https://smartdata-analysis-and-statistics.github.io/precmed/docs/articles/survival-examples.html)  
4.  [Additional examples for the `precmed`
    package](https://smartdata-analysis-and-statistics.github.io/precmed/docs/articles/additional-examples.html)  
5.  [Theoretical
    details](https://smartdata-analysis-and-statistics.github.io/precmed/docs/articles/theoretical-details.html)
