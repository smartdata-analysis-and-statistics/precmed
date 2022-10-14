<!-- README.md is generated from README.Rmd using knitr. Please edit that file -->

# precmed: Precision Medicine in R <img src='https://fromdatatowisdom.com/images/precmed_sticker.jpg' align="right" height="139" />

<!-- badges: start -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/precmed)](https://cran.r-project.org/package=precmed)
<!-- badges: end -->

A doubly robust precision medicine approach to estimate and validate
conditional average treatment effects

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
## Capabilities of `precmed` package

The `precmed` package contains functions to

- Estimate average treatment effects for count, survival and contiuous data
- Estimate conditional average treatment effect (CATE) score for count, survival and continuous data
- Cross-validation of CATE score
- Compute area between curves
- Plot validation curves


## Main functions

The main functions in the `precmed` package are:

Function name      | Description
-------------------|---------------------------------
atefit()           | Doubly robust estimator of and inference for the average treatment effect for count, survival and continuous data
catefit()          | Estimation of the conditional average treatment effect (CATE) score for count, survival and continuous data
catecv()           | Cross-validation of the conditional average treatment effect (CATE) score for count, survival or continuous outcomes
abc()              | Compute the area between curves from the `precmed` object
plot()             | Two side-by-side line plots of validation curves from the `precmed` object
boxplot()          | A set of box plots of estimated ATEs from the `precmed` object

## Background materials

1. [Zhao et al. (2013), Effectively selecting a target population for a future comparative study](https://doi.org/10.1080/01621459.2013.770705)

2. [Yadlowsky et al. (2020), Estimation and validation of ratio-based conditional average treatment effects using observational data](https://doi.org/10.1080/01621459.2020.1772080)


## Vignettes

Under construction

1. General introduction to precmed
2. Example for count outcome
3. Example for survival outcome
4. Detailed examples
5. Theoretical details
