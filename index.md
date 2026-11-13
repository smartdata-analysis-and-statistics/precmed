# precmed

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
randomized data from a clinical trial (Zhao et al. 2013) or using
real-world data (RWD) from a non-randomized study (Yadlowsky et al.
2021).

## Learning resources

For a broader introduction to personalized medicine using real-world
data, see *Comparative Effectiveness and Personalized Medicine Research
Using Real-World Data* (Debray, Nguyen, and Platt 2026). The following
chapters are particularly relevant to **precmed**:

- [Chapter 17](https://doi.org/10.1201/9781003300809-17): principles of
  CATE prediction and model development in randomized and non-randomized
  data (Hoogland et al. 2026).
- [Chapter 19](https://doi.org/10.1201/9781003300809-19): validation
  strategies and assessment of discrimination, calibration, and decision
  accuracy (Thao et al. 2026).
- [Chapter 20](https://doi.org/10.1201/9781003300809-20): visualization
  and interpretation of individualized treatment rules, validation
  curves, and responder diagnostics, with examples using **precmed**
  (Jiang et al. 2026).

The original methodological references are Yadlowsky et al. (2021) for
ratio-based CATE estimation and validation, and Zhao et al. (2013) for
target-population selection.

Two free courses complement the book and introduce **precmed**, its
practical use, and the underlying statistical methods. [An Introduction
to Precision Medicine](https://sdas.thinkific.com/courses/precmedintro)
covers the foundations through recorded webinars, exercises, and
readings. The [free companion course to the
book](https://sdas.thinkific.com/courses/rwd-book) provides background
materials, example R code, and curated references to support practical
analyses.

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

| Function | Description |
|----|----|
| [catefit](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.html)() | Estimation of the conditional average treatment effect (CATE) |
| [atefit](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.html)() | Doubly robust estimator for the average treatment effect (ATE) |
| [catecv](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.html)() | Development and cross-validation of the CATE |
| [abc](https://smartdata-analysis-and-statistics.github.io/precmed/reference/abc.html)() | Compute the area between the average treatment difference curve of competing models for the CATE (Zhao et al. 2013) |
| [plot](https://smartdata-analysis-and-statistics.github.io/precmed/reference/plot.precmed.html)() | Two side-by-side line plots of validation curves from the `precmed` object |
| [boxplot](https://smartdata-analysis-and-statistics.github.io/precmed/reference/boxplot.precmed.html)() | Plot the proportion of subjects with an estimated treatment effect no less than $`c`$ over a range of values for $`c`$ (Zhao et al. 2013). |

For more info:
<https://smartdata-analysis-and-statistics.github.io/precmed/>

## Recommended workflow

We recommend the following workflow to develop a model for estimating
the CATE in order to identify treatment effect heterogeneity:

1.  Compare up to five modelling approaches (e.g., Poisson regression,
    boosting) for estimating the CATE using cross-validation through
    [catecv](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.html).
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
      [abc](https://smartdata-analysis-and-statistics.github.io/precmed/reference/abc.html)
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

## User input

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

Debray, Thomas P. A., Tri-Long Nguyen, and Robert W. Platt, eds. 2026.
*Comparative Effectiveness and Personalized Medicine Research Using
Real-World Data*. 1st ed. Chapman & Hall/CRC.
<https://doi.org/10.1201/9781003300809>.

Hoogland, Jeroen, Thomas P. A. Debray, Guillaume Desachy, and Tri-Long
Nguyen. 2026. “Introduction to Methods for Personalizing Medicine.” In
*Comparative Effectiveness and Personalized Medicine Research Using
Real-World Data*, edited by Thomas P. A. Debray, Tri-Long Nguyen, and
Robert W. Platt, 421–63. Chapman & Hall/CRC.
<https://doi.org/10.1201/9781003300809-17>.

Jiang, Xiaotong, Gabrielle Simoneau, Fabio Pellegrini, Changyu Shen,
Carl de Moor, and Lu Tian. 2026. “Visualization and Interpretation of
Individualized Treatment Rule Results.” In *Comparative Effectiveness
and Personalized Medicine Research Using Real-World Data*, edited by
Thomas P. A. Debray, Tri-Long Nguyen, and Robert W. Platt, 513–36.
Chapman & Hall/CRC. <https://doi.org/10.1201/9781003300809-20>.

Thao, Le Thi Phuong, Orestis Efthimiou, Tibor V. Varga, Jeroen Hoogland,
and Thomas P. A. Debray. 2026. “Validation of Prediction Models for
Patient Outcomes and Individualized Treatment Effect.” In *Comparative
Effectiveness and Personalized Medicine Research Using Real-World Data*,
edited by Thomas P. A. Debray, Tri-Long Nguyen, and Robert W. Platt,
485–512. Chapman & Hall/CRC. <https://doi.org/10.1201/9781003300809-19>.

Yadlowsky, Steve, Fabio Pellegrini, Federica Lionetto, Stefan Braune,
and Lu Tian. 2021. “Estimation and Validation of Ratio-Based Conditional
Average Treatment Effects Using Observational Data.” *Journal of the
American Statistical Association* 116 (533): 335–52.
<https://doi.org/10.1080/01621459.2020.1772080>.

Zhao, Lihui, Lu Tian, Tianxi Cai, Brian Claggett, and Lee-Jen Wei. 2013.
“Effectively Selecting a Target Population for a Future Comparative
Study.” *Journal of the American Statistical Association* 108 (502):
527–39. <https://doi.org/10.1080/01621459.2013.770705>.
