
# precmed <img src="man/figures/logo.png" align="right" height="120"/>

<!-- badges: start -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/precmed)](https://cran.r-project.org/package=precmed)
[![metacran
downloads](https://cranlogs.r-pkg.org/badges/last-month/precmed)](https://cran.r-project.org/package=precmed)
[![R-CMD-check](https://github.com/smartdata-analysis-and-statistics/precmed/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/smartdata-analysis-and-statistics/precmed/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

**precmed** is an R package for estimating average treatment effects
(ATEs) and estimating and validating conditional average treatment
effects (CATEs) in randomized trials and observational studies, with
methods for binary, count and survival outcomes. The ATE describes an
effect across a population, while the CATE describes an average effect
among patients with particular observed characteristics. Together, they
support the study of treatment effect heterogeneity and precision
medicine.

Use **precmed** when you want to:

- Estimate the ATE comparing two treatments using doubly robust methods.
- Explore treatment effect heterogeneity by developing and assessing
  CATE models for binary, count and survival outcomes.
- Visualize and compare estimated effects across patient subgroups.

The methods support causal inference under appropriate assumptions; for
observational data, these include adequate control of confounding and
sufficient treatment overlap. The package combines outcome-regression
approaches, such as Poisson models fitted separately in each treatment
group, with doubly robust methods for estimating and validating
treatment effects. Its methods build on work by Zhao et al. (2013) and
Yadlowsky et al. (2021); Jiang et al. (2026) illustrate their
application and interpretation.

See [Is precmed suitable for
me?](https://smartdata-analysis-and-statistics.github.io/precmed/articles/quick-start.html#is-precmed-suitable-for-me)
for a quick overview of supported data, outcomes and analyses.

## Learning resources

For a broader introduction to personalized medicine using real-world
data, see *Comparative Effectiveness and Personalized Medicine Research
Using Real-World Data* (Debray, Nguyen, and Platt 2026). The following
chapters are particularly relevant to **precmed**:

- [Chapter 5](https://doi.org/10.1201/9781003300809-5): defining causal
  effects and understanding marginal and conditional treatment effects
  (Nguyen et al. 2026).
- [Chapter 17](https://doi.org/10.1201/9781003300809-17): principles of
  CATE prediction and model development in randomized and non-randomized
  data (Hoogland et al. 2026).
- [Chapter 19](https://doi.org/10.1201/9781003300809-19): validation
  strategies and assessment of discrimination, calibration, and decision
  accuracy (Thao et al. 2026).
- [Chapter 20](https://doi.org/10.1201/9781003300809-20): visualization
  and interpretation of individualized treatment rules, validation
  curves, and responder diagnostics, with examples using **precmed**
  (Jiang et al. 2026).

The original methodological references are Yadlowsky et al. (2021) for
ratio-based CATE estimation and validation, and Zhao et al. (2013) for
target-population selection.

Two free courses complement the book and introduce **precmed**, its
practical use, and the underlying statistical methods. [An Introduction
to Precision Medicine](https://sdas.thinkific.com/courses/precmedintro)
covers the foundations through recorded webinars, exercises, and
readings. The [free companion course to the
book](https://sdas.thinkific.com/courses/rwd-book) provides background
materials, example R code, and curated references to support practical
analyses.

<div class="learning-course-widgets">

<div class="thinkific-product-card"
data-btn-txt="Start learning for free" data-btn-txt-color="#ffffff"
data-btn-bg-color="#005a80" data-card-type="card"
data-link-type="landing_page" data-product="3112887"
data-embed-version="0.0.2" data-card-txt-color="#243746"
data-card-bg-color="#ffffff"
data-store-url="https://sdas.thinkific.com/embeds/products/show">

<div class="iframe-container">

</div>

</div>

<div class="thinkific-product-card"
data-btn-txt="Start learning for free" data-btn-txt-color="#ffffff"
data-btn-bg-color="#005a80" data-card-type="card"
data-link-type="landing_page" data-product="3591896"
data-embed-version="0.0.2" data-card-txt-color="#243746"
data-card-bg-color="#ffffff"
data-store-url="https://sdas.thinkific.com/embeds/products/show">

<div class="iframe-container">

</div>

</div>

</div>

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

Try the [quick
start](https://smartdata-analysis-and-statistics.github.io/precmed/articles/quick-start.html)
for a complete count-outcome example and validation plot.

## Citation and support

To cite the installed version of **precmed**, run `citation("precmed")`
in R, or see the [authors and citation
page](https://smartdata-analysis-and-statistics.github.io/precmed/authors.html#citation).
Please also cite the methodological references relevant to your
analysis.

For questions, bug reports, or feature requests, use [GitHub
issues](https://github.com/smartdata-analysis-and-statistics/precmed/issues).
For a bug report, please include your package and R versions, the error
message, and a minimal reproducible example using simulated or
non-sensitive data.

## Package capabilities

The package compares two treatments using outcome-specific effect
measures:

| Outcome | Effect measure | ATE example | CATE example |
|----|----|----|----|
| Binary | Risk difference in outcome probabilities (treatment 1 minus treatment 0) | [ATE](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Binary-examples.html) | [CATE](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Binary-CATE.html) |
| Count | Ratio of event rates (treatment 1 over treatment 0) | [ATE](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Count-examples.html) | [CATE](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Count-CATE.html) |
| Survival | At a specified follow-up horizon: population-average risk difference (treatment 1 minus treatment 0), ratio of restricted mean time lost, and population-level cumulative hazard ratio (treatment 1 over treatment 0) | [ATE](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Survival-examples.html) | [CATE](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Survival-CATE.html) |

Count rate ratios and survival RMTL and cumulative hazard ratios from
`atefit()` are returned on the log scale; exponentiate them to obtain
ratios. Survival risk differences are returned on the difference scale.
Count and survival CATE scores are log rate ratios and log RMTL ratios,
respectively. The ATE summarizes the population, while CATE scores
describe effects conditional on observed characteristics. Survival
subgroup assessment also reports risk differences at the follow-up
horizon.

The main functions in the **precmed** package are:

| Function | Description |
|----|----|
| [catefit](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catefit.html)() | Estimation of the conditional average treatment effect (CATE) |
| [atefit](https://smartdata-analysis-and-statistics.github.io/precmed/reference/atefit.html)() | Doubly robust estimator for the average treatment effect (ATE) |
| [catecv](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.html)() | Development and cross-validation of CATE scores |
| [abc](https://smartdata-analysis-and-statistics.github.io/precmed/reference/abc.html)() | Compute the area between the average treatment difference curve of competing models for the CATE (Zhao et al. 2013) |
| [plot](https://smartdata-analysis-and-statistics.github.io/precmed/reference/plot.precmed.html)() | Two side-by-side line plots of validation curves from the `precmed` object |
| [boxplot](https://smartdata-analysis-and-statistics.github.io/precmed/reference/boxplot.precmed.html)() | Plot the proportion of subjects with an estimated treatment effect no less than $c$ over a range of values for $c$ (Zhao et al. 2013). |

For more info:
<https://smartdata-analysis-and-statistics.github.io/precmed/>

## Recommended workflow

For count and survival outcomes, the following workflow estimates the
ATE and develops and validates a CATE model:

1.  Use `atefit()` to estimate the overall ATE between two treatment
    groups and its uncertainty using bootstrap resampling.
2.  Compare up to five modelling approaches (e.g., Poisson regression,
    boosting) for estimating the CATE using cross-validation through
    [catecv](https://smartdata-analysis-and-statistics.github.io/precmed/reference/catecv.html).
3.  Select the best modelling approach using 3 metrics:
    - Compare the steepness of the validation curves in the validation
      samples across methods using `plot()`. Two side-by-side plots are
      generated, visualizing the estimated average treatment effects in
      a series of nested subgroups. On the left side the curve is shown
      for the training set, and on the right side the curve is shown for
      the validation set. Each line in the plots represents one scoring
      method (e.g., boosting, randomForest) specified under the argument
      `score.method`.
    - The area between curves (ABC) using
      [abc](https://smartdata-analysis-and-statistics.github.io/precmed/reference/abc.html)
      quantifies a model’s ability to capture treatment effect
      heterogeneity. Higher ABC values are preferable as they indicate
      that more treatment effect heterogeneity is captured by the
      scoring method.
    - Compare the distribution of the estimated ATE across different
      levels of the CATE score percentiles using `boxplot()`.
4.  Fit the chosen modelling approach using `catefit()`; assess
    performance in an independent dataset where possible.

For binary outcomes, the binary example illustrates `atefit()`,
`catefit()` and a separate held-out subgroup assessment.

In the vignettes, we will adopt a different workflow to gradually expose
the user from simple to more complex methods.

## Vignettes

- [Understanding doubly robust
  methods](https://smartdata-analysis-and-statistics.github.io/precmed/articles/doubly-robust.html)
- [Binary ATE
  example](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Binary-examples.html)
- [Binary CATE
  example](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Binary-CATE.html)
- [Count ATE
  example](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Count-examples.html)
- [Count CATE
  example](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Count-CATE.html)
- [Survival ATE
  example](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Survival-examples.html)
- [Survival CATE
  example](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Survival-CATE.html)
- [Customizing CATE models and
  validation](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Additional-examples.html)
- [ATE
  methods](https://smartdata-analysis-and-statistics.github.io/precmed/articles/ATE-theory.html)
- [CATE development and validation
  methods](https://smartdata-analysis-and-statistics.github.io/precmed/articles/Theoretical-details.html)

# References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-debray2026comparative" class="csl-entry">

Debray, Thomas P. A., Tri-Long Nguyen, and Robert W. Platt, eds. 2026.
*Comparative Effectiveness and Personalized Medicine Research Using
Real-World Data*. 1st ed. Chapman & Hall/CRC.
<https://doi.org/10.1201/9781003300809>.

</div>

<div id="ref-hoogland2026introduction" class="csl-entry">

Hoogland, Jeroen, Thomas P. A. Debray, Guillaume Desachy, and Tri-Long
Nguyen. 2026. “Introduction to Methods for Personalizing Medicine.” In
*Comparative Effectiveness and Personalized Medicine Research Using
Real-World Data*, edited by Thomas P. A. Debray, Tri-Long Nguyen, and
Robert W. Platt, 421–63. Chapman & Hall/CRC.
<https://doi.org/10.1201/9781003300809-17>.

</div>

<div id="ref-jiang2026visualization" class="csl-entry">

Jiang, Xiaotong, Gabrielle Simoneau, Fabio Pellegrini, Changyu Shen,
Carl de Moor, and Lu Tian. 2026. “Visualization and Interpretation of
Individualized Treatment Rule Results.” In *Comparative Effectiveness
and Personalized Medicine Research Using Real-World Data*, edited by
Thomas P. A. Debray, Tri-Long Nguyen, and Robert W. Platt, 513–36.
Chapman & Hall/CRC. <https://doi.org/10.1201/9781003300809-20>.

</div>

<div id="ref-nguyen2026causaleffects" class="csl-entry">

Nguyen, Tri-Long, Jessica L. Rohmann, Tat Thang Vo, Robert W. Platt, and
Thomas P. A. Debray. 2026. “Understanding and Defining Causal Effects.”
In *Comparative Effectiveness and Personalized Medicine Research Using
Real-World Data*, edited by Thomas P. A. Debray, Tri-Long Nguyen, and
Robert W. Platt, 65–80. Chapman & Hall/CRC.
<https://doi.org/10.1201/9781003300809-5>.

</div>

<div id="ref-thao2026validation" class="csl-entry">

Thao, Le Thi Phuong, Orestis Efthimiou, Tibor V. Varga, Jeroen Hoogland,
and Thomas P. A. Debray. 2026. “Validation of Prediction Models for
Patient Outcomes and Individualized Treatment Effect.” In *Comparative
Effectiveness and Personalized Medicine Research Using Real-World Data*,
edited by Thomas P. A. Debray, Tri-Long Nguyen, and Robert W. Platt,
485–512. Chapman & Hall/CRC. <https://doi.org/10.1201/9781003300809-19>.

</div>

<div id="ref-yadlowsky2020estimation" class="csl-entry">

Yadlowsky, Steve, Fabio Pellegrini, Federica Lionetto, Stefan Braune,
and Lu Tian. 2021. “Estimation and Validation of Ratio-Based Conditional
Average Treatment Effects Using Observational Data.” *Journal of the
American Statistical Association* 116 (533): 335–52.
<https://doi.org/10.1080/01621459.2020.1772080>.

</div>

<div id="ref-zhao2013effectively" class="csl-entry">

Zhao, Lihui, Lu Tian, Tianxi Cai, Brian Claggett, and Lee-Jen Wei. 2013.
“Effectively Selecting a Target Population for a Future Comparative
Study.” *Journal of the American Statistical Association* 108 (502):
527–39. <https://doi.org/10.1080/01621459.2013.770705>.

</div>

</div>
