# Learning resources

These resources provide background on personalized medicine,
treatment-effect prediction and the methods used in **precmed**.

## Handbook

For a broader introduction to personalized medicine using real-world
data, see [*Comparative Effectiveness and Personalized Medicine Research
Using Real-World Data*](https://doi.org/10.1201/9781003300809) ([Debray
et al. 2026](#ref-debray2026comparative)). The following chapters are
particularly relevant to **precmed**:

- [Chapter 5](https://doi.org/10.1201/9781003300809-5): marginal and
  conditional treatment effects, including the ATE ([Nguyen et al.
  2026](#ref-nguyen2026causaleffects)).
- [Chapter 17](https://doi.org/10.1201/9781003300809-17): principles of
  CATE prediction and model development in randomized and non-randomized
  data ([Hoogland et al. 2026](#ref-hoogland2026introduction)).
- [Chapter 19](https://doi.org/10.1201/9781003300809-19): validation
  strategies and assessment of discrimination, calibration, and decision
  accuracy ([Thao et al. 2026](#ref-thao2026validation)).
- [Chapter 20](https://doi.org/10.1201/9781003300809-20): visualization
  and interpretation of individualized treatment rules, validation
  curves, and responder diagnostics, with examples using **precmed**
  ([Jiang et al. 2026](#ref-jiang2026visualization)).

## Free courses

These free courses can help you build the knowledge and practical skills
needed to study personalized treatment decisions. [An Introduction to
Precision Medicine](https://sdas.thinkific.com/courses/precmedintro)
introduces the foundations through recorded webinars, exercises, and
readings. The [companion course to the
book](https://sdas.thinkific.com/courses/rwd-book) provides background
materials, example R code, and curated references to support practical
analyses with **precmed** and related methods.

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

## Methodological references

The original methodological references are Yadlowsky et al.
([2021](#ref-yadlowsky2020estimation)) for ratio-based CATE estimation
and validation, and Zhao et al. ([2013](#ref-zhao2013effectively)) for
target-population selection.

Debray, Thomas P. A., Tri-Long Nguyen, and Robert W. Platt, eds. 2026.
*Comparative Effectiveness and Personalized Medicine Research Using
Real-World Data*. 1st ed. Chapman & Hall/CRC.
<https://doi.org/10.1201/9781003300809>.

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
