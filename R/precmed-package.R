#' precmed: Estimation and validation of heterogeneous treatment effects
#'
#' Tools for estimating population-average and covariate-conditional treatment
#' effects in randomized and observational studies with binary, count, or
#' time-to-event outcomes. Causal interpretation requires appropriate
#' identification and modeling assumptions.
#'
#' @details
#' Use \code{\link{atefit}()} to estimate population-average treatment effects. Use
#' \code{\link{catefit}()} to develop conditional treatment-effect scores and \code{\link{catecv}()}
#' to evaluate scores using cross-validation. Estimating a population-average
#' effect and validating a conditional treatment-effect score are distinct tasks.
#' Use \code{\link{diagnostics}()} to inspect available diagnostics and
#' \code{\link{subgroup_estimates}()} to extract subgroup results.
#'
#' @section Outcomes and effect scales:
#' For binary outcomes, \code{atefit()} estimates the population-average risk
#' difference (treatment 1 minus treatment 0). For count outcomes it estimates
#' a population-average rate ratio. For survival outcomes at a chosen horizon,
#' it estimates a restricted mean time lost ratio, cumulative hazard ratio,
#' and event risk difference (treatment 1 minus treatment 0). The CATE scores
#' describe effects conditional on covariates; their scales depend on the
#' outcome. See \code{\link{atefit}} and \code{\link{catefit}} for details.
#'
#' @section The question and the models:
#' In observational studies, treatment groups may differ in baseline
#' characteristics that also affect outcomes. A treatment model estimates
#' propensity scores, the probabilities of treatment assignment given
#' covariates. An outcome model predicts the outcome given treatment and
#' covariates. These models can be used to adjust for measured covariates when
#' estimating treatment effects. For survival outcomes, a censoring model
#' supplies probabilities used for inverse probability of censoring weighting.
#'
#' @section What does doubly robust mean?:
#' \code{\link{atefit}()} uses doubly robust estimators of population-average
#' treatment effects for binary, count, and survival outcomes. These estimators
#' combine outcome and treatment models so that, under the required causal
#' identification and regularity assumptions, they are consistent if either
#' the outcome model or the treatment model is correctly specified. Both models
#' may be correct.
#' Double robustness does not protect against unmeasured confounding, poor
#' treatment overlap, or misspecification of both models. The property is
#' asymptotic and does not guarantee negligible bias in a finite sample.
#' For survival outcomes, it also does not protect against a misspecified
#' censoring model. If an out-of-range restricted mean survival estimate
#' triggers an outcome-model fallback, the substituted estimate does not retain
#' the doubly robust guarantee. See \code{\link{atefit}} for details.
#'
#' @section Online guides and vignettes:
#' The package website provides a
#' \href{https://smartdata-analysis-and-statistics.github.io/precmed/articles/quick-start.html}{quick start},
#' an introduction to
#' \href{https://smartdata-analysis-and-statistics.github.io/precmed/articles/doubly-robust.html}{doubly robust methods},
#' and \href{https://smartdata-analysis-and-statistics.github.io/precmed/articles/index.html}{all vignettes},
#' including examples and methodological details for different outcomes.
#' These resources are available online; the vignettes are not included in the
#' CRAN package build.
#'
#' @references
#' Debray, T. P. A., Nguyen, T.-L., and Platt, R. W. (eds.) (2026).
#' \emph{Comparative Effectiveness and Personalized Medicine Research Using
#' Real-World Data}. Chapman & Hall/CRC.
#' \doi{10.1201/9781003300809}.
#'
#' @examples
#' # A small illustrative example; use more bootstrap replicates in practice.
#' dat <- simulate_binary(n = 400, seed = 2026)
#' fit <- atefit(response = "binary", data = dat,
#'               cate.model = y ~ age + previous_treatment,
#'               ps.model = trt ~ age + previous_treatment,
#'               n.boot = 20, seed = 2026)
#' coef(fit)
#' diagnostics(fit)
#'
#' @keywords internal
"_PACKAGE"
