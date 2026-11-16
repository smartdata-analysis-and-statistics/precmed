# Fit Bernoulli boosting separately in each treatment arm on a common design.
# Author: Thomas Debray <tdebray@fromdatatowisdom.com>
binary_boosting <- function(d, tree.depth, n.trees.boosting, plot.gbmperf, ...) {
  if (ncol(d$x) < 2L) stop("Binary boosting requires at least one outcome-model covariate.")
  if (!is.numeric(tree.depth) || length(tree.depth) != 1L || !is.finite(tree.depth) ||
      tree.depth < 1 || tree.depth != floor(tree.depth)) stop("tree.depth must be a positive integer.")
  if (!is.numeric(n.trees.boosting) || length(n.trees.boosting) != 1L || !is.finite(n.trees.boosting) ||
      n.trees.boosting < 1 || n.trees.boosting != floor(n.trees.boosting)) stop("n.trees.boosting must be a positive integer.")
  if (!is.logical(plot.gbmperf) || length(plot.gbmperf) != 1L || is.na(plot.gbmperf)) stop("plot.gbmperf must be TRUE or FALSE.")
  predictors <- as.data.frame(d$x[, -1, drop = FALSE])
  names(predictors) <- paste0("x", seq_len(ncol(predictors)))
  train <- data.frame(y = d$y, predictors, check.names = FALSE)
  fits <- vector("list", 2L); best <- integer(2L); predictions <- vector("list", 2L)
  for (arm in 0:1) {
    selected <- d$trt == arm
    if (sum(selected) < 5L || length(unique(d$y[selected])) != 2L)
      stop("Binary boosting requires at least five patients and both outcome values in each treatment group.")
    fits[[arm + 1L]] <- tryCatch(
      fit_boosting(y ~ ., data = train[selected, , drop = FALSE], distribution = "bernoulli",
                   interaction.depth = tree.depth, n.trees = n.trees.boosting, cv.folds = 5, ...),
      error = function(e) stop("Binary boosting failed for treatment ", arm, ": ", conditionMessage(e), call. = FALSE))
    best[[arm + 1L]] <- min(n.trees.boosting, max(10, best_boosting_iteration(fits[[arm + 1L]], "cv", plot.it = plot.gbmperf)))
    predictions[[arm + 1L]] <- as.numeric(predict(fits[[arm + 1L]], newdata = predictors,
                                                  n.trees = best[[arm + 1L]], type = "response"))
  }
  if (any(!is.finite(unlist(predictions))) || any(unlist(predictions) < 0 | unlist(predictions) > 1))
    stop("Binary boosting produced invalid predicted probabilities.")
  list(mu0 = predictions[[1]], mu1 = predictions[[2]],
       fit = list(fit0.boosting = fits[[1]], fit1.boosting = fits[[2]],
                  best0.iter = best[[1]], best1.iter = best[[2]],
                  predictor.names = colnames(d$x)[-1]))
}

# Author: Thomas Debray <tdebray@fromdatatowisdom.com>
catefitbinary <- function(data, score.method, cate.model, ps.model, ps.method = "glm", minPS = 0.01, maxPS = 0.99, higher.y = TRUE, prop.cutoff = seq(0.5, 1, length = 6), seed = NULL, ps.control = list(),
                          tree.depth = 2, n.trees.boosting = 200, plot.gbmperf = TRUE, ...) {
  if (!ps.method %in% c("glm", "lasso") || minPS < 0 || maxPS > 1 || minPS >= maxPS) stop("Invalid propensity score settings.")
  resolve_ps_control(ps.control)
  if (!is.logical(higher.y) || length(higher.y) != 1L || is.na(higher.y)) stop("higher.y must be TRUE or FALSE.")
  if (!is.numeric(prop.cutoff) || !length(prop.cutoff) || anyNA(prop.cutoff) || any(prop.cutoff <= 0 | prop.cutoff > 1)) stop("prop.cutoff must contain values in (0, 1].")
  d <- binary_data(data, cate.model, ps.model)
  if (!is.null(seed)) set.seed(seed)
  prop <- sort(prop.cutoff); out <- list(response = "binary", score.method = score.method, fit = list(), coefficients = NULL)
  for (method in score.method) {
    pred <- if (method == "boosting") binary_boosting(d, tree.depth, n.trees.boosting, plot.gbmperf, ...) else
      binary_predictions(d$y, d$trt, d$x, method)
    score <- pred$mu1 - pred$mu0
    out[[paste0("score.", method)]] <- score
    out$fit[[method]] <- pred$fit
    ranking <- if (higher.y) score else -score
    effects <- vapply(prop, function(p) {
      selected <- ranking >= unname(stats::quantile(ranking, 1 - p))
      tryCatch(unname(binary_effect(binary_subset(d, selected), ps.method, minPS, maxPS, ps.control)[3]), error = function(e) { warning("Subgroup effect unavailable: ", conditionMessage(e)); NA_real_ })
    }, numeric(1))
    names(effects) <- paste0("prop", round(prop, 2))
    out[[paste0("ate.", method)]] <- effects
  }
  structure(out, class = "catefit")
}
