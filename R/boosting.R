# Translate the existing GBM arguments to GBM3's current API. Keep the
# historical learning rate and use all observations/features for training.
fit_boosting <- function(formula, data, distribution, weights = NULL,
                         n.trees = 100, interaction.depth = 1,
                         n.minobsinnode = 10, shrinkage = 0.1,
                         bag.fraction = 0.5, train.fraction = 1,
                         cv.folds = 0, keep.data = TRUE, verbose = FALSE,
                         var.monotone = NULL, class.stratify.cv = FALSE,
                         n.cores = NULL, par.details = NULL,
                         fold.id = NULL, mFeatures = NULL) {
  if (!is.null(n.cores) && !is.null(par.details)) {
    stop("Specify only one of n.cores and par.details.")
  }
  if (is.null(par.details)) {
    par.details <- if (is.null(n.cores)) {
      getOption("gbm.parallel", gbm3::gbmParallel(num_threads = 1))
    } else {
      gbm3::gbmParallel(num_threads = n.cores)
    }
  }
  labels <- attr(stats::terms(formula, data = data), "term.labels")
  predictors <- model.frame(stats::reformulate(labels), data, na.action = stats::na.pass)
  if (is.null(mFeatures)) mFeatures <- ncol(predictors)
  if (is.null(weights)) weights <- rep(1, nrow(data))
  loss <- switch(distribution,
                 gaussian = "Gaussian", poisson = "Poisson",
                 stop("Unsupported boosting distribution: ", distribution))
  params <- gbm3::training_params(
    num_trees = n.trees, interaction_depth = interaction.depth,
    min_num_obs_in_node = n.minobsinnode, shrinkage = shrinkage,
    bag_fraction = bag.fraction, num_train = floor(train.fraction * nrow(data)),
    id = seq_len(nrow(data)), num_features = mFeatures)
  # Supply evaluated weights: model.frame otherwise resolves the symbol in
  # the caller's formula environment, where it can denote stats::weights.
  do.call(gbm3::gbmt, list(
    formula = formula, data = data, distribution = gbm3::gbm_dist(loss),
    weights = weights, train_params = params,
    cv_folds = if (cv.folds == 0) 1 else cv.folds,
    keep_gbm_data = keep.data, is_verbose = verbose,
    var_monotone = var.monotone, cv_class_stratify = class.stratify.cv,
    fold_id = fold.id, par_details = par.details))
}

best_boosting_iteration <- function(object, method, plot.it = TRUE) {
  performance <- gbm3::gbmt_performance(object, method = method)
  if (plot.it) plot(performance)
  as.numeric(performance)
}

#' Boosting controls and GBM3 backend
#'
#' Boosting uses \code{gbm3::gbmt()} for count, continuous, and survival
#' outcomes. Existing public argument names are translated internally.
#'
#' @details Additional arguments passed through \code{...} can include
#' \code{shrinkage} (default \code{0.1}), \code{n.minobsinnode} (\code{10}),
#' \code{bag.fraction} (\code{0.5}), \code{train.fraction} (\code{1}),
#' \code{var.monotone}, \code{mFeatures}, \code{keep.data},
#' \code{class.stratify.cv}, and \code{fold.id}.
#' Set \code{n.cores} to select the number of GBM3 threads, or supply
#' \code{par.details = gbm3::gbmParallel(...)}; do not supply both.
#' With neither argument, the \code{gbm.parallel} option is respected and
#' the fallback is one thread.
#'
#' GBM3's learning-rate and training-sample defaults are overridden to retain
#' precmed's previous GBM settings. Fits use the \code{GBMFit} class rather
#' than \code{gbm}; use the generic \code{predict()} method for prediction.
#' Predictions exclude exposure offsets, preserving the previous rate-scale
#' behavior. Results may differ from GBM even with the same settings and seed.
#' @name boosting-control
NULL
