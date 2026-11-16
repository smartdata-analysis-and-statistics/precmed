# Predict risk-difference scores using models fitted in the training split.
# Author: Thomas Debray <tdebray@fromdatatowisdom.com>
binary_score_new <- function(fit, x, method) {
  if (method == "twoReg") {
    beta <- fit$coefficients
    return(drop(x %*% (beta[, 2] - beta[, 1])))
  }
  if (method == "logistic") {
    base <- cbind(x, trt = 0, if (ncol(x) > 1L) x[, -1, drop = FALSE] * 0)
    treated <- cbind(x, trt = 1, if (ncol(x) > 1L) x[, -1, drop = FALSE])
    beta <- fit$control$coefficients
    return(stats::plogis(drop(treated %*% beta)) - stats::plogis(drop(base %*% beta)))
  }
  predictors <- as.data.frame(x[, -1, drop = FALSE])
  names(predictors) <- paste0("x", seq_len(ncol(predictors)))
  p0 <- predict(fit$fit0.boosting, newdata = predictors, n.trees = fit$best0.iter, type = "response")
  p1 <- predict(fit$fit1.boosting, newdata = predictors, n.trees = fit$best1.iter, type = "response")
  as.numeric(p1 - p0)
}

# Build both design matrices using terms and transformation parameters learned in training.
# Author: Thomas Debray <tdebray@fromdatatowisdom.com>
binary_cv_design <- function(data, train, formula) {
  frame <- stats::model.frame(formula, data = data[train, , drop = FALSE], na.action = stats::na.pass)
  terms <- attr(frame, "terms")
  xtrain <- stats::model.matrix(terms, frame)
  xlevels <- stats::.getXlevels(terms, frame)
  valid <- stats::model.frame(stats::delete.response(terms), data = data[!train, , drop = FALSE],
                              na.action = stats::na.pass, xlev = xlevels)
  xvalid <- stats::model.matrix(stats::delete.response(terms), valid, contrasts.arg = attr(xtrain, "contrasts"))
  if (!identical(colnames(xtrain), colnames(xvalid)) || any(!is.finite(xtrain)) || any(!is.finite(xvalid)))
    stop("Training and validation model matrices are incompatible or non-finite.")
  list(train = xtrain, valid = xvalid)
}

# Repeated treatment-stratified validation for binary risk-difference scores.
# Author: Thomas Debray <tdebray@fromdatatowisdom.com>
catecvbinary <- function(data, score.method, cate.model, ps.model, ps.method = "glm",
                         minPS = 0.01, maxPS = 0.99, higher.y = TRUE,
                         prop.cutoff = seq(0.5, 1, length = 6),
                         prop.multi = c(0, 1/3, 2/3, 1), abc = TRUE,
                         train.prop = 3/4, cv.n = 10, tree.depth = 2,
                         n.trees.boosting = 200, B = 3, Kfold = 5,
                         seed = NULL, plot.gbmperf = TRUE, verbose = 0, ps.control = list()) {
  check_score_method("binary", score.method)
  if (anyDuplicated(score.method)) stop("score.method must not contain duplicates.")
  if (length(ps.method) != 1L || is.na(ps.method) || !ps.method %in% c("glm", "lasso") ||
      !is.numeric(minPS) || !is.numeric(maxPS) || length(minPS) != 1L || length(maxPS) != 1L || !is.finite(minPS) || !is.finite(maxPS) ||
      minPS < 0 || maxPS > 1 || minPS >= maxPS) stop("Invalid propensity score settings.")
  resolve_ps_control(ps.control)
  if (!is.logical(higher.y) || length(higher.y) != 1L || is.na(higher.y) ||
      !is.logical(abc) || length(abc) != 1L || is.na(abc)) stop("higher.y and abc must be TRUE or FALSE.")
  if (!is.numeric(prop.cutoff) || !length(prop.cutoff) || any(!is.finite(prop.cutoff)) ||
      any(prop.cutoff <= 0 | prop.cutoff > 1) || anyDuplicated(prop.cutoff))
    stop("prop.cutoff must contain distinct values in (0, 1].")
  if (!is.numeric(prop.multi) || length(prop.multi) < 3L || any(!is.finite(prop.multi)) ||
      prop.multi[1] != 0 || prop.multi[length(prop.multi)] != 1 || any(diff(prop.multi) <= 0))
    stop("prop.multi must increase strictly from 0 to 1 and contain at least three values.")
  if (!is.numeric(train.prop) || length(train.prop) != 1L || !is.finite(train.prop) || train.prop <= 0 || train.prop >= 1 ||
      !is.numeric(cv.n) || length(cv.n) != 1L || !is.finite(cv.n) || cv.n < 1 || cv.n != floor(cv.n))
    stop("train.prop must be in (0, 1) and cv.n must be a positive integer.")
  if (!is.numeric(verbose) || length(verbose) != 1L || is.na(verbose) || !verbose %in% 0:2) stop("verbose must be 0, 1 or 2.")
  prop <- sort(prop.cutoff); bi <- prop[prop < 1]
  if (abc && length(prop) < 2L) stop("ABC requires at least two distinct prop.cutoff values.")
  d <- binary_data(data, cate.model, ps.model)
  narm <- tabulate(d$trt + 1L, 2L)
  ntrain <- floor(narm * train.prop)
  if (any(ntrain < 1 | ntrain >= narm)) stop("Each treatment group must have observations in both training and validation sets.")
  if (!is.null(seed)) set.seed(seed)
  result <- stats::setNames(vector("list", length(score.method)), paste0("ate.", score.method))
  # Author: Thomas Debray <tdebray@fromdatatowisdom.com>
  matrix_cv <- function(rows, labels) {
    matrix(NA_real_, rows, cv.n, dimnames = list(labels, paste0("cv", seq_len(cv.n))))
  }
  for (method in score.method) {
    entry <- list()
    for (split in c("train", "valid")) {
      entry[[paste0("ate.est.", split, ".high.cv")]] <- matrix_cv(length(prop), paste0("prop", round(prop, 2)))
      entry[[paste0("ate.est.", split, ".low.cv")]] <- matrix_cv(length(bi), paste0("prop", round(1 - bi, 2)))
      entry[[paste0("ate.est.", split, ".group.cv")]] <- matrix_cv(length(prop.multi) - 1L, paste0("prop", round(prop.multi[-1], 2)))
    }
    if (abc) entry$abc.valid <- rep(NA_real_, cv.n)
    result[[paste0("ate.", method)]] <- entry
  }
  messages <- data.frame(method = character(), split = character(), iteration = integer(), group = character(), subgroup = character(), type = character(), measure = character(), message = character(), stringsAsFactors = FALSE)
  support <- data.frame(method = character(), split = character(), iteration = integer(), group = character(), subgroup = character(), n0 = integer(), n1 = integer(), events0 = numeric(), events1 = numeric(), stringsAsFactors = FALSE)
  # Author: Thomas Debray <tdebray@fromdatatowisdom.com>
  record <- function(method, split, iteration, group, subgroup, type, message) {
    split <- if (identical(split, "train")) "training" else if (identical(split, "valid")) "validation" else split
    messages[nrow(messages) + 1L, ] <<- list(method, split, iteration, group, subgroup, type, "risk difference", message)
  }
  # Author: Thomas Debray <tdebray@fromdatatowisdom.com>
  effect <- function(part, selected, method, split, iteration, group, label, report = TRUE) {
    counts <- c(sum(selected & part$trt == 0), sum(selected & part$trt == 1))
    if (report) support[nrow(support) + 1L, ] <<- list(method, if (split == "train") "training" else "validation", iteration, group, label, counts[1], counts[2], sum(part$y[selected & part$trt == 0]), sum(part$y[selected & part$trt == 1]))
    tryCatch({
      if (any(counts == 0)) stop("subgroup lacks one or both treatment groups")
      value <- unname(binary_effect(binary_subset(part, selected), ps.method, minPS, maxPS, ps.control)[3])
      if (!is.finite(value)) stop("subgroup estimate is not finite")
      value
    }, error = function(e) {
      if (report) {
        msg <- paste0(method, ", ", split, " iteration ", iteration, ", ", group, " ", label, ": ", conditionMessage(e), "; returning NA.")
        record(method, split, iteration, group, label, "error", msg)
        warning(msg, call. = FALSE)
      }
      NA_real_
    })
  }
  result$overall.ate.train <- result$overall.ate.valid <- rep(NA_real_, cv.n)
  result$splits <- vector("list", cv.n)
  for (i in seq_len(cv.n)) {
    if (verbose >= 1) message("Binary validation iteration ", i, "/", cv.n)
    train <- rep(FALSE, length(d$y))
    for (a in 0:1) {
      ids <- which(d$trt == a)
      train[ids[sample.int(length(ids), ntrain[a + 1L])]] <- TRUE
    }
    result$splits[[i]] <- list(train = which(train), valid = which(!train))
    parts <- list(train = binary_subset(d, train), valid = binary_subset(d, !train))
    designs <- tryCatch(list(x = binary_cv_design(data, train, cate.model), z = binary_cv_design(data, train, ps.model)), error = function(e) {
      msg <- paste0("Binary validation iteration ", i, ": design construction failed: ", conditionMessage(e))
      record(NA_character_, "train", i, "design", NA_character_, "error", msg)
      warning(msg, call. = FALSE)
      NULL
    })
    if (is.null(designs)) next
    for (s in c("train", "valid")) {
      parts[[s]]$x <- designs$x[[s]]
      parts[[s]]$z <- designs$z[[s]]
    }
    for (s in c("train", "valid")) result[[paste0("overall.ate.", s)]][i] <- effect(parts[[s]], rep(TRUE, length(parts[[s]]$y)), NA_character_, s, i, "overall", "1")
    for (method in score.method) {
      fitted <- tryCatch({
        tr <- parts$train
        if (method == "twoReg") binary_tworeg_dr(tr, ps.method, minPS, maxPS, ps.control, B, Kfold)$fit else if (method == "logistic") binary_predictions(tr$y, tr$trt, tr$x, "logistic")$fit else binary_boosting(tr, tree.depth, n.trees.boosting, plot.gbmperf)$fit
      }, error = function(e) {
        msg <- paste0(method, ", training iteration ", i, ": score fitting failed: ", conditionMessage(e))
        record(method, "train", i, "score", NA_character_, "error", msg)
        warning(msg, call. = FALSE)
        NULL
      })
      if (is.null(fitted)) next
      for (s in c("train", "valid")) {
        part <- parts[[s]]
        score <- binary_score_new(fitted, part$x, method)
        if (any(!is.finite(score))) { record(method, s, i, "score", NA_character_, "error", "Non-finite score."); warning("Non-finite binary score; returning NA.", call. = FALSE); next }
        rankscore <- if (higher.y) score else -score
        entry <- result[[paste0("ate.", method)]]
        for (j in seq_along(prop)) {
          selected <- rankscore >= unname(stats::quantile(rankscore, 1 - prop[j]))
          entry[[paste0("ate.est.", s, ".high.cv")]][j, i] <- effect(part, selected, method, s, i, "high", as.character(prop[j]))
          if (prop[j] < 1) entry[[paste0("ate.est.", s, ".low.cv")]][j, i] <- effect(part, !selected, method, s, i, "low", as.character(1 - prop[j]))
        }
        cuts <- unname(stats::quantile(rankscore, prop.multi))
        for (j in seq_len(length(prop.multi) - 1L)) {
          selected <- rankscore >= cuts[j] & if (j == length(prop.multi) - 1L) rankscore <= cuts[j + 1L] else rankscore < cuts[j + 1L]
          entry[[paste0("ate.est.", s, ".group.cv")]][j, i] <- effect(part, selected, method, s, i, "group", paste(prop.multi[j], prop.multi[j + 1L], sep = "-"))
        }
        if (s == "valid" && abc) {
          grid <- seq(min(prop), max(prop), length.out = 100)
          vals <- vapply(grid, function(p) effect(part, rankscore >= unname(stats::quantile(rankscore, 1 - p)), method, s, i, "ABC", as.character(p), FALSE), numeric(1))
          ref <- if (higher.y) vals - result$overall.ate.valid[i] else result$overall.ate.valid[i] - vals
          if (all(is.finite(ref))) entry$abc.valid[i] <- auc(grid, ref, type = "spline") else {
            msg <- paste0(method, ", validation iteration ", i, ": ABC unavailable because a required estimate is not finite.")
            record(method, s, i, "ABC", NA_character_, "error", msg); warning(msg, call. = FALSE)
          }
        }
        result[[paste0("ate.", method)]] <- entry
      }
    }
  }
  result$props <- list(prop.onlyhigh = prop, prop.bi = bi, prop.multi = prop.multi)
  result$higher.y <- higher.y; result$abc <- abc; result$cv.n <- cv.n; result$response <- "binary"
  result$formulas <- list(cate.model = cate.model, ps.model = ps.model, trt_labels = d$labels)
  result$diagnostic.messages <- messages; result$subgroup.support <- support
  structure(result, class = "precmed")
}
