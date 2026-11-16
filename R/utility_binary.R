# Author: Thomas Debray <tdebray@fromdatatowisdom.com>
binary_data <- function(data, cate.model, ps.model) {
  cm <- model.frame(cate.model, data, na.action = stats::na.pass)
  pm <- model.frame(ps.model, data, na.action = stats::na.pass)
  y <- model.response(cm); trt <- model.response(pm)
  if (is.null(y) || !(is.numeric(y) || is.logical(y)) || anyNA(y) || !all(y %in% 0:1))
    stop("The binary outcome must be numeric or logical and coded 0/1 without missing values.")
  if (is.null(trt) || anyNA(trt) || length(unique(trt)) != 2L)
    stop("Treatment must contain exactly two non-missing groups.")
  labels <- if (is.logical(trt)) c("FALSE", "TRUE") else if (all(trt %in% 0:1)) c("0", "1") else as.character(sort(unique(trt)))
  t <- as.numeric(as.character(trt) == labels[2])
  x <- model.matrix(cate.model, cm); z <- model.matrix(ps.model, pm)
  if (anyNA(x) || anyNA(z) || any(!is.finite(x)) || any(!is.finite(z))) stop("Missing or non-finite covariates are not allowed.")
  if (ncol(x) == 0L || !"(Intercept)" %in% colnames(x) || !"(Intercept)" %in% colnames(z)) stop("Binary models must include an intercept.")
  if (!is.null(model.offset(cm))) stop("Offsets are not supported for binary outcomes.")
  list(y = as.numeric(y), trt = t, x = x, z = z, labels = labels)
}

# Author: Thomas Debray <tdebray@fromdatatowisdom.com>
binary_predictions <- function(y, trt, x, method = "twoReg", interactions = TRUE) {
  if (length(unique(trt)) != 2L) stop("Both treatment groups are required.")
  # Author: Thomas Debray <tdebray@fromdatatowisdom.com>
  predict_fit <- function(design, index, new) {
    fit <- stats::glm.fit(design[index, , drop = FALSE], y[index], family = stats::binomial())
    if (!fit$converged || fit$rank < ncol(design)) stop("Binary outcome model failed to converge or is rank deficient.")
    list(pred = stats::plogis(drop(new %*% fit$coefficients)), fit = fit)
  }
  if (method == "twoReg") {
    f0 <- predict_fit(x, trt == 0, x); f1 <- predict_fit(x, trt == 1, x)
  } else {
    design <- cbind(x, trt = trt, if (interactions && ncol(x) > 1L) x[, -1, drop = FALSE] * trt)
    d0 <- design; d1 <- design
    d0[, ncol(x) + 1L] <- 0; d1[, ncol(x) + 1L] <- 1
    if (ncol(design) > ncol(x) + 1L) {
      idx <- (ncol(x) + 2L):ncol(design)
      d0[, idx] <- 0; d1[, idx] <- x[, -1, drop = FALSE]
    }
    f0 <- predict_fit(design, rep(TRUE, length(y)), d0)
    f1 <- list(pred = stats::plogis(drop(d1 %*% f0$fit$coefficients)), fit = f0$fit)
  }
  list(mu0 = f0$pred, mu1 = f1$pred, fit = list(control = f0$fit, treated = f1$fit))
}

# Cross-fit nuisance models and project each augmented potential outcome on x.
# Author: Thomas Debray <tdebray@fromdatatowisdom.com>
binary_tworeg_dr <- function(d, ps.method, minPS, maxPS, ps.control, B, Kfold) {
  if (length(B) != 1L || !is.numeric(B) || !is.finite(B) || B < 1 || B != floor(B) ||
      length(Kfold) != 1L || !is.numeric(Kfold) || !is.finite(Kfold) || Kfold < 2 || Kfold != floor(Kfold))
    stop("B must be a positive integer and Kfold must be an integer of at least 2.")
  if (any(tabulate(d$trt + 1L, 2L) < Kfold))
    stop("Binary twoReg requires at least Kfold observations in each treatment group.")
  if (qr(d$x)$rank < ncol(d$x)) stop("Binary twoReg projection design is rank deficient.")
  coefs <- array(NA_real_, c(ncol(d$x), 2L, B))
  for (b in seq_len(B)) {
    fold <- integer(length(d$y))
    for (arm in 0:1) {
      idx <- which(d$trt == arm)
      fold[idx] <- sample(rep(seq_len(Kfold), length.out = length(idx)))
    }
    mu0 <- mu1 <- p <- numeric(length(d$y))
    for (k in seq_len(Kfold)) {
      train <- fold != k; valid <- !train
      for (arm in 0:1) {
        selected <- train & d$trt == arm
        fit <- stats::glm.fit(d$x[selected, , drop = FALSE], d$y[selected], family = stats::binomial())
        if (!fit$converged || fit$rank < ncol(d$x) || any(!is.finite(fit$coefficients)))
          stop("Binary twoReg outcome model failed to converge or is rank deficient in a cross-fitting fold.")
        pred <- stats::plogis(drop(d$x[valid, , drop = FALSE] %*% fit$coefficients))
        if (arm == 0) mu0[valid] <- pred else mu1[valid] <- pred
      }
      args <- list(trt = d$trt[train], x.ps = d$z[train, , drop = FALSE],
                   xnew = d$z[valid, , drop = FALSE], minPS = minPS, maxPS = maxPS)
      if (ps.method == "lasso") args$solver <- resolve_ps_control(ps.control)$solver
      p[valid] <- do.call(if (ps.method == "glm") glm.simplereg.ps else glm.ps, args)
    }
    if (any(!is.finite(p) | p <= 0 | p >= 1)) stop("Propensity scores must be strictly between 0 and 1.")
    pseudo0 <- mu0 + (1 - d$trt) * (d$y - mu0) / (1 - p)
    pseudo1 <- mu1 + d$trt * (d$y - mu1) / p
    coefs[, 1, b] <- stats::lm.fit(d$x, pseudo0)$coefficients
    coefs[, 2, b] <- stats::lm.fit(d$x, pseudo1)$coefficients
  }
  beta <- apply(coefs, c(1, 2), mean)
  dimnames(beta) <- list(colnames(d$x), c("control", "treated"))
  list(score = drop(d$x %*% (beta[, 2] - beta[, 1])),
       fit = list(coefficients = beta, B = B, Kfold = Kfold))
}

# Author: Thomas Debray <tdebray@fromdatatowisdom.com>
binary_effect <- function(d, ps.method, minPS, maxPS, ps.control, interactions = TRUE, diagnostic = FALSE) {
  if (length(unique(d$trt)) != 2L) stop("Both treatment groups are required.")
  p <- if (ps.method == "glm") glm.simplereg.ps(d$trt, d$z, minPS = minPS, maxPS = maxPS, return.raw = diagnostic) else
    glm.ps(d$trt, d$z, minPS = minPS, maxPS = maxPS, solver = resolve_ps_control(ps.control)$solver, return.raw = diagnostic)
  if (any(!is.finite(p) | p <= 0 | p >= 1)) stop("Propensity scores must be strictly between 0 and 1.")
  m <- binary_predictions(d$y, d$trt, d$x, "logistic", interactions)
  r1 <- mean(m$mu1 + d$trt * (d$y - m$mu1) / p)
  r0 <- mean(m$mu0 + (1 - d$trt) * (d$y - m$mu0) / (1 - p))
  result <- c(risk0 = r0, risk1 = r1, risk.difference = r1 - r0)
  if (diagnostic) attr(result, "overlap") <- treatment_overlap(d$trt, p, minPS, maxPS)
  result
}

# Author: Thomas Debray <tdebray@fromdatatowisdom.com>
binary_subset <- function(d, i) {
  list(y = d$y[i], trt = d$trt[i], x = d$x[i, , drop = FALSE], z = d$z[i, , drop = FALSE])
}

