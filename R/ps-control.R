# PS interaction screening options are distinct from fitting method and GBM ... .
resolve_ps_control <- function(ps.control) {
  if (!is.list(ps.control) || (length(ps.control) &&
      (is.null(names(ps.control)) || anyNA(names(ps.control)) ||
       any(!nzchar(names(ps.control))) || anyDuplicated(names(ps.control)) ||
       any(!names(ps.control) %in% "solver")))) {
    stop("ps.control must be a named list containing solver.")
  }
  list(solver = match.arg(if (length(ps.control)) ps.control$solver else "C", c("C", "R")))
}

# Materialize only selected interactions. Both paths preserve pair ordering,
# sample variances, strict SMD > .1, and the original variance stabilizer.
ps_interactions <- function(x, trt, xnew = NULL, solver = c("C", "R")) {
  solver <- match.arg(solver)
  x <- as.matrix(x)
  if (!is.numeric(x) || any(!is.finite(x)) || !is.numeric(trt) ||
      length(trt) != nrow(x) || anyNA(trt) || any(!trt %in% c(0, 1))) {
    stop("PS screening requires finite numeric predictors and binary treatment.")
  }
  if (!is.null(xnew) && (!is.matrix(xnew) || !is.numeric(xnew) ||
      ncol(xnew) != ncol(x) || any(!is.finite(xnew)))) {
    stop("PS prediction design must be finite and have the same columns as training data.")
  }
  p <- ncol(x)
  if (p > 1 && any(tabulate(trt + 1L, nbins = 2) < 2))
    stop("Interaction screening requires at least two observations per treatment group.")
  select_pair <- function(v) {
    a <- v[trt == 1]; b <- v[trt == 0]
    variance <- stats::var(a) + stats::var(b)
    sqrt(2) * abs(mean(a) - mean(b)) / sqrt(variance + 1e-6) > .1 && variance > 0
  }
  if (p < 2) {
    pairs <- matrix(integer(), 2, 0)
  } else if (solver == "C") {
    storage.mode(x) <- "double"
    screened <- .Call(precmed_ps_screen, x, as.double(trt))
    pairs <- screened[[1L]]
    keep <- rep(TRUE, ncol(pairs))
    for (k in which(screened[[2L]])) {
      keep[k] <- select_pair(x[, pairs[1L, k]] * x[, pairs[2L, k]])
    }
    pairs <- pairs[, keep, drop = FALSE]
  } else {
    # Retain the original dense implementation as the R reference backend.
    pairs <- utils::combn(p, 2)
    all <- matrix(NA_real_, nrow(x), ncol(pairs))
    for (k in seq_len(ncol(pairs))) all[, k] <- x[, pairs[1L, k]] * x[, pairs[2L, k]]
    if (any(!is.finite(all))) stop("PS interaction products must be finite.")
    mu1 <- apply(all[trt == 1, , drop = FALSE], 2, mean)
    mu0 <- apply(all[trt == 0, , drop = FALSE], 2, mean)
    v1 <- apply(all[trt == 1, , drop = FALSE], 2, stats::var)
    v0 <- apply(all[trt == 0, , drop = FALSE], 2, stats::var)
    smd <- sqrt(2) * abs(mu1 - mu0) / sqrt(v1 + v0 + 1e-6)
    pairs <- pairs[, smd > .1 & v1 + v0 > 0, drop = FALSE]
  }
  train <- matrix(NA_real_, nrow(x), ncol(pairs))
  prediction <- if (!is.null(xnew)) matrix(NA_real_, nrow(xnew), ncol(pairs))
  for (k in seq_len(ncol(pairs))) {
    train[, k] <- x[, pairs[1L, k]] * x[, pairs[2L, k]]
    if (!is.null(xnew)) prediction[, k] <- xnew[, pairs[1L, k]] * xnew[, pairs[2L, k]]
  }
  list(train = train, prediction = prediction, pairs = pairs)
}
