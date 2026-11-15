# Options forwarded to randomForestSRC for both treatment-arm forests.
resolve_rf_control <- function(rf.control) {
  allowed <- c("mtry", "nodesize", "nsplit", "nodedepth", "splitrule")
  if (!is.list(rf.control) || (length(rf.control) &&
      (is.null(names(rf.control)) || anyNA(names(rf.control)) ||
       any(!nzchar(names(rf.control))) || anyDuplicated(names(rf.control)) ||
       any(!names(rf.control) %in% allowed))))
    stop("rf.control must be a named list containing only mtry, nodesize, nsplit, nodedepth, or splitrule; use n.trees.rf for the number of trees.")
  for (name in intersect(names(rf.control), c("mtry", "nodesize", "nsplit", "nodedepth"))) {
    x <- rf.control[[name]]
    minimum <- if (name == "nsplit") 0 else 1
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x < minimum || x != floor(x))
      stop("rf.control$", name, " must be an integer >= ", minimum, ".")
  }
  if ("splitrule" %in% names(rf.control) &&
      (!is.character(rf.control$splitrule) || length(rf.control$splitrule) != 1L ||
       is.na(rf.control$splitrule) || !nzchar(rf.control$splitrule)))
    stop("rf.control$splitrule must be a non-empty character string.")
  rf.control
}

fit_rf <- function(formula, data, n.trees.rf, rf.control) {
  do.call(randomForestSRC::rfsrc, c(list(formula = formula, data = data, ntree = n.trees.rf), rf.control))
}
