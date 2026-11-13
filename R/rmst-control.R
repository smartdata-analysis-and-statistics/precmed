# Keep RMST integration options separate from model-fitting arguments.
resolve_rmst_control <- function(rmst.control) {
  if (!is.list(rmst.control) || (length(rmst.control) &&
      (is.null(names(rmst.control)) || anyNA(names(rmst.control)) ||
       any(!nzchar(names(rmst.control))) || anyDuplicated(names(rmst.control)) ||
       any(!names(rmst.control) %in% "solver")))) {
    stop("rmst.control must be a named list containing solver.")
  }
  solver <- if (length(rmst.control)) rmst.control$solver else "C"
  list(solver = match.arg(solver, c("C", "R")))
}
