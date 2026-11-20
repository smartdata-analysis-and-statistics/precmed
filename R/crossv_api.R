# Rematch the previous positional argument order before resolving named aliases.
# Author: Thomas Debray <tdebray@fromdatatowisdom.com>
match_catecv_call <- function(fun, call, caller, previous) {
  matched <- as.list(match.call(fun, call, envir = caller))[-1L]
  supplied <- names(as.list(call)[-1L])
  modern <- c("outcome_model", "outcome_method", "benefit_model", "benefit_method")
  old_named <- any(c("cate.model", "score.method") %in% supplied) &&
    !any(modern %in% supplied)
  candidate <- if (!"outcome_model" %in% supplied && "outcome_model" %in% names(matched))
    tryCatch(eval(matched$outcome_model, caller), error = function(e) NULL)
  old_positional <- is.character(candidate) && !any(grepl("~", candidate, fixed = TRUE))
  if (old_named || old_positional)
    return(as.list(match.call(previous, call, envir = caller))[-1L])
  matched
}

# Translate public score configurations to the shared fitting routines.
# Author: Thomas Debray <tdebray@fromdatatowisdom.com>
resolve_catecv_call <- function(args, response, caller) {
  args <- resolve_cate_aliases(args, censoring = response == "survival", caller = "catecv()")
  old <- c(cate.model = "outcome_model", score.method = "benefit_method",
           initial.predictor.method = "outcome_method", init.model = "benefit_model")
  for (name in intersect(names(old), names(args))) {
    replacement <- old[[name]]
    if (sum(names(args) == name) > 1L || replacement %in% names(args))
      stop("Specify only ", replacement, "; ", name, " is deprecated.", call. = FALSE)
    if (name == "init.model") {
      warning(name, " is deprecated for catecv(); use outcome_model for the outcome formula and benefit_model for a separate score formula.", call. = FALSE)
    } else if (name == "score.method") {
      warning(name, " is deprecated for catecv(); use outcome_method and benefit_method instead.", call. = FALSE)
    } else {
      warning(name, " is deprecated for catecv(); use ", replacement, " instead.", call. = FALSE)
    }
  }

  initial <- if ("initial.predictor.method" %in% names(args))
    eval(args[["initial.predictor.method"]], caller)
  if (response == "survival") {
    initial <- normalize_survival_initial_method(initial)
    if (identical(initial, "random_forest")) initial <- "randomForest"
  }
  if (response == "binary" && !is.null(initial) && !identical(initial, "logistic"))
    stop("Binary two-regression uses logistic outcome models.", call. = FALSE)

  legacy_score <- "score.method" %in% names(args)
  if (legacy_score) {
    scores <- eval(args[["score.method"]], caller)
    check_score_method(response, scores)
    if ("outcome_method" %in% names(args))
      stop("Specify benefit_method to pair scoring methods with outcome_method.", call. = FALSE)
    config <- list(scores = scores, initial = initial)
  } else {
    methods <- if ("benefit_method" %in% names(args)) eval(args[["benefit_method"]], caller)
    outcome <- if ("outcome_method" %in% names(args)) eval(args[["outcome_method"]], caller)
    if (response == "continuous" && (is.null(methods) || identical(methods, "t_learner")))
      stop("Continuous-outcome cross-validation is not implemented.", call. = FALSE)
    config <- resolve_catefit_arguments(response, TRUE, outcome, methods, NULL, NULL,
                                        initial, !"initial.predictor.method" %in% names(args))
  }

  names(args)[names(args) == "outcome_model"] <- "cate.model"
  names(args)[names(args) == "ps_model"] <- "ps.model"
  names(args)[names(args) == "ps_method"] <- "ps.method"
  names(args)[names(args) == "censoring_model"] <- "ipcw.model"
  names(args)[names(args) == "censoring_method"] <- "ipcw.method"
  args[c("outcome_method", "benefit_method", "initial.predictor.method")] <- NULL
  if (!"cate.model" %in% names(args)) stop("outcome_model must be supplied.", call. = FALSE)
  if (!"ps.model" %in% names(args)) stop("ps_model must be supplied.", call. = FALSE)
  args[["score.method"]] <- config$scores
  if (!is.null(config$initial) && response != "binary")
    args[["initial.predictor.method"]] <- config$initial
  if (length(config$initial_methods) && response %in% c("count", "survival"))
    args[["initial.methods"]] <- config$initial_methods
  args
}

# Construct the benefit design with transformations estimated from training rows.
# Author: Thomas Debray <tdebray@fromdatatowisdom.com>
cv_benefit_design <- function(data, valid, formula, x_train, x_valid) {
  if (is.null(formula)) return(list(train = x_train, valid = x_valid))
  train <- !seq_len(nrow(data)) %in% valid
  design <- binary_cv_design(data, train, formula, benefit = TRUE, valid_rows = valid)
  lapply(design, function(x) x[, colnames(x) != "(Intercept)", drop = FALSE])
}
