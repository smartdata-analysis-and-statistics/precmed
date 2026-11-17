#' Plot estimated conditional treatment-effect scores
#'
#' Display the distribution of scores fitted by \code{catefit()}. These are
#' in-sample predictions; their spread alone does not establish treatment-effect
#' heterogeneity or validate the score.
#'
#' @param x An object of class \code{"catefit"}.
#' @param methods Optional character vector of fitted scoring methods to show,
#'   in the desired order. The default shows all available scores.
#' @param type Either \code{"histogram"} (default) or \code{"density"}.
#'   Density curves require at least two distinct finite values per method.
#' @param bins Positive integer number of histogram bins; default 30.
#' @param null.value Position of a dashed reference line on the score scale;
#'   default 0. Use \code{NULL} or \code{NA} to omit the line.
#' @param theme A ggplot2 theme, defaulting to the package plot theme.
#' @param ... Currently unused.
#' @return A ggplot object. Non-finite scores are omitted with a warning.
#' @details Binary scores are risk differences, count scores are log rate ratios,
#'   survival scores are log restricted mean time lost ratios, and continuous
#'   scores are mean differences. The treatment contrast is treatment 1 versus
#'   treatment 0 (or treatment 1 minus treatment 0 for differences). The
#'   methods may impose different model assumptions; see \code{\link{catefit}()}.
#'   Original treatment labels are shown when stored in the fit; older objects
#'   without labels use treatment 1 and treatment 0.
#'   Histogram panels show patient counts for each method. Density curves are
#'   normalized separately and overlaid for comparison.
#' @examples
#' \donttest{
#' d <- simulate_binary(n = 300, seed = 1)
#' fit <- catefit("binary", d, y ~ age, "logistic", "two_regression", trt ~ age,
#'                prop.cutoff = 1, B = 1, Kfold = 2, seed = 1, higher.y = FALSE)
#' plot(fit)
#' }
#' @author Thomas Debray \email{tdebray@fromdatatowisdom.com}
#' @export
plot.catefit <- function(x, methods = NULL, type = c("histogram", "density"),
                         bins = 30, null.value = 0, theme = precmed_plot_theme(), ...) {
  type <- match.arg(type)
  available <- sub("^score\\.", "", setdiff(names(x)[startsWith(names(x), "score.")], "score.method"))
  if (!length(available)) stop("The catefit object has no scores to plot.")
  if (is.null(methods)) methods <- available
  if (!is.character(methods) || !length(methods) || anyNA(methods) ||
      anyDuplicated(methods) || any(!methods %in% available))
    stop("methods must contain distinct available scoring methods: ", paste(available, collapse = ", "), ".")
  if (!is.numeric(bins) || length(bins) != 1L || !is.finite(bins) || bins < 1 || bins != floor(bins))
    stop("bins must be a positive integer.")
  if (!is.null(null.value) && (length(null.value) != 1L ||
      !(is.numeric(null.value) || (is.logical(null.value) && is.na(null.value))) ||
      (!is.na(null.value) && !is.finite(null.value))))
    stop("null.value must be a finite number, NA, or NULL.")
  scores <- lapply(methods, function(method) x[[paste0("score.", method)]])
  if (any(vapply(scores, function(z) !is.numeric(z) || !length(z), logical(1))))
    stop("Scores must be non-empty numeric vectors.")
  if (any(vapply(scores, function(z) any(!is.finite(z)), logical(1))))
    warning("Non-finite scores were omitted.", call. = FALSE)
  scores <- lapply(scores, function(z) z[is.finite(z)])
  if (any(lengths(scores) == 0L)) stop("Every selected method must have at least one finite score.")
  if (type == "density" && any(vapply(scores, function(z) length(unique(z)) < 2L, logical(1))))
    stop("Density plots require at least two distinct finite scores for each method; use type = 'histogram'.")
  dat <- data.frame(score = unlist(scores, use.names = FALSE),
                    method = factor(rep(methods, lengths(scores)), levels = methods))
  labels <- x$treatment.labels
  if (is.null(labels) || length(labels) != 2L || anyNA(labels) || any(!nzchar(as.character(labels))) ||
      identical(as.character(labels), c("0", "1"))) labels <- c("treatment 0", "treatment 1")
  labels <- as.character(labels)
  difference <- paste(labels[2], "-", labels[1])
  ratio <- paste(labels[2], "vs", labels[1])
  xlabel <- switch(x$response, binary = paste0("Estimated risk difference (", difference, ")"),
                   count = paste0("Estimated log rate ratio (", ratio, ")"),
                   survival = paste0("Estimated log RMTL ratio (", ratio, ")"),
                   continuous = paste0("Estimated mean difference (", difference, ")"),
                   "Estimated CATE score")
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$score))
  if (type == "histogram") {
    p <- p + ggplot2::geom_histogram(bins = bins, fill = "#7293AD", color = "white") +
      ggplot2::facet_wrap(ggplot2::vars(.data$method)) + ggplot2::labs(y = "Patients")
  } else {
    p <- p + ggplot2::geom_density(ggplot2::aes(color = .data$method, fill = .data$method), alpha = 0.25) +
      ggplot2::labs(y = "Density", color = "Method", fill = "Method")
  }
  if (!is.null(null.value) && !is.na(null.value))
    p <- p + ggplot2::geom_vline(xintercept = null.value, linetype = "dashed", color = "#52616B")
  p + ggplot2::labs(x = xlabel) + theme
}
