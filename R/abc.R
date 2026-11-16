#' Compute the area between curves from the \code{"precmed"} object
#'
#' Compute the area between curves (ABC) for each scoring method in the \code{"precmed"} object.
#' This should be run only after results of \code{\link{catecv}()} have been obtained.
#'
#' @param x An object of class \code{"precmed"}.
#' @param ... Additional arguments (currently unused).
#'
#'
#' @return Returns a matrix of numeric values with number of columns equal to the number cross-validation
#' iteration and number of rows equal to the number of scoring methods in \code{x}.
#'
#' @details The ABC is the area between a validation curve and the overall ATE in the validation set.
#' It is calculated for each scoring method separately. Higher ABC values are preferable as they
#' indicate that more treatment effect heterogeneity is captured by the scoring method.
#' Negative values of ABC are possible if segments of the validation curve cross the overall ATE line.
#' The ABC is calculated with the \code{\link{auc}()} in \code{utility.R} with a natural
#' cubic spline interpolation. The calculation of the ABC is always based on validation curves based on
#' 100 proportions equally spaced from \code{min(prop.cutoff)} to \code{max(prop.cutoff)}.
#'
#' The ABC is a metric to help users select the best scoring method in terms of capturing treatment
#' effect heterogeneity in the data. It should be used in complement to the visual inspection of
#' the validation curves in the validation set in \code{\link{plot}()}.
#'
#' @references Zhao, L., Tian, L., Cai, T., Claggett, B., & Wei, L. J. (2013).
#' \emph{Effectively selecting a target population for a future comparative study.
#' Journal of the American Statistical Association, 108(502), 527-539.}
#'
#' Jiang, X., Simoneau, G., Pellegrini, F., Shen, C., de Moor, C., & Tian, L. (2026).
#' Visualization and interpretation of individualized treatment rule results. In
#' T. P. A. Debray, T.-L. Nguyen, & R. W. Platt (Eds.),
#' \emph{Comparative Effectiveness and Personalized Medicine Research Using Real-World Data}
#' (Chapter 20, pp. 513--536). Chapman & Hall/CRC.
#' \doi{10.1201/9781003300809-20}
#'
#' @seealso \code{\link{catecv}()} function and \code{\link{plot}()}, \code{\link{boxplot}()} methods for
#' \code{"precmed"} objects.
#'
#' @examples
#' # The result was created with:
#' # catecvCountExample <- catecv(
#' #   response = "count", data = countExample, score.method = "poisson",
#' #   cate.model = y ~ age + female + previous_treatment + previous_cost +
#' #     previous_number_relapses + offset(log(years)),
#' #   ps.model = trt ~ age + previous_treatment, higher.y = FALSE,
#' #   cv.n = 5, seed = 999, verbose = 0)
#' # Load the pre-fitted result instead of repeating the calculation.
#' data(catecvCountExample)
#' # Rows are scoring methods; columns are validation iterations.
#' abc(catecvCountExample)
#'
#' @export
abc <- function(x, ...) UseMethod("abc", x)

#' @export
#' @rdname abc
abc.default <- function(x, ...){
  warning(paste("ABC does not know how to handle object of class ",
                class(x),
                "and can only be used on class precmed"))
}

#' Compute the area between curves from the \code{"precmed"} object
#'
#' Compute the area between curves (ABC) for each scoring method in the \code{"precmed"} object.
#' This should be run only after results of \code{\link{catecv}()} have been obtained.
#'
#' @param x An object of class \code{"precmed"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Returns a matrix of numeric values with number of columns equal to the number cross-validation
#' iteration and number of rows equal to the number of scoring methods in \code{x}.
#'
#' @details The ABC is the area between a validation curve and the overall ATE in the validation set.
#' It is calculated for each scoring method separately. Higher ABC values are preferable as they
#' indicate that more treatment effect heterogeneity is captured by the scoring method.
#' Negative values of ABC are possible if segments of the validation curve cross the overall ATE line.
#' The ABC is calculated with the \code{\link{auc}()} in \code{utility.R} with a natural
#' cubic spline interpolation. The calculation of the ABC is always based on validation curves based on
#' 100 proportions equally spaced from \code{min(prop.cutoff)} to \code{max(prop.cutoff)}.
#'
#' The ABC is a metric to help users select the best scoring method in terms of capturing treatment
#' effect heterogeneity in the data. It should be used in complement to the visual inspection of
#' the validation curves in the validation set in \code{\link{plot}()}.
#'
#' @references Zhao, L., Tian, L., Cai, T., Claggett, B., & Wei, L. J. (2013).
#' \emph{Effectively selecting a target population for a future comparative study.
#' Journal of the American Statistical Association, 108(502), 527-539.}
#'
#' Jiang, X., Simoneau, G., Pellegrini, F., Shen, C., de Moor, C., & Tian, L. (2026).
#' Visualization and interpretation of individualized treatment rule results. In
#' T. P. A. Debray, T.-L. Nguyen, & R. W. Platt (Eds.),
#' \emph{Comparative Effectiveness and Personalized Medicine Research Using Real-World Data}
#' (Chapter 20, pp. 513--536). Chapman & Hall/CRC.
#' \doi{10.1201/9781003300809-20}
#'
#' @seealso \code{\link{catecv}()} function and \code{\link{plot}()}, \code{\link{boxplot}()} methods for
#' \code{"precmed"} objects.
#'
#' @examples
#' # The result was created with:
#' # catecvCountExample <- catecv(
#' #   response = "count", data = countExample, score.method = "poisson",
#' #   cate.model = y ~ age + female + previous_treatment + previous_cost +
#' #     previous_number_relapses + offset(log(years)),
#' #   ps.model = trt ~ age + previous_treatment, higher.y = FALSE,
#' #   cv.n = 5, seed = 999, verbose = 0)
#' # Load the pre-fitted result instead of repeating the calculation.
#' data(catecvCountExample)
#' # Rows are scoring methods; columns are validation iterations.
#' abc(catecvCountExample)
#'
#' @export
#'
#' @importFrom stringr str_extract
abc.precmed <- function(x, ...) {

  # Check the value of abc (must have area between curves values)
  if (!x$abc) {
    stop("Area between curves (ABC) must have been calculated in x by setting abc = TRUE in catecv().")
  }

  score.method <- str_extract(names(x), "(?<=^ate\\.).*$")
  score.method <- score.method[!is.na(score.method)]
  m <- length(score.method)

  # Value of ABC
  value.abc <- matrix(NA, nrow = m, ncol = x$cv.n)
  for (i in seq(m)) {
    value.abc[i,] <- x[[paste0("ate.", score.method[i])]]$abc.valid
  }
  rownames(value.abc) <- score.method
  colnames(value.abc) <- paste0("cv", 1:(x$cv.n))

  return(value.abc)
}

