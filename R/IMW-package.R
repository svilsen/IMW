#' @keywords internal 
"_PACKAGE"

#' @title Incremental Moment Windows
#'  
#' @encoding UTF-8
#' 
#' @importFrom Rcpp evalCpp
#' @importFrom utils tail
#' @importFrom graphics plot points abline par
#' 
#' @useDynLib IMW
#' 
#' @name IMW-package
NULL

#' @title An imw-object
#'
#' @description An imw-object is a list of class \code{imw} containing the following:
#' \describe{
#'     \item{\code{stats}}{Matrix: The rolling mean, variance, skewness, and kurtosis of data.}
#'     \item{\code{lagged}}{Vector: A vector of the incremental moments, as well as the number of observations used to calculate these moments, at the beginning of the last window.}
#'     \item{\code{total}}{Vector: A vector of the incremental moments, as well as the number of observations used to calculate these moments, at the end of the last window.}
#'     \item{\code{x}}{Matrix: The data used to estimate the moments.}
#'     \item{\code{k}}{Numeric: The size of the window used to calculate the four moments.}
#' }
#'
#' @name imw-object
#' @rdname imw-object
NULL