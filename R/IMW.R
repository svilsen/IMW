#' @title Incremental Moment Windows
#' 
#' @description Incrementally updated moments for a fixed windows. 
#' 
#' @param x A vector, or matrix with a single column, of new data.
#' @param k A numeric indicating the lag/windows.
#' 
#' @return An object of class \code{imw}, containing the following:
#' \describe{
#'     \item{\code{stats}}{A matrix whose columns contain the rolling mean, variance, skewness, and kurtosis.}
#'     \item{\code{totalMoments}}{The incremental moments calculated using all the supplied information.}
#'     \item{\code{windowMoments}}{The incremental moments calculated using the supplied information up until the beginning of the last window (need if more information is added later, using the \link{update_imw} function).}
#'     \item{\code{x}}{A matrix of the data supplied to the function (gets updated when using the \link{update_imw} function).}
#'     \item{\code{k}}{The size of the window used to calculate the statistics.}
#' }
#'  
#' @examples x <- cumsum(rnorm(100))
#' k <- 10
#' 
#' imw(x, k)
#'  
#' @export 
imw <- function(x, k) {
    if (!is.matrix(x)) {
        x <- matrix(x, ncol = 1)
    }
    
    if (!is.integer(k)) {
        k <- as.integer(k)
    }
    
    res <- imw_cpp(x, k)
    res$x <- x
    res$k <- k
    
    class(res) <- "imw"
    return(res)
}


#' @export
mean.imw <- function(x, ...) {
    res <- x$stats[, 1]
    return(res)
}

#' @export
variance <- function(x, ...) {
    UseMethod("variance")
}

#' @export
variance.imw <- function(x, ...) {
    dots <- list(...)
    if (is.null(dots$type)) {
        type <- 1
    }
    else {
        type <- dots$type
    }
    
    if (type == 1) {
        res <- x$stats[, 2]
    }
    else if (type == 2) {
        k <- dim(x$k)[1]
        res <- (k - 1) * x$stats[, 2] / k
    }
    
    return(res)
}

#' @export
skewness <- function(x, ...) {
    UseMethod("skewness")
}

#' @export
skewness.imw <- function(x, ...) {
    dots <- list(...)
    if (is.null(dots$type)) {
        type <- 1
    }
    else {
        type <- dots$type
    }
    
    if (type == 1) {
        res <- x$stats[, 3]
    }
    else if (type == 2) {
        k <- x$k
        if (k < 3) {
            stop("The window size 'k' needs to be at least 3.")
        }
        
        res <- x$stats[, 3] * sqrt(k * (k - 1)) / (k - 2)
    }
    else if (type == 3) {
        k <- x$k
        res <- x$stats[, 3] * (1.0 - 1.0 / k)^(3 / 2)
    }
    else {
        stop("The supplied 'type' was not valid.")
    }
    
    
    return(res)
}

#' @export
kurtosis <- function(x, ...) {
    UseMethod("kurtosis")
}

#' @export
kurtosis.imw <- function(x, ...) {
    dots <- list(...)
    if (is.null(dots$type)) {
        type <- 1
    }
    else {
        type <- dots$type
    }
    
    if (type == 1) {
        res <- x$stats[, 4]
    }
    else if (type == 2) {
        k <- x$k
        
        if (k < 4) {
            stop("The window size 'k' needs to be at least 4.")
        }
        
        res <- ((k + 1) * x$stats[, 4] + 6.0) * (k - 1) / ((k - 2) * (k - 3))
        
    }
    else if (type == 3) {
        k <- x$k
        res <- (x$stats[, 4] + 3) * (1.0 + 1.0 / k)^2 - 3.0
    }
    else {
        stop("The supplied 'type' was not valid.")
    }
    
    return(res)
}

#' @export
as.matrix.imw <- function(x, ...) {
    res <- x$stats
    colnames(res) <- c("Mean", "Variance", "Skewness", "Kurtosis")
    return(res)
}

#' @title Update
#' 
#' @description Update an object of class \link{imw} with new information.
#' 
#' @param object An object of class \code{imw}.
#' @param x_new A vector, or single column matrix, of new data.
#' 
#' @return An object of class \code{imw}.
#' 
#' @examples N <- 100
#' N_new <- 20
#' k <- 10
#' 
#' y <- cumsum(rnorm(N + N_new))
#' x <- head(y, N)
#' x_new <- tail(y, N_new)
#' 
#' imw_x <- imw(x, k)
#' 
#' update_imw(imw_x, x_new = x_new)
#' 
#' @export 
update_imw <- function(object, x_new) {
    if (!is.matrix(x_new)) {
        x_new <- matrix(x_new, ncol = 1)
    }
    
    N <- length(object$x)
    x <- object$x
    k <- object$k
    
    x_c <- rbind(tail(x, k), x_new)
    res <- imw_update_cpp(x = x_c, k = k, t = object$totalMoments, l = object$windowMoments)
    
    res$stats <- rbind(object$stats, res$stats)
    res$x <- rbind(x, x_new)
    res$k <- k
    
    class(res) <- "imw"
    return(res)
}

