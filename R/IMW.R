#' @title Incremental Moment Windows
#' 
#' @description Incrementally updated moments for a fixed window size. 
#' 
#' @param x A vector of data.
#' @param k A numeric defining the window size.
#' 
#' @return An object of class \code{imw}, containing the following:
#' \describe{
#'     \item{\code{stats}}{A matrix which columns contain the rolling mean, variance, skewness, and kurtosis.}
#'     \item{\code{total}}{The incremental moments calculated using the supplied information.}
#'     \item{\code{window}}{The incremental moments calculated using the supplied information up until the beginning of the last window (need if more information is recieved and the moments needs to be updated, see \link{uimw} function).}
#'     \item{\code{x}}{A matrix of the supplied data (gets replaced when using the \link{uimw} function).}
#'     \item{\code{k}}{The size of the window used to calculate the four moments.}
#' }
#'  
#' @examples N <- 100
#' x <- cumsum(rnorm(N))
#' 
#' k <- 10
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

#' @title Mean
#' 
#' @description Computes the mean.
#' 
#' @param x An \link{imw} object.
#' @param ... Additional arguments; not used in this instance.
#' 
#' @return The estimated mean.
#' 
#' @rdname mean.imw
#' @method mean imw
#' 
#' @examples N <- 100
#' x <- cumsum(rnorm(N))
#' 
#' k <- 10
#' mw <- imw(x, k)
#' 
#' mean(mw)
#' 
#' @export
mean.imw <- function(x, ...) {
    res <- x$stats[, 1]
    return(res)
}

#' @title Variance
#' 
#' @description Computes the variance.
#' 
#' @param x An \link{imw} object.
#' @param ... Additional arguments, see details.
#' 
#' @details The only additional argument is \code{type} which takes the values \code{1} and \code{2} corresponding to the population and sample variance, respectively.  
#' 
#' @return The estimated variance.
#' 
#' @examples N <- 100
#' x <- cumsum(rnorm(N))
#' 
#' k <- 10
#' mw <- imw(x, k)
#' 
#' variance(mw)
#' 
#' @export
variance <- function(x, ...) {
    UseMethod("variance")
}

#' @title Variance
#' 
#' @description Computes the variance.
#' 
#' @param x An \link{imw} object.
#' @param ... Additional arguments, see details.
#' 
#' @details The only additional argument is \code{type} which takes the values \code{1} and \code{2} corresponding to the population and sample variance, respectively.  
#' 
#' @return The estimated variance.
#' 
#' @rdname variance.imw
#' @method variance imw
#' 
#' @examples N <- 100
#' x <- cumsum(rnorm(N))
#' 
#' k <- 10
#' mw <- imw(x, k)
#' 
#' variance(mw)
#' 
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

#' @title Skewness
#' 
#' @description Computes the skewness.
#' 
#' @param x An \link{imw} object.
#' @param ... Additional arguments, see details.
#' 
#' @details The only additional argument is \code{type} which takes the values \code{1, 2, 3}, corresponding to the three methods discussed by Joanes and Gill (1998).  
#' 
#' @return The estimated skewness.
#' 
#' @references D. N. Joanes and C. A. Gill (1998), Comparing measures of sample skewness and kurtosis. The Statistician, 47, 183–189.
#' 
#' @examples N <- 100
#' x <- cumsum(rnorm(N))
#' 
#' k <- 10
#' mw <- imw(x, k)
#' 
#' skewness(mw)
#' 
#' @export
skewness <- function(x, ...) {
    UseMethod("skewness")
}

#' @title Skewness
#' 
#' @description Computes the skewness.
#' 
#' @param x An \link{imw} object.
#' @param ... Additional arguments, see details.
#' 
#' @details The only additional argument is \code{type} which takes the values \code{1, 2, 3}, corresponding to the three methods discussed by Joanes and Gill (1998).  
#' 
#' @return The estimated skewness.
#' 
#' @references D. N. Joanes and C. A. Gill (1998), Comparing measures of sample skewness and kurtosis. The Statistician, 47, 183–189.
#' 
#' @rdname skewness.imw
#' @method skewness imw
#' 
#' @examples N <- 100
#' x <- cumsum(rnorm(N))
#' 
#' k <- 10
#' mw <- imw(x, k)
#' 
#' skewness(mw)
#' 
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

#' @title Kurtosis
#' 
#' @description Computes the kurtosis.
#' 
#' @param x An \link{imw} object.
#' @param ... Additional arguments, see details.
#' 
#' @details The only additional argument is \code{type} which takes the values \code{1, 2, 3}, corresponding to the three methods discussed by Joanes and Gill (1998).  
#' 
#' @return The estimated kurtosis.
#' 
#' @references D. N. Joanes and C. A. Gill (1998), Comparing measures of sample skewness and kurtosis. The Statistician, 47, 183–189.
#' 
#' @examples N <- 100
#' x <- cumsum(rnorm(N))
#' 
#' k <- 10
#' mw <- imw(x, k)
#' 
#' kurtosis(mw)
#' 
#' @export
kurtosis <- function(x, ...) {
    UseMethod("kurtosis")
}

#' @title Kurtosis
#' 
#' @description Computes the kurtosis.
#' 
#' @param x An \link{imw} object.
#' @param ... Additional arguments, see details.
#' 
#' @details The only additional argument is \code{type} which takes the values \code{1, 2, 3}, corresponding to the three methods discussed by Joanes and Gill (1998).  
#' 
#' @return The estimated kurtosis.
#' 
#' @references D. N. Joanes and C. A. Gill (1998), Comparing measures of sample skewness and kurtosis. The Statistician, 47, 183–189.
#' 
#' @rdname kurtosis.imw
#' @method kurtosis imw
#' 
#' @examples N <- 100
#' x <- cumsum(rnorm(N))
#' 
#' k <- 10
#' mw <- imw(x, k)
#' 
#' kurtosis(mw)
#'
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

#' @title Converts to matrix
#' 
#' @description Converts an object of \link{imw} to matrix.
#' 
#' @param x An \link{imw} object.
#' @param ... Additional arguments; not used in this instance.
#' 
#' @rdname as.matrix.imw
#' @method as.matrix imw
#' 
#' @return A matrix.
#' 
#' @examples N <- 100
#' x <- cumsum(rnorm(N))
#' 
#' k <- 10
#' mw <- imw(x, k)
#' 
#' as.matrix(mw)
#'
#' @export
as.matrix.imw <- function(x, ...) {
    res <- x$stats
    colnames(res) <- c("Mean", "Variance", "Skewness", "Kurtosis")
    return(res)
}

#' @title Update Incremental Moment Windows
#' 
#' @description Update an object of class \link{imw} with new information.
#' 
#' @param object An object of class \code{imw}.
#' @param x_new A vector, or column matrix, of new information.
#' 
#' @return An object of class \code{imw}.
#' 
#' @examples # Initial data
#' N <- 100
#' x <- cumsum(rnorm(N))
#' 
#' k <- 10
#' mw <- imw(x, k)
#' 
#' # Additional data
#' M <- 20
#' y <- cumsum(c(tail(x, 1), rnorm(M)))[-1]
#' 
#' uimw(mw, y)
#' 
#' @export 
uimw <- function(object, x_new) {
    if (!is.matrix(x_new)) {
        x_new <- matrix(x_new, ncol = 1)
    }
    
    N <- length(object$x)
    x <- object$x
    k <- object$k
    
    x_c <- rbind(tail(x, k), x_new)
    res <- uimw_cpp(x = x_c, k = k, t = object$totalMoments, l = object$windowMoments)
    
    res$stats <- rbind(object$stats, res$stats)
    res$x <- rbind(x, x_new)
    res$k <- k
    
    class(res) <- "imw"
    return(res)
}

#' @title Plot of rolling statistics
#' 
#' @param x An imw-object.
#' 
#' @param ... Additional arguments passed to the internal plot functions.
#' 
#' @return NULL
#' 
#' @rdname plot.imw
#' @method plot imw
#' 
#' @examples N <- 100
#' x <- cumsum(rnorm(N))
#' 
#' k <- 10
#' mw <- imw(x, k)
#' 
#' plot(mw)
#'
#' @export
plot.imw <- function(x, ...) {
    par(mfrow = c(2, 2), mar = c(5, 5, 0.5, 0.5))
    plot(x$x, pch = 16, xlab = "t", ylab = expression("x"[" t"]), ...)
    points(mean(x), pch = 16, col = "red", type = "l", lwd = 3)
    
    plot(variance(x), pch = 16, xlab = "t", ylab = expression("Variance(x"[" t "]*")"), ...)
    abline(h = mean(variance(x), na.rm = TRUE), col = "red", lwd = 3)
    
    plot(skewness(x), pch = 16, xlab = "t", ylab = expression("Skewness(x"[" t "]*")"), ...)
    abline(h = 0, col = "red", lwd = 3)
    
    plot(kurtosis(x), pch = 16, xlab = "t", ylab = expression("Kurtosis(x"[" t "]*")"), ...)
    abline(h = 0, col = "red", lwd = 3)
    par(mfrow = c(1, 1), mar = c(6, 6, 2, 2))
    
    return(invisible(NULL))
}






