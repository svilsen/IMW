#' @title Incremental Moment Windows
#' 
#' @description Rolling estimates of moments given a fixed window size. 
#' 
#' @param x Vector, or column matrix: The data used to estimate the moments.
#' @param k Numeric: The window size.
#' 
#' @return An \link{imw-object}.
#' 
#' @examples N <- 100
#' x <- cumsum(rnorm(N))
#' 
#' k <- 10
#' imw(x, k)
#'  
#' @export 
imw <- function(x, k) {
    #
    if (!is.matrix(x)) {
        x <- matrix(x, ncol = 1)
    }
    
    #
    if (!is.integer(k)) {
        k <- as.integer(k)
    }
    
    #
    res <- imw_cpp(x, k)
    
    #
    res[["stats"]] <- do.call("cbind", res[["stats"]])
    res[["x"]] <- x
    res[["k"]] <- k
    
    #
    class(res) <- "imw"
    return(res)
}

#' @title Rolling Mean
#' 
#' @description Returns the rolling mean of an \link{imw-object}.
#' 
#' @param x An \link{imw-object}.
#' @param ... Additional arguments; not used in this instance.
#' 
#' @return The estimated rolling mean.
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
    res <- x[["stats"]][, 1]
    return(res)
}

#' @title Rolling Variance
#' 
#' @description Returns the rolling variance of an \link{imw-object}.
#' 
#' @param x An \link{imw-object}.
#' @param ... Additional arguments, see details.
#' 
#' @details The only additional argument is \code{type} which takes the values \code{1} and \code{2} corresponding to the population and sample variance, respectively.  
#' 
#' @return The estimated rolling variance.
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

#' @title Rolling Variance
#' 
#' @description Returns the rolling variance of an \link{imw-object}.
#' 
#' @param x An \link{imw-object}.
#' @param ... Additional arguments, see details.
#' 
#' @details The only additional argument is \code{type} which takes the values \code{1} and \code{2} corresponding to the population and sample variance, respectively.  
#' 
#' @return The estimated rolling variance.
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
    
    m <- x[["stats"]][, 2]
    if (type == 1) {
        res <- m
    }
    else if (type == 2) {
        k <- x[["k"]]
        res <- (k - 1) * m / k
    }
    else {
        stop("'type' has to be 1 or 2.")
    }
    
    return(res)
}

#' @title Rolling Skewness
#' 
#' @description Returns the rolling skewness of an \link{imw-object}.
#' 
#' @param x An \link{imw-object}.
#' @param ... Additional arguments, see details.
#' 
#' @details The only additional argument is \code{type} which takes the values \code{1, 2, 3}, corresponding to the three methods discussed by Joanes and Gill (1998).  
#' 
#' @return The estimated rolling skewness.
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

#' @title Rolling Skewness
#' 
#' @description Returns the rolling skewness of an \link{imw-object}.
#' 
#' @param x An \link{imw-object}.
#' @param ... Additional arguments, see details.
#' 
#' @details The only additional argument is \code{type} which takes the values \code{1, 2, 3}, corresponding to the three methods discussed by Joanes and Gill (1998).  
#' 
#' @return The estimated rolling skewness.
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
    
    m <- x[["stats"]][, 3]
    k <- x[["k"]]
        
    if (type == 1) {
        res <- m
    }
    else if (type == 2) {
        if (k < 3) {
            stop("The window size 'k' needs to be at least 3.")
        }
        
        res <- m * sqrt(k * (k - 1)) / (k - 2)
    }
    else if (type == 3) {
        res <- m * (1.0 - 1.0 / k)^(3 / 2)
    }
    else {
        stop("'type' has to be 1, 2, or 3.")
    }
    
    return(res)
}

#' @title Rolling Kurtosis
#' 
#' @description Returns the rolling kurtosis of an \link{imw-object}.
#' 
#' @param x An \link{imw-object}.
#' @param ... Additional arguments, see details.
#' 
#' @details The only additional argument is \code{type} which takes the values \code{1, 2, 3}, corresponding to the three methods discussed by Joanes and Gill (1998).  
#' 
#' @return The estimated rolling kurtosis.
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

#' @title Rolling Kurtosis
#' 
#' @description Returns the rolling kurtosis of an \link{imw-object}.
#' 
#' @param x An \link{imw-object}.
#' @param ... Additional arguments, see details.
#' 
#' @details The only additional argument is \code{type} which takes the values \code{1, 2, 3}, corresponding to the three methods discussed by Joanes and Gill (1998).  
#' 
#' @return The estimated rolling kurtosis.
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
    
    m <- x[["stats"]][, 4]
    k <- x[["k"]]
    
    if (type == 1) {
        res <- m
    }
    else if (type == 2) {
        if (k < 4) {
            stop("The window size 'k' needs to be at least 4.")
        }
        
        res <- ((k + 1) * m + 6.0) * (k - 1) / ((k - 2) * (k - 3))
        
    }
    else if (type == 3) {
        res <- (m + 3) * (1.0 + 1.0 / k)^2 - 3.0
    }
    else {
        stop("The supplied 'type' was not valid.")
    }
    
    return(res)
}

#' @title Matrix Conversion
#' 
#' @description Converts an object of \link{imw-object} to a matrix.
#' 
#' @param x An \link{imw-object}.
#' @param ... Additional arguments; not used in this instance.
#' 
#' @rdname as.matrix.imw
#' @method as.matrix imw
#' 
#' @return A numeric matrix of rolling statistics.
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
    res <- x[["stats"]]
    return(res)
}

#' @title Update Incremental Moment Windows
#' 
#' @description Update an \link{imw-object} with new data.
#' 
#' @param object An \link{imw-object}.
#' @param xnew Vector, or column matrix: The new data used to update the estimated moments.
#' 
#' @return An \link{imw-object}.
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
uimw <- function(object, xnew) {
    if (!is.matrix(xnew)) {
        xnew <- matrix(xnew, ncol = 1)
    }
    
    x <- object[["x"]]
    k <- object[["k"]]
    N <- length(x)
    
    xc <- rbind(tail(x, k), xnew)
    res <- uimw_cpp(x = xc, k = k, l = object[["lagged"]], t = object[["total"]])
    
    res[["stats"]] <- do.call("cbind", res[["stats"]])
    res[["stats"]] <- rbind(object[["stats"]], res[["stats"]])
    res[["x"]] <- rbind(x, xnew)
    res[["k"]] <- k
    
    class(res) <- "imw"
    return(res)
}

#' @title Plot of rolling statistics
#' 
#' @description Shows the four rolling statistics of an \link{imw-object}.
#' 
#' @param x An \link{imw-object}.
#' @param ... Additional arguments passed to the internal plot functions.
#' 
#' @return An invisible \code{NULL}.
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
    ##
    #
    oldpar <- par(no.readonly = TRUE) 
    on.exit(par(oldpar))
    
    ##
    #
    par(mfrow = c(2, 2), mar = c(5, 5, 0.5, 0.5))
    
    #
    plot(x[["x"]], pch = 16, xlab = "t", ylab = expression("x"[" t"]), ...)
    points(mean(x), pch = 16, col = "red", type = "l", lwd = 3)
    
    #
    plot(variance(x), pch = 16, xlab = "t", ylab = expression("Variance(x"[" t "]*")"), ...)
    abline(h = mean(variance(x), na.rm = TRUE), col = "red", lwd = 3)
    
    #
    plot(skewness(x), pch = 16, xlab = "t", ylab = expression("Skewness(x"[" t "]*")"), ...)
    abline(h = 0, col = "red", lwd = 3)
    
    #
    plot(kurtosis(x), pch = 16, xlab = "t", ylab = expression("Kurtosis(x"[" t "]*")"), ...)
    abline(h = 0, col = "red", lwd = 3)
    
    ##
    #
    return(invisible(NULL))
}

