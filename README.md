# IMW
The `IMW`-package implements a method for calculating the rolling statistics mean, variance, skewness, and kurtosis. Furthermore, it provides a function for updating these statistics as new information is received, allowing for online/batch calculation of these statistics. 

## Installation

The `IMW`-package depends on `R` (>= 4.0.1), `Rcpp` (>= 1.0.4.6), and `RcppArmadillo`. As the package is not available on CRAN, `remotes` is needed to install the package from github. 

From R, run the following commands:  

```r
install.packages("Rcpp")
install.packages("RcppArmadillo")

if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

remotes::install_github("svilsen/IMW")
```

## Usage
In the following, a series of 100 observations is generated randomly, and the rolling mean, variance, skewness, and kurtosis are calculated using a window size of 10. Subsequently, a new batch of 20 observations are observed, and the rolling statistics are updated given this additional information.

```r
## Initial set-up
# Observed data
N <- 100
x <- cumsum(rnorm(N))

# Window size
k <- 10

# Calculating rolling statistics using window size 'k'
mw <- imw(x, k)

## Updating moments
# Additional data
M <- 20
y <- cumsum(c(tail(x, 1), rnorm(M)))[-1]

# Update the moments using new information
umw <- uimw(mw, y)

#
plot(umw)
```

## License

This project is licensed under the MIT License.

