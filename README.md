# IMW
The `IMW`-package implements a method for calculating the rolling statistics mean, variance, skewness, and kurtosis. Furthermore, it provides a function for updating these statistics as new information is revieved, allowing for online/batch calculation of these statistics. 

## Installation

The `IMW`-package depends on `R` (>= 4.0.1), `Rcpp` (>= 1.0.4.6), and `RcppArmadillo`. As the package is not available on CRAN, devtools is needed to install the package from github. 

From R, run the following commands:  

```r
install.packages("Rcpp")
install.packages("RcppArmadillo")

istall.packages("devtools")
devtools::install_github("svilsen/IMW")
```

## Usage
In the following, a series of size 100 is generated randomly, then the rolling statistics of mean, variance, skewness, and kurtosis are calculated using a window size of size 10. Afterwhich new batch of 20 observations are observed, and the rolling statistics are updated given the additional information.

```r
## Data is observed
N <- 100
x <- cumsum(rnorm(N))

## Window size
k <- 10

## Calculating rolling statistics using window size 'k'
x_imw <- imw(x, k)

## New data is observed
N_new <- 20
x_new <- cumsum(c(tail(x, 1), rnorm(N_new)))[-1]

## Update the moments using new information
x_imw_new <- update_imw(x_imw, x_new)
```

## License

This project is licensed under the MIT License.

