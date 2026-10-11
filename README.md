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
A series of 100 observations is generated randomly, and the rolling mean, variance, skewness, and kurtosis are calculated using a window size of 10. 

```r
## Initial set-up
# 
set.seed(1234)

# Observed data
N <- 100
x <- cumsum(rnorm(N))

# Window size
k <- 10

# Calculating rolling statistics using window size 'k'
mw <- imw(x, k)

#
plot(mw, xlim = c(1, 150))
```
<img width="1533" height="773" alt="Rplot01" src="https://github.com/user-attachments/assets/7df6a8df-26c1-426c-9e99-b799492ce0b0" />

Subsequently, a new batch of 50 observations are observed, and the rolling statistics are updated given this additional information.
```r
## Updating moments
# Additional data
M <- 50
y <- cumsum(c(tail(x, 1), rnorm(M)))[-1]

# Update the moments using new information
umw <- uimw(mw, y)

#
plot(umw, xlim = c(1, 150))
```
<img width="1533" height="779" alt="Rplot02" src="https://github.com/user-attachments/assets/236ba467-fc42-44ae-9bda-505ec0d7947d" />

## License

This project is licensed under the MIT License.

