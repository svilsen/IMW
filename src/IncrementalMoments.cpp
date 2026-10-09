// 
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include "IncrementalMoments.hpp"
using namespace IMW;

//[[Rcpp::export()]] 
Rcpp::List imw_cpp(const arma::colvec & x, const int & k) {
    //
    int N = x.size();
    
    IncrementalMoments total;
    IncrementalMoments lagged; 
    
    arma::mat R(N, 4);
    
    // 
    for (int n = 0; n < k; n++) {
        double x_n = x[n];
        total.add(x_n);
        
        R(n, 0) = R(n, 1) = R(n, 2) = R(n, 3) = NA_REAL;
    }
    
    //
    R(k - 1, 0) = total.mean();
    R(k - 1, 1) = total.variance();
    R(k - 1, 2) = total.skewness(); 
    R(k - 1, 3) = total.kurtosis();
    
    // 
    for (int n = k; n < N; n++) {
        double x_n = x[n];
        total.add(x_n);
        
        double x_n_k = x[n - k];
        lagged.add(x_n_k);
        
        IncrementalMoments difference = total - lagged;
        
        R(n, 0) = difference.mean();
        R(n, 1) = difference.variance();
        R(n, 2) = difference.skewness(); 
        R(n, 3) = difference.kurtosis();
    }
    
    //
    return Rcpp::List::create(
        Rcpp::Named("stats") = R,
        Rcpp::Named("totalMoments") = total.moments(),
        Rcpp::Named("windowMoments") = lagged.moments()
    );
} 

//[[Rcpp::export()]] 
Rcpp::List uimw_cpp(const arma::colvec & x, const int & k, const arma::colvec & t, const arma::colvec & l) {
    //
    int N = x.size();
    
    IncrementalMoments total(t[0], t[1], t[2], t[3], t[4]);
    IncrementalMoments lagged(l[0], l[1], l[2], l[3], l[4]);
    
    //
    arma::mat R(N - k, 4);
    for (int n = k; n < N; n++) {
        double x_n = x[n];
        //
        total.add(x_n);
        IncrementalMoments delta_ = total - lagged;
        
        //
        double x_n_k = x[n - k];
        lagged.add(x_n_k);
        IncrementalMoments delta = total - lagged;
        
        //
        R(n - k, 0) = delta.mean();
        R(n - k, 1) = delta.variance();
        R(n - k, 2) = delta.skewness(); 
        R(n - k, 3) = delta.kurtosis();
    }
    
    //
    return Rcpp::List::create(
        Rcpp::Named("stats") = R,
        Rcpp::Named("totalMoments") = total.moments(),
        Rcpp::Named("windowMoments") = lagged.moments()
    );
} 
