// 
#include "Rcpp.h"
#include "IncrementalMoments.hpp"

//
using namespace IMW;

//[[Rcpp::export()]] 
Rcpp::List imw_cpp(const std::vector<double> & x, const int & k) {
    ////
    //
    int N = x.size();
    
    //
    IncrementalMoments total;
    IncrementalMoments lagged; 
    
    //
    std::vector<double> mean(N);
    std::vector<double> variance(N);
    std::vector<double> skewness(N);
    std::vector<double> kurtosis(N);
    
    // 
    for (int n = 0; n < k; n++) {
        double x_n = x[n];
        total.add(x_n);
        
        mean[n] = variance[n] = skewness[n] = kurtosis[n] = NA_REAL;
    }
    
    ////
    //
    mean[k - 1] = total.mean();
    variance[k - 1] = total.variance();
    skewness[k - 1] = total.skewness(); 
    kurtosis[k - 1] = total.kurtosis();
    
    // 
    for (int n = k; n < N; n++) {
        //
        double x_n_k = x[n - k];
        double x_n = x[n];
        
        //
        lagged.add(x_n_k);
        total.add(x_n);
        
        //
        IncrementalMoments difference = total - lagged;
        
        //
        mean[n] = difference.mean();
        variance[n] = difference.variance();
        skewness[n] = difference.skewness(); 
        kurtosis[n] = difference.kurtosis();
    }
    
    ////
    //
    Rcpp::List stats = Rcpp::List::create(
        Rcpp::Named("Mean") = mean, 
        Rcpp::Named("Variance") = variance, 
        Rcpp::Named("Skewness") = skewness, 
        Rcpp::Named("Kurtosis") = kurtosis
    );

    //
    return Rcpp::List::create(
        Rcpp::Named("stats") = stats,
        Rcpp::Named("lagged") = lagged.moments(),
        Rcpp::Named("total") = total.moments()
    );
} 

//[[Rcpp::export()]] 
Rcpp::List uimw_cpp(const std::vector<double> & x, const int & k, const std::vector<double> & l, const std::vector<double> & t) {
    ////
    //
    int N = x.size();
    
    //
    IncrementalMoments lagged(l[0], l[1], l[2], l[3], l[4]);
    IncrementalMoments total(t[0], t[1], t[2], t[3], t[4]);
    
    ////
    //
    std::vector<double> mean(N - k);
    std::vector<double> variance(N - k);
    std::vector<double> skewness(N - k);
    std::vector<double> kurtosis(N - k);
    
    //
    for (int n = k; n < N; n++) {
        //
        double x_n_k = x[n - k];
        double x_n = x[n];
        
        //
        lagged.add(x_n_k);
        total.add(x_n);
        
        //
        IncrementalMoments difference = total - lagged;
        
        //
        mean[n - k] = difference.mean();
        variance[n - k] = difference.variance();
        skewness[n - k] = difference.skewness(); 
        kurtosis[n - k] = difference.kurtosis();
    }
    
    ////
    //
    Rcpp::List stats = Rcpp::List::create(
        Rcpp::Named("Mean") = mean, 
        Rcpp::Named("Variance") = variance, 
        Rcpp::Named("Skewness") = skewness, 
        Rcpp::Named("Kurtosis") = kurtosis
    );
    
    //
    return Rcpp::List::create(
        Rcpp::Named("stats") = stats,
        Rcpp::Named("lagged") = lagged.moments(),
        Rcpp::Named("total") = total.moments()
    );
} 
