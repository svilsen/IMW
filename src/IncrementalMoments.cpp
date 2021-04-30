#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

class IncrementalMoments
{
public:
    // 
    int N;
    double M1, M2, M3, M4;
    
    //
    IncrementalMoments() {
        N = 0;
        M1 = M2 = M3 = M4 = 0.0;
    }
    
    IncrementalMoments(const double & _M1, const double & _M2, const double & _M3, const double & _M4, const int & _N) {
        N = _N;
        M1 = _M1; 
        M2 = _M2;
        M3 = _M3;
        M4 = _M4;
    }
    
    //
    void add(const double & x) {
        int M = N;
        N++;
        
        double D = x - M1;
        double D_N = D / N;
        double D_N_2 = D_N * D_N;
        double tau = D * D_N * M;
        
        M4 += tau * D_N_2 * (N * N - 3.0 * N + 3) + 6 * D_N_2 * M2 - 4 * D_N * M3;
        M3 += tau * D_N * (N - 2) - 3.0 * D_N * M2;
        M2 += tau;
        M1 += D_N;
    }
    
    // 
    double mean() {
        return M1;
    }
    
    double variance() {
        return M2 / (N - 1.0);
    }
    
    double skewness() {
        return sqrt(double(N)) * M3 / pow(M2, 1.5);
    }
    
    double kurtosis() {
        return double(N) * M4 / (M2 * M2) - 3.0;
    }
    
    arma::colvec moments() {
        arma::colvec x(5);
        x[0] = M1;
        x[1] = M2;
        x[2] = M3;
        x[3] = M4;
        x[4] = N;
        
        return x;
    }
};


//
IncrementalMoments operator+(const IncrementalMoments A, const IncrementalMoments B) {
    //
    IncrementalMoments C;
    
    //
    C.N = A.N + B.N;
    C.M1 = (A.N * A.M1 + B.N * B.M1) / C.N;
    
    //
    double D = B.M1 - A.M1;
    double D_2 = D * D;
    C.M2 = A.M2 + B.M2 + D_2 * A.N * B.N / C.N;
    
    double D_3 = D * D_2;
    C.M3 = A.M3 + B.M3 + 
        D_3 * A.N * B.N * (A.N - B.N) / (C.N * C.N) +
        3.0 * D * (A.N * B.M2 - B.N * A.M2) / C.N;
    
    double D_4 = D_2 * D_2;
    C.M4 = A.M4 + B.M4 + D_4 * A.N * B.N * (A.N * A.N - A.N * B.N + B.N * B.N) / (C.N * C.N * C.N) + 
        6.0 * D_2 * (A.N * A.N * B.M2 + B.N * B.N * A.M2) / (C.N * C.N) + 
        4.0 * D * (A.N * B.M3 - B.N * A.M3) / C.N;
    
    return C;
}

IncrementalMoments operator-(const IncrementalMoments A, const IncrementalMoments B) {
    IncrementalMoments C;
    
    //
    C.N = A.N - B.N;
    C.M1 = (A.N * A.M1 - B.N * B.M1) / C.N;
    
    //
    double D = C.M1 - B.M1;
    double D_2 = D * D;
    C.M2 = 
        A.M2 - 
        B.M2 - 
        D_2 * B.N * C.N / A.N;
    
    double D_3 = D * D_2;
    C.M3 = 
        A.M3 - 
        B.M3 - 
        D_3 * B.N * C.N * (B.N - C.N) / (A.N * A.N) -
        3.0 * D * (B.N * C.M2 - C.N * B.M2) / A.N;
    
    double D_4 = D_2 * D_2;
    C.M4 = 
        A.M4 - 
        B.M4 - 
        D_4 * B.N * C.N * (B.N * B.N - B.N * C.N + C.N * C.N) / (A.N * A.N * A.N) - 
        6.0 * D_2 * (B.N * B.N * C.M2 + C.N * C.N * B.M2) / (A.N * A.N) - 
        4.0 * D * (B.N * C.M3 - C.N * B.M3) / A.N;
    
    return C;
}

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
Rcpp::List imw_update_cpp(const arma::colvec & x, const int & k, const arma::colvec & t, const arma::colvec & l) {
    //
    int N = x.size();
    
    IncrementalMoments total(t[0], t[1], t[2], t[3], t[4]);
    IncrementalMoments lagged(l[0], l[1], l[2], l[3], l[4]);
    
    //
    arma::mat R(N - k, 4);
    for (int n = k; n < N; n++) {
        double x_n = x[n];
        total.add(x_n);
        
        double x_n_k = x[n - k];
        lagged.add(x_n_k);
        
        IncrementalMoments difference = total - lagged;
        
        R(n - k, 0) = difference.mean();
        R(n - k, 1) = difference.variance();
        R(n - k, 2) = difference.skewness(); 
        R(n - k, 3) = difference.kurtosis();
    }
    
    //
    return Rcpp::List::create(
        Rcpp::Named("stats") = R,
        Rcpp::Named("totalMoments") = total.moments(),
        Rcpp::Named("windowMoments") = lagged.moments()
    );
} 
