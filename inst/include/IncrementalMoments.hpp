#ifndef IncrementalMoments_H
#define IncrementalMoments_H

//
#include "Rcpp.h"

//
namespace IMW {
    class IncrementalMoments {
    public:
        // Data
        int N;
        double M1, M2, M3, M4;
        
        // Constructors
        IncrementalMoments() : N(0), M1(0.0), M2(0.0), M3(0.0), M4(0.0) { }
        IncrementalMoments(const double & _M1, const double & _M2, const double & _M3, const double & _M4, const int & _N) : 
            N(_N), M1(_M1), M2(_M2), M3(_M3), M4(_M4) { } 
        
        // Operator
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
        
        // Moments
        std::vector<double> moments() {
            std::vector<double> x(5);
            x[0] = M1;
            x[1] = M2;
            x[2] = M3;
            x[3] = M4;
            x[4] = N;
            
            return x;
        }
        
        // Statistics
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
    };
    
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
}

#endif //IncrementalMoments_H