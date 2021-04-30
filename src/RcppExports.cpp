// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// imw_cpp
Rcpp::List imw_cpp(const arma::colvec& x, const int& k);
RcppExport SEXP _IMW_imw_cpp(SEXP xSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(imw_cpp(x, k));
    return rcpp_result_gen;
END_RCPP
}
// imw_update_cpp
Rcpp::List imw_update_cpp(const arma::colvec& x, const int& k, const arma::colvec& t, const arma::colvec& l);
RcppExport SEXP _IMW_imw_update_cpp(SEXP xSEXP, SEXP kSEXP, SEXP tSEXP, SEXP lSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type l(lSEXP);
    rcpp_result_gen = Rcpp::wrap(imw_update_cpp(x, k, t, l));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_IMW_imw_cpp", (DL_FUNC) &_IMW_imw_cpp, 2},
    {"_IMW_imw_update_cpp", (DL_FUNC) &_IMW_imw_update_cpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_IMW(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}