// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// v3semC
NumericMatrix v3semC(NumericVector par, NumericVector PAR);
RcppExport SEXP PMCMC_v3semC(SEXP parSEXP, SEXP PARSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type PAR(PARSEXP);
    rcpp_result_gen = Rcpp::wrap(v3semC(par, PAR));
    return rcpp_result_gen;
END_RCPP
}
