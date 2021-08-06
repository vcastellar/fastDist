// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// euclidean
NumericMatrix euclidean(NumericMatrix Ar, NumericMatrix Br);
RcppExport SEXP _fastDist_euclidean(SEXP ArSEXP, SEXP BrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Ar(ArSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Br(BrSEXP);
    rcpp_result_gen = Rcpp::wrap(euclidean(Ar, Br));
    return rcpp_result_gen;
END_RCPP
}
// manhattan
NumericMatrix manhattan(NumericMatrix Ar, NumericMatrix Br);
RcppExport SEXP _fastDist_manhattan(SEXP ArSEXP, SEXP BrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Ar(ArSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Br(BrSEXP);
    rcpp_result_gen = Rcpp::wrap(manhattan(Ar, Br));
    return rcpp_result_gen;
END_RCPP
}
// minkowsky
NumericMatrix minkowsky(NumericMatrix Ar, NumericMatrix Br, double p);
RcppExport SEXP _fastDist_minkowsky(SEXP ArSEXP, SEXP BrSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Ar(ArSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Br(BrSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(minkowsky(Ar, Br, p));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fastDist_euclidean", (DL_FUNC) &_fastDist_euclidean, 2},
    {"_fastDist_manhattan", (DL_FUNC) &_fastDist_manhattan, 2},
    {"_fastDist_minkowsky", (DL_FUNC) &_fastDist_minkowsky, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_fastDist(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
