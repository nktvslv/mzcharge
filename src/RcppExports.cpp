// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// find_centroids
List find_centroids(const NumericVector& spectrum, int min_width, double min_intensity);
RcppExport SEXP _mzcharge_find_centroids(SEXP spectrumSEXP, SEXP min_widthSEXP, SEXP min_intensitySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type spectrum(spectrumSEXP);
    Rcpp::traits::input_parameter< int >::type min_width(min_widthSEXP);
    Rcpp::traits::input_parameter< double >::type min_intensity(min_intensitySEXP);
    rcpp_result_gen = Rcpp::wrap(find_centroids(spectrum, min_width, min_intensity));
    return rcpp_result_gen;
END_RCPP
}
// charge_singlespec
DataFrame charge_singlespec(const NumericVector& spectrum, int min_width, double min_intensity, double mz_tol, double intensity_tol, int max_charge, int num_iso);
RcppExport SEXP _mzcharge_charge_singlespec(SEXP spectrumSEXP, SEXP min_widthSEXP, SEXP min_intensitySEXP, SEXP mz_tolSEXP, SEXP intensity_tolSEXP, SEXP max_chargeSEXP, SEXP num_isoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type spectrum(spectrumSEXP);
    Rcpp::traits::input_parameter< int >::type min_width(min_widthSEXP);
    Rcpp::traits::input_parameter< double >::type min_intensity(min_intensitySEXP);
    Rcpp::traits::input_parameter< double >::type mz_tol(mz_tolSEXP);
    Rcpp::traits::input_parameter< double >::type intensity_tol(intensity_tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_charge(max_chargeSEXP);
    Rcpp::traits::input_parameter< int >::type num_iso(num_isoSEXP);
    rcpp_result_gen = Rcpp::wrap(charge_singlespec(spectrum, min_width, min_intensity, mz_tol, intensity_tol, max_charge, num_iso));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mzcharge_find_centroids", (DL_FUNC) &_mzcharge_find_centroids, 3},
    {"_mzcharge_charge_singlespec", (DL_FUNC) &_mzcharge_charge_singlespec, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_mzcharge(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
