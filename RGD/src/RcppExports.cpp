// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// GD
arma::vec GD(arma::vec beta_0, const double tau, const arma::mat X, const arma::mat Y, double eta_0, const double alpha);
RcppExport SEXP _RGD_GD(SEXP beta_0SEXP, SEXP tauSEXP, SEXP XSEXP, SEXP YSEXP, SEXP eta_0SEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type beta_0(beta_0SEXP);
    Rcpp::traits::input_parameter< const double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type eta_0(eta_0SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(GD(beta_0, tau, X, Y, eta_0, alpha));
    return rcpp_result_gen;
END_RCPP
}
// runif_in_pball
arma::mat runif_in_pball(const int n, const int d, const int p, const double r);
RcppExport SEXP _RGD_runif_in_pball(SEXP nSEXP, SEXP dSEXP, SEXP pSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(runif_in_pball(n, d, p, r));
    return rcpp_result_gen;
END_RCPP
}
// RGD
arma::vec RGD(const arma::mat X, const arma::mat Y, const double tau, const double iter, double eta_0, const double alpha);
RcppExport SEXP _RGD_RGD(SEXP XSEXP, SEXP YSEXP, SEXP tauSEXP, SEXP iterSEXP, SEXP eta_0SEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const double >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< double >::type eta_0(eta_0SEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(RGD(X, Y, tau, iter, eta_0, alpha));
    return rcpp_result_gen;
END_RCPP
}
// cplexcoef
Rcpp::List cplexcoef(const arma::mat X, const arma::mat Y, const double tau);
RcppExport SEXP _RGD_cplexcoef(SEXP XSEXP, SEXP YSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(cplexcoef(X, Y, tau));
    return rcpp_result_gen;
END_RCPP
}
// picksamples
Rcpp::List picksamples(const arma::mat X, const arma::mat Y, const double n_ratio);
RcppExport SEXP _RGD_picksamples(SEXP XSEXP, SEXP YSEXP, SEXP n_ratioSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type n_ratio(n_ratioSEXP);
    rcpp_result_gen = Rcpp::wrap(picksamples(X, Y, n_ratio));
    return rcpp_result_gen;
END_RCPP
}
