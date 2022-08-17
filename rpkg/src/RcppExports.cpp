// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <vector>
#include <Eigen/Dense>
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace std;
using namespace Eigen;
using namespace Rcpp;
using namespace RcppEigen;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// getPred
VectorXd getPred(const MatrixXd& C, const MatrixXd& X);
RcppExport SEXP _rpkg_getPred(SEXP CSEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MatrixXd& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(getPred(C, X));
    return rcpp_result_gen;
END_RCPP
}
// getRMSE
double getRMSE(const VectorXd& v1, const VectorXd& v2);
RcppExport SEXP _rpkg_getRMSE(SEXP v1SEXP, SEXP v2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd& >::type v1(v1SEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type v2(v2SEXP);
    rcpp_result_gen = Rcpp::wrap(getRMSE(v1, v2));
    return rcpp_result_gen;
END_RCPP
}
// getMAE
double getMAE(const VectorXd& v1, const VectorXd& v2);
RcppExport SEXP _rpkg_getMAE(SEXP v1SEXP, SEXP v2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd& >::type v1(v1SEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type v2(v2SEXP);
    rcpp_result_gen = Rcpp::wrap(getMAE(v1, v2));
    return rcpp_result_gen;
END_RCPP
}
// getRho
double getRho(const VectorXd& v1, const VectorXd& v2);
RcppExport SEXP _rpkg_getRho(SEXP v1SEXP, SEXP v2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd& >::type v1(v1SEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type v2(v2SEXP);
    rcpp_result_gen = Rcpp::wrap(getRho(v1, v2));
    return rcpp_result_gen;
END_RCPP
}
// getNEDist
MatrixXd getNEDist(const MatrixXd& X);
RcppExport SEXP _rpkg_getNEDist(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MatrixXd& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(getNEDist(X));
    return rcpp_result_gen;
END_RCPP
}
// Smplx
VectorXd Smplx(const MatrixXd& X, const RowVectorXd& y, const MatrixXd& dmat, const vector<size_t>& idxLib, const vector<size_t>& idxPrd, size_t nns);
RcppExport SEXP _rpkg_Smplx(SEXP XSEXP, SEXP ySEXP, SEXP dmatSEXP, SEXP idxLibSEXP, SEXP idxPrdSEXP, SEXP nnsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const RowVectorXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type dmat(dmatSEXP);
    Rcpp::traits::input_parameter< const vector<size_t>& >::type idxLib(idxLibSEXP);
    Rcpp::traits::input_parameter< const vector<size_t>& >::type idxPrd(idxPrdSEXP);
    Rcpp::traits::input_parameter< size_t >::type nns(nnsSEXP);
    rcpp_result_gen = Rcpp::wrap(Smplx(X, y, dmat, idxLib, idxPrd, nns));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rpkg_getPred", (DL_FUNC) &_rpkg_getPred, 2},
    {"_rpkg_getRMSE", (DL_FUNC) &_rpkg_getRMSE, 2},
    {"_rpkg_getMAE", (DL_FUNC) &_rpkg_getMAE, 2},
    {"_rpkg_getRho", (DL_FUNC) &_rpkg_getRho, 2},
    {"_rpkg_getNEDist", (DL_FUNC) &_rpkg_getNEDist, 1},
    {"_rpkg_Smplx", (DL_FUNC) &_rpkg_Smplx, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_rpkg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}