#ifndef FUNC_CPP_H
#define FUNC_CPP_H

#include <vector>
#include <Eigen/Dense>
#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]

// Forward declarations of functions
std::vector<size_t> sortID(const std::vector<double>& vec, const std::vector<size_t>& idx);

Eigen::VectorXd getPred(const Eigen::MatrixXd& C, const Eigen::MatrixXd& X);

double getRMSE(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2);

double getMAE(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2);

double getRho(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2);

Eigen::MatrixXd getNEDist(const Eigen::MatrixXd& X);

Eigen::VectorXd Smplx(const Eigen::MatrixXd& X, const Eigen::RowVectorXd& y, const Eigen::MatrixXd& dmat, const std::vector<size_t>& idxLib, const std::vector<size_t>& idxPrd, size_t nns);

#endif

