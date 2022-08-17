#include "func_cpp.h"

using namespace std;
using namespace Eigen;
using namespace Rcpp;
using namespace RcppEigen;

// Sort index by ascending order of vec's values
std::vector<size_t> sortID(const std::vector<double>& vec, const std::vector<size_t>& idx) {
    std::vector<size_t> tmp = idx;
    sort(tmp.begin(), tmp.end(), [&vec](size_t i1, size_t i2) { return vec[i1] < vec[i2]; });
    return tmp;
}

// Return prediction of mapping
// C: coefficient matrix
// X: data matrix
// [[Rcpp::export]]
Eigen::VectorXd getPred(const Eigen::MatrixXd& C, const Eigen::MatrixXd& X) {
    return C.cwiseProduct(X).rowwise().sum();
}

// Return RMSE between v1 and v2
// [[Rcpp::export]]
double getRMSE(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2) {
    return sqrt(pow((v1 - v2).array(), 2.0).mean());
}

// Return MAE between v1 and v2
// [[Rcpp::export]]
double getMAE(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2) {
    return abs((v1 - v2).array()).mean();
}

// Return correlation coefficient between v1 and v2
// [[Rcpp::export]]
double getRho(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2) {
    Eigen::ArrayXd a1 = v1.array() - v1.mean(), a2 = v2.array() - v2.mean();
    return (a1 * a2).sum() / sqrt(pow(a1, 2.0).sum() * pow(a2, 2.0).sum());
}

// Calculate Normal Euclidean Distance (NE-Dist) matrix
// [[Rcpp::export]]
Eigen::MatrixXd getNEDist(const Eigen::MatrixXd& X) {
    size_t m = X.rows();
    Eigen::MatrixXd dist = Eigen::MatrixXd::Zero(m, m);

    for (size_t i = 0; i < m - 1; ++i) {
        for (size_t j = i + 1; j < m; ++j) {
            auto d = sqrt((X.row(i) - X.row(j)).squaredNorm());
            dist(i, j) = dist(j, i) = d;
        }
    }

    return dist;
}

// Prediction with Simplex Projection
// [[Rcpp::export]]
Eigen::VectorXd Smplx(const Eigen::MatrixXd& X, const Eigen::RowVectorXd& y, const Eigen::MatrixXd& dmat,
               const std::vector<size_t>& idxLib, const std::vector<size_t>& idxPrd, size_t nns) {
    // Implementation of NNS (Nearest Neighbour Points) search
    auto find_nn_pts = [&idxLib, &nns](size_t tar, const std::vector<double>& dist) {
        std::vector<size_t> idx = sortID(dist, idxLib);
        std::vector<size_t> nn_pts;

        for (size_t i : idx) {
            if (i != tar) nn_pts.push_back(i);
            if (nn_pts.size() == nns) break;
        }

        double maxD = nn_pts.back();

        // Check whether same distance exists after nns'th state
        for (size_t i = nns; i < idx.size(); ++i) {
            if (dist[idx[i]] != maxD) break;
            nn_pts.push_back(idx[i]);
        }

        return nn_pts;
    };

    size_t xm = X.rows();
    std::vector<double> pred;

    for (size_t tar : idxPrd) {
        std::vector<double> dist(xm);
        Map<Eigen::VectorXd>(&dist[0], xm) = dmat.row(tar);

        std::vector<size_t> nn_pts = find_nn_pts(tar, dist);
        double min_dist = dist[nn_pts[0]], tmp = 0.0;
        size_t len = nn_pts.size();
        Eigen::VectorXd weights = Eigen::VectorXd::Zero(len);

        if (min_dist == 0.0) {
            for (size_t i = 0; i < len; ++i) {
                if (dist[nn_pts[i]] == min_dist) {
                    weights(i) = 1.0;
                    tmp += y(nn_pts[i]);
                } else {
                    break;
                }
            }
        } else {
            for (size_t i = 0; i < len; ++i) {
                weights(i) = exp(-dist[nn_pts[i]] / min_dist);
                tmp += weights(i) * y(nn_pts[i]);
            }
        }

        pred.push_back(tmp / weights.sum());
    }

    return Map<Eigen::VectorXd>(&pred[0], pred.size());
}

