####
#### R function file for "Flexibility of species interaction and ecological stability"
#### func_common.R: contains common functions throughout the analyses
####
####
#### ver 0.0.1: Initially written on 220210319 by K.Kawatsu
#### ver 0.0.2: Updated on 20210326
#### ver 0.0.3: Updated on 20210329
#### ver 0.0.4: Updated on 20210402
#### ver 0.1.0: Updated on 20210617
#### ver 0.1.1: Updated on 20210621
#### ver 0.1.2: Updated on 20210702 (added rho, MAE & RMSE function, causality test function)
#### ver 0.2.0: Updated on 20210706
#### ver 0.3.0: Updated on 20211117
####

#' R wrapper function of getRMSE in func_cpp.cpp
#'
#' \code{get_rmse} returns the RMSE (Root Mean Squared Error) between x & y
#'
#' @param x Numeric vector
#' @param y Numeric vector
#' @export
get_rmse <- function(x, y) {
    idx <- complete.cases(cbind(x, y))
    res <- ifelse(any(idx), getRMSE(x[idx], y[idx]), NA)
    return(res)
}

#' R wrapper function of getMAE in func_cpp.cpp
#'
#' \code{get_mae} returns the MAE (Mean Absolute Error) between x & y
#'
#' @inheritParams get_rmse
#' @export
get_mae <- function(x, y) {
    idx <- complete.cases(cbind(x, y))
    res <- ifelse(any(idx), getMAE(x[idx], y[idx]), NA)
    return(res)
}

#' R wrapper function of getRho in func_cpp.cpp
#'
#' \code{get_rho} returns the correlation coefficient between x & y
#' 
#' @inheritParams get_rmse
#' @export
get_rho <- function(x, y) {
    idx <- complete.cases(cbind(x, y))
    res <- ifelse(any(idx), getRho(x[idx], y[idx]), NA)
    return(res)
}

#' Calculation of CV (Coefficient of Variance)
#'
#' \code{get_cv} returns the CV between x & y
#'
#' @param v Numeric vector
#' @export
get_cv <- function(v) sd(v, na.rm = TRUE) / abs(mean(v, na.rm = TRUE))

#' R wrapper function of getPred in func_cpp.cpp
#'
#' @param C Model coefficient matrix
#' @param X Data matrix
#' @export
get_pred <- function(C, X) return(getPred(C, X))

#' Shift a vector forward or backward depending on the lag value
#'
#' \code{shift} shifts an input vector v backward or forward
#' if lag < 0: backward and otherwise forward
#'
#' @param v: input vector
#' @param lag: lag value for shifting v
#' @export
shift <- function(v, lag) do.call(if_else(lag < 0, "lag", "lead"), args = list(v, abs(lag)));

#' R wrapper function of getNEDist in func_cpp.cpp
#'
#' \code{get_NED} returns the Normal Euclidean Distance (NED) of the matrix X
#'
#' @param X Data matrix
#' @export
get_ned <- function(X) {
    if (!is.matrix(X)) X <- matrix(X, ncol = 1)
    return(getNEDist(X))
}

#' Normalization function
#'
#' @param x Numeric vector, to be normalized
#' @param method Character, which switches the normalization method
#'      if method == "0-1", the x is normalized to range between 0-1
#'      otherwise, x is normalized with mean 0 and unit variance
#' @export
normalize <- function(x, method = "0-1") {
    if (method == "0-1") {
        x_ <- (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    } else {
        x_ <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
    }

    return(x_)
}
