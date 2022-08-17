#' rpkg: an R package for "Flexibility of species interactions and ecological stability"
#'  This package includes functions for causality test and Jacobian-estimation
#'
#' @docType package
#' @name rpkg
#'
#' @useDynLib rpkg, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom dplyr lag
#' @importFrom dplyr lead
#' @importFrom dplyr if_else
#' @importFrom dplyr arrange
#' @importFrom dplyr slice
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom tibble tibble
#' @importFrom tibble is_tibble
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom stringr str_c
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach %:%
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#' @importFrom glmnet cv.glmnet
NULL

