####
#### R functions file for "Flexibility of species interaction and ecological stability"
#### func_ntsa.R: contains functions related to nonlinear time-series analysis
####
#### ver 0.0.1: Initially written on 20210706 by K.Kawatsu
#### ver 0.0.2: Updated on 20210714
#### ver 0.1.0: Updated on 20211118
####

#' Find knot positions in data
#'
#' \code{find_knot} returns a 2-column matrix, where start and end point of an i'th knot in X
#'
#' @param df data.frame or tibble that contains the time-series data
#' @param key_col Numeric or character, which contains the lables of the different data group
#' @param time_col Numeric or character, which specifies the time-point column
#' @param diff Numeric, which sets the value of increment step in 'time_col'
#' @export
find_knot <- function(df, key_col = NULL, time_col = 1, diff = 1) {
    if (is.null(key_col)) {
        tmp <- df %>% pull(time_col)
        bwd <- tmp - lag(tmp)
        fwd <- lead(tmp) - tmp
        knot <- cbind(s = which(is.na(bwd) | bwd != diff), e = which(is.na(fwd) | fwd != diff))
    } else {
        knot <- foreach(k = unique(pull(df, key_col)), .combine = rbind) %do% {
            tmp <- which(pull(df, key_col) == k)
            matrix(c(s = min(tmp), e = max(tmp)), nrow = 1)
        }
    }

    return(knot)
}

#' Generate valid time-index from knot data
#'
#' \code{gen_valid_idx} returns a vector, which contains positions of valid rows in knot
#' 
#' @param knot vector or 2-column matrix, which specifies the knot position in data matrix
#' @param idx_all a vector contains valid position in the whole data matrix
gen_valid_idx <- function(knot, idx_all) {
    tmp <- foreach(i = 1:nrow(knot), .combine = c) %do% knot[i, 1]:knot[i, 2]
    return(tmp[tmp %in% idx_all])
}

#' Embedding matrix generation with multivariate and multi-timelags
#'
#' \code{gen_emat} returns an embedding matrix with multivariate and multi-timelags
#'
#' @param df vector, matrix, data.frame or tibble,
#'      which contains time series to be embedded
#' @param cols Numeric or character vector,
#'      which selects the df's column for the embedding matrix reconstruction
#' @param lags Integer vector, which sets
#'      the time delay value for each coordinate in embedding matrix
#' @param knot vector or 2-column matrix, which specifies the knot
#'      positions in the data
#' @export
gen_emat <- function(df, cols, lags, knot = matrix(c(1, nrow(df)), nrow = 1)) {
    if (!is.matrix(knot)) knot <- matrix(knot, nrow = 1)
    if (ncol(knot) != 2) stop("Inappropriate style knot!")
    if (!is_tibble(df)) df <- df %>% as.data.frame() %>% setNames(str_c("x", 1:ncol(.))) %>% as_tibble()

    emat <- matrix(NA, nrow = max(knot) - min(knot) + 1, ncol = length(cols))

    for (i in 1:nrow(knot)) {
        idx <- knot[i, 1]:knot[i, 2]

        for (j in 1:length(cols)) {
            tmp <- df %>% slice(idx) %>% pull(cols[j])
            emat[idx, j] <- shift(tmp, lags[j])
        }
    }

    return(emat)
}

#' Cross-mapping with the information of nearest-neighbour points
#'
#' \code{xmap} returns the prediction skill of ref cross-mapping tar
#' cross-mapping skill is calculated with lmdSmplx in ands_ntsa.cpp
#'
#' @param df vector, matrix, data.frame or tibble,
#'      which contains time series to be analyzed
#' @param ref Numeric or character, which sets the library variable for xmap
#' @param tar Numeric or character, which sets the target variable for xmap
#' @param lib Vector or 2-column matrix, which sets the knot position in library data
#' @param E Integer, which is the dimension of ref's embedding matrix
#' @param Tp Integer, which is the prediction-time horizon in ref xmap tar
#' @param libSize Integer, which sets the size of library data used
#' @param randomLibs Logical
#'      if TRUE, library data is assinged from random boot-strap sampling
#' @param numSamples Integer, which is the boot-strap sampling size
#'      (works when randomLibs = TRUE)
#' @param replace Logical
#'      if TRUE, random sampling can select same data (works when randomLibs = TRUE)
#' @param seed Integer, which sets the RND seed (works when randomLibs = TRUE)
#' @export
xmap <- function(df, ref, tar, lib = matrix(c(1, nrow(df)), nrow = 1), E, Tp = 0,
                 libSize = 10, randomLibs = TRUE, numSamples = 100, replace = FALSE, seed = NULL) {
    ## Check whether lib is provided appropriately
    if (!is.matrix(lib)) lib <- matrix(lib, nrow = 1)
    if (ncol(lib) != 2) stop("Inappropriate lib!")

    X <- df %>% gen_emat(cols = rep(ref, E), lags = 0:-(E - 1), knot = lib)
    y <- df %>% gen_emat(cols = tar, lags = Tp, knot = lib)
    dmat <- get_ned(X)
    idx <- cbind(y, X) %>% complete.cases() %>% which()

    if (randomLibs) {
        set.seed(seed)
        if (!replace) libSize <- min(length(idx), libSize)

        op <- foreach(i = 1:numSamples, .combine = rbind) %do% {
            idx_i <- sample(idx, libSize, replace = replace)
            y_ <- Smplx(X, y, dmat, idx_i - 1, idx - 1, E + 1)
            tibble(ref = ref, tar = tar, E = E, Tp = Tp, libSize = libSize,
                   rho = get_rho(y[idx], y_), mae = get_mae(y[idx], y_), rmse = get_rmse(y[idx], y_))
        }
    } else {
        y_ <- Smplx(X, y, dmat, idx - 1, idx - 1, E + 1)
        op <- tibble(ref = ref, tar = tar, E = E, Tp = Tp, libSize = libSize,
                     rho = get_rho(y[idx], y_), mae = get_mae(y[idx], y_), rmse = get_rmse(y[idx], y_))
    }

    return(op)
}

#' @export
func_smplx <- function(X, y, dmat, idx_l, idx_p, nns) Smplx(X, y, dmat, idx_l, idx_p, nns);

#' Causality test with CCM (Convergent Cross Mapping)
#'
#' \code{ccm_test} returns the (tibble) result of causality test with CCM
#' Algorithm is adopted from Sugihara et al. 2012, Science
#'
#' @inheritParams xmap
#' @param libMin Numeric, which setsthe minumum Library size
#' @param surr A list that consists of the surrogate data of the ref's time series
#' @export
ccm_test <- function(df, ref, tar, lib = matrix(c(1, nrow(df)), nrow = 1), E, Tp = 0,
                     libMin = E + 2, numSamples = 100, seed = NULL, surr = NULL) {
    ## Check whether lib is provided appropriately
    if (!is.matrix(lib)) lib <- matrix(lib, nrow = 1)
    if (ncol(lib) != 2) stop("Inappropriate lib!")

    ## CCM analysis with original data
    max_op <- xmap(df, ref, tar, lib, E, Tp, randomLibs = FALSE) %>% mutate(data = "max")
    min_op <- xmap(df, ref, tar, lib, E, Tp, libSize = libMin, numSamples = numSamples) %>% mutate(data = "min")
    op <- rbind(max_op, min_op) %>% select(ref, tar, E, Tp, data, rho, mae, rmse)

    ## CCM analysis with surrogate data
    if (!is.null(surr)) {
        sur_op <- foreach(i = 1:ncol(surr), .combine = rbind) %do% {
            df %>% mutate(!!ref := surr[, i]) %>%
                xmap(ref, tar, lib, E, Tp, randomLibs = FALSE) %>%
                mutate(data = "surr") %>% select(ref, tar, E, Tp, data, rho, mae, rmse)
        }

        op <- rbind(op, sur_op)
    }

    return(op)
}

#' Grid search of best embedding dimension with simplex projection
#'
#' \code{find_best_dim} returns the best embedding dimension
#'
#' @param df vector, matrix, data.frame or tibble,
#'      which contains time series to be analyzed
#' @param cols Numeric or character vector,
#'      which specifies the df's column for the analysis
#' @param lib Vector or 2-column matrix, which sets the knot position in library data
#' @param pred Vector or 2-column matrix, which sets the knot position in test data
#' @param range Integer vector, which sets the range of dimension search
#' @param Tps Integer vector, which sets the prediction time horizon for each variable
#' @param criterion String, which switch the cost function (Rho, MAE & RMSE)
#'      for the dimension search (default = "rmse")
#' @param both Logical, which switches the back/forward prediction or
#'      forward only for the dimension search (defualt = TRUE)
#' @param threadNo Integer, which sets the thread number of parallel computing
#' @export
find_best_dim <- function(df, cols, lib = matrix(c(1, nrow(df)), nrow = 1), pred = NULL, range = 1:10,
                          Tps = rep(1, length(cols)), criterion = "rmse", both = TRUE, threadNo = detectCores()) {
    ## Check whether lib is provided appropriately
    if (!is.matrix(lib)) lib <- matrix(lib, nrow = 1)
    if (ncol(lib) != 2) stop("Inappropriate lib!")

    if (is.null(pred)) {
        knot <- lib
        pred <- lib
    } else {
        knot <- rbind(lib, pred)
    }

    ## Implementation of Simplex Projection algorithm
    smplx <- function(col, E, Tp) {
        X <- df %>% gen_emat(cols = rep(col, E), lags = 0:(sign(Tp) * -(E - 1)), knot = knot)
        y <- df %>% gen_emat(cols = col, lags = Tp, knot = knot)
        dmat <- get_ned(X)
        y_ <- rep(NA, nrow(X))

        idx_all <- cbind(y, X) %>% complete.cases() %>% which()
        idx_l <- gen_valid_idx(lib, idx_all)
        idx_p <- gen_valid_idx(pred, idx_all)
        
        tryCatch({
            y_[idx_p] <- Smplx(X, y, dmat, idx_l - 1, idx_p - 1, E + 1)
            return(tibble(E = E, rho = -get_rho(y, y_), mae = get_mae(y, y_), rmse = get_rmse(y, y_)))
        }, error = function(e) return(tibble(E = E, rho = NA, mae = NA, rmse = NA)))
    }

    ## Set environment for parallel computing
    cl <- makeCluster(threadNo, type = "PSOCK")
    registerDoParallel(cl)
    on.exit(stopCluster(cl))

    ## Find best embedding dimension for each variable
    op <- foreach(i = 1:length(cols), .combine = rbind) %:% foreach(E = range, .combine = rbind) %dopar% {
        fwd <- smplx(cols[i], E, Tps[i])

        if (both) {
            bwd <- smplx(cols[i], E, -Tps[i])
            fwd <- (fwd + bwd) / 2
        }

        fwd %>% mutate(var = cols[i], Tp = Tps[i])
    }

    op <- op %>% group_by(var) %>% arrange(!!!rlang::syms(criterion)) %>% slice(1) %>%
        mutate(rho = -rho) %>% ungroup() %>% select(var, Tp, E, rho, mae, rmse)
    return(op)
}

#' S-map with elastic net (a.k.a. Regularized S-map)
#'
#' \code{smap_net} returns the (tibble) result of Regularized S-map
#' Algorithm is adopted from Cenci et al. (2019), MEE
#'
#' @inheritParams find_best_dim
#' @param X A matrix, which is constructed by the function \code{\link{gen_emat}}
#' @param col Integer, which specifies the position of prediction target in X
#' @param Tp Integer, which is the prediction-time horizon
#' @param range Numeric vector, which sets the range of nonlinear parameter theta for S-map
#' @param seed Integer, which sets the RNG seed
#' @param s String, which switches the criterion of the best theta in elastic net
#'      (for more detail, see original function \code{\link[glmnet]{cv.glmnet}})
#' @param lambda Numeric vector, which sets the sequence of penalty strength for glmnet
#'      if NULL (default), then the adaptive search with cv.glmnet is adopted
#'      (for more detail, see original function \code{\link[glmnet]{cv.glmnet}})
#' @param alpha Numeric, which controls the Elastic net parameter
#'      \code{alpha = 0} is the ridge penalty and \code{alpha = 1} is the lasso penalty
#'      (for more detail, see original function \code{\link[glmnet]{cv.glmnet}})
#' @export
smap_net <- function(X, col = 1, lib = matrix(c(1, nrow(X)), nrow = 1), pred = NULL, Tp = 1, threadNo = detectCores(),
                     range = seq(0, 10, 1), seed = NULL, criterion = "rmse", s = "lambda.1se", lambda = NULL, alpha = 0.0) {
    ## Check whether lib is provided appropriately
    if (!is.matrix(lib)) lib <- matrix(lib, nrow = 1)
    if (ncol(lib) != 2) stop("Inappropriate lib style!")

    if (is.null(pred)) {
        knot <- lib
        pred <- lib
    } else {
        knot <- rbind(lib, pred)
    }

    dmat <- get_ned(X);

    ## Preparation for Elastic-net analysis with R package 'glmnet'
    dim <- dim(X)
    cnames <- str_c("C", c(1:dim[2], 0))
    y <- X %>% gen_emat(cols = col, lags = Tp, knot = knot)
    idx_all <- apply(dmat, 1, function(r) {sum(is.na(r)) != dim[1] - 1})
    idx_all <- which(complete.cases(y) & idx_all)
    idx_l <- gen_valid_idx(lib, idx_all)
    idx_p <- gen_valid_idx(pred, idx_all)

    dbar <- rowMeans(dmat[, idx_l], na.rm = TRUE)

    ## Set environment for parallel computing
    cl <- makeCluster(threadNo, type = "PSOCK")
    registerDoParallel(cl)
    on.exit(stopCluster(cl))

    set.seed(seed)
    seeds <- sample(32768, dim[1], replace = TRUE)

    ## grid-search for best theta
    op <- foreach(theta = range, .combine = rbind) %dopar% {
        tryCatch({
            ## Make weight matrix
            coef <- matrix(NA, nrow = dim[1], ncol = dim[2] + 1)
            wmat <- exp(-theta * dmat / dbar)
            diag(wmat) <- 0

            ## Sequential estimation of Jacobian
            coef[idx_p, ] <- foreach(t = idx_p, .combine = rbind) %do% {
                set.seed(seeds[t])
                fit <- cv.glmnet(x = X[idx_l, ], y = y[idx_l], weights = wmat[t, idx_l], lambda = lambda, alpha = alpha)
                coef(fit, s = s) %>% as.numeric() %>% .[c(2:(dim[2] + 1), 1)]
            }

            y_ <- getPred(coef, cbind(X, 1))
            output <- tibble(obs = as.numeric(y), pred = y_) %>% list()
            coef <- coef %>% as.data.frame() %>% setNames(cnames) %>% as_tibble() %>% list()

            tibble(theta = theta, rho = -get_rho(y, y_), mae = get_mae(y, y_), rmse = get_rmse(y, y_),
                   output = output, coef = coef)
        }, error = function(e) {
            output <- tibble(obs = as.numeric(y), pred = NA) %>% list()
            coef <- coef %>% as.data.frame() %>% setNames(cnames) %>% as_tibble() %>% list()
            tibble(theta = theta, rho = NA, mae = NA, rmse = NA, output = output, coef = coef)
        })
    } %>% arrange(!!!rlang::syms(criterion)) %>% dplyr::slice_head() %>% mutate(rho = -rho)

    return(op)
}

