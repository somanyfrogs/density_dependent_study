####
#### R function file for "Flexibility of species interaction and ecological stability"
#### func_surrogate.R: contains functions related to surrogate-data generation
####
#### ver 0.0.1: Initially written on 20200227 by K.Kawatsu
#### ver 0.0.2: Updated on 20200313
#### ver 0.1.0: Updated on 20210706
#### ver 0.1.1: Updated on 20211015
#### ver 0.2.0: Updated on 20211118
####

#' Calculate recurrence matrix
#'
#' \code{gen_rec_mat} returns (binary) 0-1 matrix based on distance matrix
#'
#' @param dmat A distance matrix
#' @param s Numeric, which sets neighbourhood radius based on the proportion 's' of dmat
gen_rec_mat <- function(dmat, s = 0.100) {
    rad <- quantile(dmat, probs = s, na.rm = TRUE)
    ifelse(dmat <= rad, 1, 0)
}

#' Find the 'twin' status in embedding matrix
#'
#' \code{find_twins} returns a list preserving dynamical 'twins' for each time point
#' Algorithm is adopted from Thiel et al. (2006), Europhys Lett and Ushio et al. (2018)
#' Nature, with an extension that includes NN (nearest neighbour) twin method
#'
#' @param dmat A distance matrix
#' @param twin_op Character, which can be "RM" or "NN"
#'      if "RM", twin search is based on recurrence matrix,
#'      otherwise is based on NN twin method
#' @param time_idx Vector of time index
#' @param phase_lock Logical, if TRUE, twin search is limited to the points having
#'      time index same with the target (a.k.a. phase-locked twin surrogate, Ushio et al. 2018)
#' @param period Integer, the phase period (only works when phase_lock == TRUE)
#' @param lim, Integer, the minimum No. of twins (only works when twin_op == "RM")
#' @param range Vector, the range of threshold in recurrence matrix
#'      (only works when twin_op == "RM", see also \code{\link{gen_rec_mat}})
#' @param nn Integer, the number of nearest neighbour points
#'      (only works when twin_op != "RM")
#' @export
find_twins <- function(dmat, twin_op = "RM", time_idx, phase_lock = TRUE,
                       period = 12, lim = 10, range, nn) {
    twins <- as.list(1:nrow(dmat))
    idx <- which(apply(dmat, 1, function(r) !all(is.na(r))))

    if (twin_op == "RM") {          ## twin search follows to recurrence matrix algorithm
        for (s in range) {
            rmat <- gen_rec_mat(dmat, s = s)

            for (i in idx) {
                tmp <- NULL
                for (j in idx) if (all(rmat[i, ] == rmat[j, ], na.rm = TRUE)) tmp <- c(tmp, j)
                if (phase_lock) tmp <- tmp[(time_idx[tmp] - time_idx[i]) %% period == 0]
                twins[[i]] <- tmp
            }

            tws_no <- length(unlist(twins)) - length(twins)
            if (tws_no >= lim) break;
        }
    } else if (twin_op == "NN") {   ## twin search follows to nearest neighbour algorithm
        for (i in idx) {
            tmp <- order(dmat[i, ], na.last = NA)
            twin_end <- nn
            if (phase_lock) tmp <- tmp[(time_idx[tmp] - tmp_idx[i]) %% period == 0]

            if (length(tmp) >= twin_end) {
                repeat {
                    if (dmat[i, tmp[twin_end]] < dmat[i, tmp[twin_end + 1]]) break;
                    twin_end <- twin_end + 1
                }

                tmp <- tmp[1:twin_end]
            }
        }
    } else {
        warning("Inappropriate method for 'twin_op', NULL twin list returns")
    }

    return(twins)
}

#' Calculate twin probabilities
#'
#' \code{twin_probs} returns a list of swapping probabilities of twins for each state
#'
#' @inheritParams find_twins
#' @param twins A list, that preserves twins for each state (calculated by \code{\link{find_twins}})
#' @param theta Numeric, the locality parameter (works if twin_op == "NN")
#' @export
twin_probs <- function(twins, dmat, theta = 1, twin_op) {
    probs <- NULL

    if (twin_op == "RM") {
        probs <- lapply(twins, function(v) {rep(1 / length(v), length(v))})
    } else if (twin_op == "NN") {
        probs <- foreach(i = 1:length(twins)) %do% {
            di <- dmat[i, twins[i]]
            wi <- exp(-theta * di) / mean(di)
            wi / sum(wi)
        }
    } else {
        warning("Inappropriate method for 'twin_op', NULL probability list returns")
    }

    return(probs)
}

#' Transition procedure to the next point from twins
#'
#' \code{point_nex} proceeds the time series with swapping between original and twin states
#' If valid data does not exist, this function returns 0
#'
#' @inheritParams twin_probs
#' @param curr Integer, the position of current state
#' @param probs List of twin probabilitities (see \code{\link{twin_probs}})
#' @param knot Vector or 2-column matrix, which specifies the knot positions in the data
#' @export
point_nex <- function(curr, twins, knot, probs) {
    if (curr == 0) return(0)
    if (!is.matrix(knot)) knot <- matrix(knot, nrow = 1)
    if (ncol(knot) != 2) stopd("Inappropriate knot!")

    ## replace with 0 if the next point is the end of the knot
    if (length(twins[[curr]]) > 1) {
        nex <- sample(x = twins[[curr]], size = 1, prob = probs[[curr]])
    } else {
        nex <- curr
    }

    nex <- if_else(nex %in% knot[, 2], 0, nex + 1)
    return(nex)
}

#' Surrogate generation with twin-surrogate
#'
#' \code{gen_surr_ts} returns a column-wide piled surrogate data with TS algorithm
#'
#' @inheritParams find_twins
#' @param X data matrix
#' @param col Integer, which specifies the column position of surrogate target variable in 'X'
#' @param knot Vector or 2-column matrix, which specifies the knot positions in the data
#' @param iter Integer, No. of surrogate iteration
#' @param seed, Integer, which sets RNG seed
#' @param theta Numeric, which determines the locality parameter (works if twin_op == 'simplex')
#' @export
gen_surr_ts <- function(X, col, knot = matrix(c(1, nrow(X)), nrow = 1), iter = 100, seed = NULL,
                        twin_op = "RM", time_idx, phase_lock = TRUE, period = 12, lim = 10,
                        range = c(0.125, seq(0.12, 0.05, -0.01), seq(0.15, 0.20, 0.01), 0.04), theta = 1) {
    if (!is.matrix(knot)) knot <- matrix(knot, nrow = 1)
    if (ncol(knot) != 2) stop("Inappropriate knot style!")

    dmat <- get_ned(X)
    idx <- X %>% complete.cases() %>% which()
    twins <- find_twins(dmat, twin_op = twin_op, time_idx = time_idx, phase_lock = phase_lock,
                        period = period, lim = lim, range = range)

    if (max(sapply(twins, length)) == 1) {
        twins <- find_twins(dmat, twin_op = "NN", time_idx = time_idx,
                            phase_lock = phase_lock, period = period, nn = ncol(X) + 1)
        probs <- twin_probs(twins, dmat, theta = theta, twin_op = "NN")
    } else {
        probs <- twin_probs(twins, dmat, theta = theta, twin_op = twin_op)
    }

    set.seed(seed)
    surr <- foreach(i = 1:nrow(knot), .combine = rbind) %do% {
        tmp <- foreach(j = 1:iter, .combine = cbind) %do% {
            if (phase_lock) {
                ids <- sample(which((time_idx - time_idx[knot[i, 1]]) %% period == 0), 1)
            } else {
                ids <- sample(1:nrow(X)[-knot[, 2]], 1)
            }

            for (k in 2:(knot[i, 2] - knot[i, 1] + 1)) ids <- c(ids, point_nex(ids[k - 1], twins, knot, probs))

            ## return the surrogate if it reaches to the original data length
            if (all(ids != 0)) X[ids, col]
        }

        ## if cannot yield surrogate, original data is used instad
        ## also if the size of surrogate data does no reach the iteration number
        ## then meet the condition by sampling the surrpgate column with replace TRUE
        if (is.null(tmp)) tmp <- X[knot[i, 1]:knot[i, 2], rep(col, iter)];
        if (iter == 1) tmp <- matrix(tmp, ncol = 1)
        if (ncol(tmp) < iter) tmp <- tmp[, sample(1:ncol(tmp), iter, replace = TRUE)]

        tmp
    }

    return(surr)
}
