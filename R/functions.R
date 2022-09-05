####
#### R functions file for "Flexibility of species interactino and ecological stability"
#### functions.R: contains functions for the downstream analyses included in the paper
####
#### ver 0.0.1: Initially written on 20210319 by K.Kawatsu
#### ver 0.0.2: Updated on 20210326
#### ver 0.0.3: Updated on 20210329
#### ver 0.0.4: Updated on 20210402
#### ver 0.1.0: Updated on 20210617
#### ver 0.1.1: Updated on 20210621
#### ver 0.2.0: Updated on 20210702
#### ver 0.3.0: Updated on 20211117
#### ver 0.4.0: Updated on 20220808

### load dependent libraries
library(tidyverse)
library(broom)
library(broom.mixed)
library(foreach)
library(doParallel)
library(lme4)
library(glmnet)
library(rpkg)
library(patchwork)

### Make density-dependent coupled logistic map
### r, A: intrinsic growth rate and interaction matrix
gen_dl_model <- function(state, len, r, A, beta, h = 0.01) {
    fun <- function(x) {
        A_ <- A
        A_[1, 2] <- A_[1, 2] * x[1]^beta
        return(as.numeric((r + A_ %*% x) * x))
    }

    tmp <- foreach(t = 1:len, .combine = rbind) %do% {
        k1 <- h * fun(state)
        k2 <- h * fun(state + k1 / 2)
        k3 <- h * fun(state + k2 / 2)
        k4 <- h * fun(state + k3)

        state <- (state + (k1 + 2 * k2 + 2 * k3 + k4) / 6)
    }

    cnames <- c("time", str_c("x", 1:length(state)))
    tmp <- data.frame(1:len, tmp) %>% setNames(cnames) %>% as_tibble()
    return(tmp)
}

### Make density-dependent coupled logistic map
### r, A: intrinsic growth rate and interaction matrix
gen_dl_map <- function(state, len, r, A, beta, minTh, maxTh, p_noise) {
    tmp <-
        foreach(t = 1:len, .combine = rbind, .init = state) %:%
        when(all(between(r, minTh, maxTh))) %do% {
            A_ <- A
            A_[1, 2] <- A_[1, 2] * state[1]^beta
            state <- as.numeric((r + A_ %*% state) * state) * runif(length(state), 1 - p_noise, 1 + p_noise)
        }

    cnames <- c("time", str_c("x", 1:length(state)))
    tmp <- data.frame(1:(len + 1), tmp) %>% setNames(cnames) %>% as_tibble()
    return(tmp)
}

### Make state-structured 2-prey-1-predator model
gen_pp_map <- function(state, len, rc, rp, C, P, beta) {
    tmp <- foreach(i = 1:len, .combine = rbind) %do% {
        nex <- as.numeric((rc + C %*% state) * state)

        P_ <- P
        P_[1, 3] <- P[1, 3] * nex[1]^beta
        P_[2, 3] <- P[2, 3] * nex[2]^beta
        state <- as.numeric((rp + P_ %*% nex) * nex)
    }

    cnames <- c("time", str_c("x", 1:length(state)))
    tmp <- data.frame(1:len, tmp) %>% setNames(cnames) %>% as_tibble()
    return(tmp)
}

### Shuffle target interaction coefficients
shuffle_coef <- function(C_orig, dens, rho, sigma = NULL, seed = NULL) {
    set.seed(seed)

    C_hat <- loess(C_orig ~ dens, span = 1) %>% predict(dens)
    mu <- mean(C_hat, na.rm = TRUE)
    if (is.null(sigma)) sigma <- sd(C_hat, na.rm = TRUE)
    C_surr <- C_orig - C_hat + mu + rho * (C_hat - mu) + sqrt(1 - rho^2) * rnorm(length(C_hat), 0, sigma)

    return(C_surr)
}

### ggplot theme for figures
theme_st <- function(just = NULL, pos = NULL, hl = 5.5) {
    mytheme <-
        theme_light() +
        theme(axis.title = element_text(size = 8),
              axis.text = element_text(size = 6),
              legend.background = element_blank(),
              legend.title = element_blank(),
              legend.text = element_text(size = 6),
              legend.key.size = unit(0.25, "cm"),
              plot.tag = element_text(size = 8),
              strip.text = element_text(size = 6, color = "black",
                                        margin = margin(0.5 * hl, 0.5 * hl, 0.5 * hl, 0.5 * hl)))

    if (!is.null(just)) mytheme <- mytheme + theme(legend.justification = just)
    if (!is.null(pos)) mytheme <- mytheme + theme(legend.position = pos)

    return(mytheme)
}
