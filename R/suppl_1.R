####
#### suppl_1.R: R code for Figure S1-2 in "Flexibility of species intraction and ecological stability"
#### analysis of coupled multi-species Logistic Map
####
#### ver 1.0.0: initially written on 20210621 by K.Kawats
#### ver 1.1.0: updated on 20210718 by K.Kawatsu
#### ver 2.0.0: updated on 20220808 by K.Kawatsu

### load source code
source("R/functions.R")

### Figure S1
r <- c(0.5, -1.5)
A <- matrix(c(0.0, -1.4,
              0.8,  0.0), nrow = 2, byrow = TRUE)

curr1 <- c(2, 2)
curr2 <- c(-r[2] / A[2, 1], r[1] / (-A[1, 2] * (-r[2] / A[2, 1])^-0.25) + 0.03)
df1 <- gen_dl_model(state = curr1, len = 7500, r = r, A = A, beta =  0.25, h = 0.01) %>% rename(V = x1, P = x2)
df2 <- gen_dl_model(state = curr2, len = 7500, r = r, A = A, beta = -0.25, h = 0.01) %>% rename(V = x1, P = x2)

### make Figure S1
gp_S1a <- df1 %>% pivot_longer(-time, names_to = "species", values_to = "density") %>%
    ggplot(aes(x = time, y = density, color = species)) +
    geom_line(size = 0.25) +
    scale_color_brewer(palette = "Set1") + xlab("Time") + ylab("Density") + labs(tag = expression((italic(a))))

gp_S1b <- tibble(V = seq(0, 5, 0.01), A12 = -1.4 * V^0.25) %>% ggplot(aes(x = V, y = A12)) +
    geom_hline(yintercept = 0, size = 0.25, linetype = 2) + geom_line() +
    scale_y_continuous(breaks = seq(-5, 0, 1), limit = c(-5, 0)) +
    xlab(expression(italic(V))) + ylab(expression(italic(A[VP])))

gp_S1c <- df2 %>% pivot_longer(-time, names_to = "species", values_to = "density") %>%
    ggplot(aes(x = time, y = density, color = species)) +
    geom_line(size = 0.25) +
    scale_color_brewer(palette = "Set1") + xlab("Time") + ylab("Density") + labs(tag = expression((italic(b))))

gp_S1d <- tibble(V = seq(0, 5, 0.01), A12 = -1.4 * V^-0.25) %>% ggplot(aes(x = V, y = A12)) +
    geom_hline(yintercept = 0, size = 0.25, linetype = 2) + geom_line() +
    scale_y_continuous(breaks = seq(-5, 0, 1), limit = c(-5, 0)) +
    xlab(expression(italic(V))) + ylab(expression(italic(A[VP])))

fig_S1 <- (gp_S1a + gp_S1b + plot_layout(nrow = 1, width = c(2, 1))) /
    (gp_S1c + gp_S1d + plot_layout(nrow = 1, width = c(2, 1))) & theme_st(c(1, 1), c(1, 1))
ggsave("fig/fig_S1.eps", fig_S1, device = cairo_ps, fallback_resolution = 600, family = "Times", width = 10, height = 8, units = "cm")

### Figure S2a
seed <- 333
size <- 5
p_noise <- 0.00
o_noise <- 0.00
prop <- 0.5

### demonstration of estimating per capita interaction strengths with coupled multi-species Logistic Map
system.time(suppl1 <- foreach(b12 = seq(-0.75, 0.75, 0.75), .combine = rbind) %do% {
    ## set demographic parameters
    set.seed(seed)
    r <- rnorm(size, 3.65, 0.05)
    A <- matrix(rnorm(size^2, 0, 0.05) * rbinom(size^2, 1, prop), nrow = size)
    diag(A) <- -r
    A[1, 2] <- -0.3

    ## generate time-series data
    df <- gen_dl_map(rep(0.01, size), 250, r = r, A = A, beta = b12, minTh = 0, maxTh = Inf, p_noise)

    if (nrow(df) > 200) {
        df <- df %>% filter(time > 50)
        df_ <- df %>% mutate(across(.cols = starts_with("x"), .fns = ~{.x * runif(length(.x), 1 - o_noise, 1 + o_noise)}))
        S <- foreach(i = which(A[1, ] != 0), .combine = c) %do% str_c("x", i);
        E <- df_ %>% find_best_dim(cols = "x1") %>% pull(E)

        emat <- df_ %>% gen_emat(cols = S, lags = rep(0, length(S)))
        if (E > length(S)) emat <- foreach(i = 1:(E - length(S)), .combine = cbind, .init = emat) %do% lag(df_$x1, i);

        ## S-map analysis for indirect estimation
        dmat <- NULL; # dmat <- df_ %>% get_mvd(cols = S, E = max(length(S), E), Tp = 1, lmax = 4)
        op1 <- smap_net(emat = emat, Tp = 1, range = seq(10, 40, 0.2), dmat = dmat, seed = seed)

        ## S-map analysis for direct estimation
        emat[, 2:length(S)] <- emat[, 1] * emat[, 2:length(S)]
        # dmat <- df_ %>% mutate(across(.cols = all_of(S[-1]), .fns = ~{.x * x1})) %>% get_mvd(cols = S, E = max(length(S), E), Tp = 1, lmax = 4)
        # op2 <- smap_net(emat = emat, Tp = 1, range = seq(10, 40, 0.2), dmat = dmat, seed = seed)
        op2 <- find_best_theta(emat = emat, iterLim = 100, tsLength = 100, gSeed = seed, sigma = 1.0, minTemp = 1e-10, maxTemp = 1e-2)

        tibble(a12 = A[1, 2], b12 = b12, x1 = df$x1, indirect = op1$coef[[1]]$C2 / emat[, 1], direct = op2$coef[[1]]$C2)
    } else {
        tibble(a12 = A[1, 2], b12 = b12, x1 = df$x1, indirect = NA, direct = NA)
    }
})

suppl1 %>% write_csv("data/suppl1.csv")

### Figure S2b
seed <- 333
size <- 5
prop <- 0.5
params <- expand.grid(iter = 1:100, o_noise = seq(0.001, 0.051, length = 21), p_noise = seq(0.001, 0.051, length = 21))
seeds <- sample(32768, nrow(params), replace = TRUE)

system.time(suppl2 <- foreach(i = 1:nrow(params), .combine = rbind) %do% {
    ## set demographic parameters
    set.seed(seeds[i])
    r <- rnorm(size, 3.65, 0.05)
    A <- matrix(rnorm(size^2, 0, 0.05) * rbinom(size^2, 1, prop), nrow = size)
    diag(A) <- -r
    A[1, 2] <- -0.3
    b12 <- rnorm(1, 0.0, 0.2)
    p_noise <- params$p_noise[i]
    o_noise <- params$o_noise[i]

    ## generate time-series data
    df <- gen_dl_map(rep(0.01, size), 250, r = r, A = A, beta = b12, minTh = 0, maxTh = Inf, p_noise)

    tryCatch({
        df_ <- df %>% filter(time > 50) %>%
            mutate(across(.cols = starts_with("x"), .fns = ~{.x * runif(length(.x), 1 - o_noise, 1 + o_noise)}))
        S <- foreach(i = which(A[1, ] != 0), .combine = c) %do% str_c("x", i);
        E <- df_ %>% find_best_dim(cols = "x1") %>% pull(E)

        emat <- df_ %>% gen_emat(cols = S, lags = rep(0, length(S)))
        if (E > length(S)) emat <- emat %>% cbind(gen_emat(df_, cols = rep("x1", E - length(S) + 1), 0:-(E - length(S))))

        ## S-map analysis for indirect estimation
        op1 <- smap_net(emat = emat, Tp = 1, range = seq(0, 50, 1), seed = seeds[i]) %>% pull(coef) %>% .[[1]]

        ## S-map analysis for direct estimation
        emat[, 2:length(S)] <- emat[, 1] * emat[, 2:length(S)]
        op2 <- smap_net(emat = emat, Tp = 1, range = seq(0, 50, 1), seed = seeds[i]) %>% pull(coef) %>%.[[1]]

        coef_t <- df %>% filter(time > 50) %>% mutate(x1 = A[1, 2] * x1^b12) %>% pull(x1)
        tibble(b12 = b12, p_noise = p_noise, o_noise = o_noise, indirect = get_rmse(coef_t, op1$C2 / emat[, 1]), direct = get_rmse(coef_t, op2$C2))
    }, error = function(e) {tibble(b12 = b12, p_noise = p_noise, o_noise = o_noise, indirect = NA, direct = NA)})
})

suppl2 %>% write_csv("data/suppl2.csv")

### make Figure S2
suppl1 <- read_csv("data/suppl1.csv")

gp_S2a <- suppl1 %>%
    pivot_longer(cols = c(4, 5), names_to = "method", values_to = "coef") %>%
    ggplot(aes(x = x1, color = factor(b12))) +
    facet_wrap(~method) +
    geom_point(aes(y = coef), size = 0.25) +
    geom_line(aes(y = a12 * x1^b12)) +
    scale_color_brewer(palette = "Set2") +
    ylim(-1.2, 1.2) +
    xlab(expression(paste("Density of species ", italic(x)[1](italic(t))))) +
    ylab(expression(paste("Per capita coefficients ", italic(f)[12]))) +
    labs(tag = expression((italic(a)))) +
    theme_st()

suppl2 <- read_csv("data/suppl2.csv")

gp_S2b <- suppl2 %>% group_by(p_noise, o_noise) %>%
    group_modify(~summarize(., across(.cols = everything(), .fns = ~mean(.x, na.rm = TRUE)))) %>%
    ungroup() %>% pivot_longer(cols = c(direct, indirect), names_to = "method", values_to = "RMSE") %>%
    ggplot(aes(x = p_noise, y = o_noise, fill = RMSE)) + facet_wrap(~method) +
    geom_raster(interpolate = TRUE) +
    scale_fill_viridis(direction = -1) +
    scale_x_continuous(expand = rep(5e-4, 2)) +
    scale_y_continuous(expand = rep(5e-4, 2)) +
    xlab(expression(paste("Process noise (", sigma[p], ")"))) +
    ylab(expression(paste("Observation error (", sigma[o], ")"))) +
    labs(tag = expression((italic(b)))) +
    theme_st()

fig_S2 <- gp_S2a / gp_S2b
ggsave("fig/fig_S2.eps", fig_S2, device = cairo_ps, fallback_resolution = 600, family = "Times", width = 10, height = 10, units = "cm")

