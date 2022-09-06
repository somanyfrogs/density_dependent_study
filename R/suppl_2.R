####
#### suppl_2.R: R code for Figure S3-4 in "Temporal flexibility of species interactions and ecological stability"
#### analysis of maiduru fish monitoring data
####
#### ver 0.0.1: Initially written on 20210718 by K.Kawatsu
#### ver 0.0.2: Updated on 20210723
#### ver 0.1.0: Updated on 20211120
#### ver 0.2.0: Updated on 20211215
#### ver 0.3.0: Updated on 20220808
#### ver 0.3.1: Updated on 20220906

### load source code
source("R/functions.R")

### load time-series dataset of maiduru fish monitoring
### prepare 'data' directory containing maiduru_obs.csv under your current workspace
fish_orig <- read_csv("data/maiduru_obs.csv")

### Figure S3
### make monthly averaged data
fish <- fish_orig %>%
    mutate(date = as.Date(date, origin = as.Date("2002-01-01")), year = str_sub(date, 1, 4), month = str_sub(date, 6, 7)) %>%
    select(-date) %>% group_by(year, month) %>%
    group_modify(~summarize(., across(.cols = everything(), .fns = ~mean(.x, na.rm = TRUE)))) %>% ungroup()

### species selection
fish_sum <- read_csv("data/fish_sum.csv")
col_alpha <- fish_sum %>% filter(Obs >= nrow(fish) * 0.5) %>% pull(Sp)
col_beta <- fish_sum %>% filter(Obs >= nrow(fish) * 0.4) %>% pull(Sp) %>% .[-(1:13)]
fish <- fish %>% mutate(time = 1:nrow(.), across(.cols = all_of(c(col_alpha, col_beta)), .fns = normalize))

### Figure S3
gp_S3a <- fish %>% select(all_of(c("time", col_alpha))) %>%
    pivot_longer(-time, names_to = "Species", values_to = "Obs") %>%
    mutate(Species = as.numeric(factor(Species, levels = col_alpha))) %>%
    ggplot(aes(x = time, y = Obs)) +
    facet_wrap(~Species, ncol = 3) +
    geom_path(size = 0.25) +
    scale_x_continuous(breaks = seq(1, 201, 48), labels = seq(2002, 2018, 4)) +
    xlab("Year") + ylab("Observation (normalized)") + labs(tag = expression(bold(A)))

gp_S3b <- fish %>% select(all_of(c("time", col_beta))) %>%
    pivot_longer(-time, names_to = "Species", values_to = "Obs") %>%
    mutate(Species = as.numeric(factor(Species, levels = col_beta)) + 13) %>%
    ggplot(aes(x = time, y = Obs)) +
    facet_wrap(~Species, ncol = 3) +
    geom_path(size = 0.25) +
    scale_x_continuous(breaks = seq(1, 201, 48), labels = seq(2002, 2018, 4)) +
    xlab("Year") + ylab("Observation (normalized)") + labs(tag = expression(bold(B)))

fig_S3 <- (gp_S3a + gp_S3b + plot_layout(ncol = 1, height = c(5, 2))) & theme_st()

## save figure S3 in 'fig' directory under your current workspace
ggsave("fig/fig_S3.eps", fig_S3, device = cairo_ps, fallback_resolution = 1200, family = "Helvetica", width = 15, height = 20, units = "cm")

### Figure S4
### Jacobian-estimation analysis and simulation
### load time-series dataset of maiduru fish monitoring
fish_orig <- read_csv("data/maiduru_obs.csv")

### make monthly averaged data
fish <- fish_orig %>%
    mutate(date = as.Date(date, origin = as.Date("2002-01-01")), year = str_sub(date, 1, 4), month = str_sub(date, 6, 7)) %>%
    select(-date) %>% group_by(year, month) %>%
    group_modify(~summarize(., across(.cols = everything(), .fns = ~mean(.x, na.rm = TRUE)))) %>% ungroup()

### Species selection (observed more than 40% of total data length (201))
ccm_op <- read_csv("data/ccm_op.csv")
cols <- read_csv("data/fish_sum.csv") %>% filter(Obs > 201 * 0.4) %>% pull(Sp)
fish <- fish %>% mutate(across(.cols = all_of(cols), .fns = ~normalize(.x, method = "0-1")))

seed <- 123
rhos <- seq(0.09, 1, 0.1)
iter <- 100

seeds <- sample(32768, iter * length(rhos), replace = TRUE)

system.time(smap_op_02 <- foreach(rec = cols, .combine = rbind) %do% {
    ccm_r <- ccm_op %>% filter(ref == rec & tar %in% cols) %>% filter(p_min <= 0.05 & p_sur <= 0.05)

    tryCatch({
        ccm_r <- ccm_r %>% group_by(tar) %>% group_modify(~filter(., rho == max(rho))) %>% ungroup()
        dons <- ccm_r %>% pull(tar)
        S <- dons %>% length()
        E <- ccm_r %>% pull(E) %>% unique()

        ## make product embedding
        emat <- fish %>% gen_emat(cols = c(dons, rep(rec, max(1, E - S))), lags = c(pull(ccm_r, Tp) + 1, 0:-max(0, E - (S + 1))))
        emat[, 1:S] <- emat[, 1:S] * emat[, 1 + S]
        op <- emat %>% smap_net(Tp = 1, range = seq(0, 20, 1), s = "lambda.min", seed = seed, alpha = 0.05)

        ## density-dependent analysis
        coef <- op$coef[[1]]
        c_means <- coef %>% summarize(across(.cols = everything(), .fns = ~mean(.x, na.rm = TRUE))) %>% as.numeric()
        c_sds <- coef %>% summarize(across(.cols = everything(), .fns = ~sd(.x, na.rm = TRUE))) %>% as.numeric()
        cnames <- names(coef)[1:S]
        cv_o <- coef %>% as.matrix() %>% get_pred(cbind(emat, 1)) %>% get_cv()

        op_dda <- foreach(i = 1:S, .combine = rbind) %do% {
            glm_op <- broom::tidy(glm(pull(coef, i) ~ emat[, S + 1]))

            tibble(rec = !!rec, theta = op$theta, rho = op$rho, don = dons[i], c_mean = c_means[i], c_sd = c_sds[i],
                   dd_est = glm_op$estimate[2], dd_sd = glm_op$std.error[2],
                   p_value1 = t.test(pull(coef, i), na.action = "na.or.complete")$p.value,
                   p_value2 = glm_op$p.value[2])
        }

        cl <- makeCluster(detectCores(), type = "PSOCK")
        registerDoParallel(cl)
        on.exit(stopCluster(cl))

        ## sensitivity analysis
        op_sim <- foreach(rho = rhos, .combine = rbind) %:%
                  foreach(j = 1:iter, .combine = rbind, .packages = c("foreach", "tidyverse", "rpkg")) %dopar% {
            s <- seeds[iter * (which(rhos == rho) - 1) + j]
            surr <- coef %>% mutate(across(.cols = all_of(cnames), .fns = ~shuffle_coef(.x, dens = emat[, 1 + S], rho = rho, seed = s)))
            cv_all <- surr %>% as.matrix() %>% get_pred(cbind(emat, 1)) %>% get_cv()

            foreach(i = 1:S, .combine = rbind) %do% {
                cv <- coef %>% mutate(!!cnames[i] := pull(surr, cnames[i])) %>% as.matrix() %>% get_pred(cbind(emat, 1)) %>% get_cv()
                tibble(rec = rec, don = dons[i], rho = rho, sensitivity = cv / cv_o)
            } %>% rbind(tibble(rec = rec, don = "all", rho = rho, sensitivity = cv_all / cv_o))
        }

        tibble(rec = rec, dda = list(op_dda), sim = list(op_sim))
    }, error = function(e) {
        op_dda <- tibble(rec = rec, theta = NA, rho = NA, don = cols[cols != rec], c_mean = NA, c_sd = NA,
                         dd_est = NA, dd_sd = NA, p_value1 = NA, p_value2 = NA)
        op_sim <- tibble(rec = rec, don = cols[cols != rec], rho = NA, sensitivity = NA)
        tibble(rec = rec, dda = list(op_dda), sim = list(op_sim))
    })
})

## save density-dependence analysis of 18 fish species in 'data' directory
foreach(i = 1:nrow(smap_op_02), .combine = rbind) %do% {smap_op_02$dda[[i]]} %>% write_csv("data/dda_op_02.csv")
foreach(i = 1:nrow(smap_op_02), .combine = rbind) %do% {smap_op_02$sim[[i]]} %>% write_csv("data/sim_op_02.csv")

### make Figure S4
dda_op <- read_csv("data/dda_op_02.csv")
sim_op <- read_csv("data/sim_op_02.csv")
cols <- read_csv("data/fish_sum.csv") %>% filter(Obs >= 201 * 0.4) %>% pull(Sp)

gp_S4a <- dda_op %>% mutate(rec = factor(rec, levels = rev(cols)), don = factor(don, levels = cols)) %>%
    mutate(sign = case_when(c_mean < 0 & p_value1 <= 0.05 ~ "neg", c_mean > 0 & p_value1 <= 0.05 ~ "pos", TRUE ~ "zero")) %>%
    ggplot(aes(x = don, y = rec, color = sign)) +
    geom_point(aes(size = abs(c_mean))) +
    scale_color_manual(values = c(RColorBrewer::brewer.pal(3, "Set1")[1:2], "grey25"), guide = "none") +
    scale_size_continuous(guide = "none") +
    scale_x_discrete(labels = 1:18, position = "top") +
    scale_y_discrete(labels = 18:1) +
    xlab("Donor species' ID") + ylab("Recipient species' ID") + labs(tag = expression(bold(A)))

gp_S4b <-  dda_op %>% ggplot(aes(x = c_mean, y = dd_est)) +
    geom_hline(yintercept = 0.0, size = 0.25, linetype = 2) +
    geom_vline(xintercept = 0.0, size = 0.25, linetype = 2) +
    geom_errorbar(aes(xmin = c_mean - c_sd, xmax = c_mean + c_sd)) +
    geom_errorbar(aes(ymin = dd_est - dd_sd, ymax = dd_est + dd_sd)) +
    geom_point() + xlim(-0.6, 1.4) + ylim(-1.0, 0.8) +
    xlab(expression(paste("Mean per capita effect ", bar(italic(A))[italic(ij)]))) +
    ylab(expression(paste("Density dependence of ", tilde(italic(A[ij]))))) +
    labs(tag = expression(bold(B)))

tmp <- dda_op %>% mutate(sign = if_else(c_mean < 0.0, "Harm", "Bene")) %>%
    filter(!is.na(sign)) %>%
    ggplot(aes(x = sign, y = dd_est)) +
    geom_hline(yintercept = 0.0, size = 0.25, linetype = 2) +
    geom_boxplot(width = 0.5) +
    ylim(-0.6, 0.6) +
    theme_st() + theme(axis.title = element_blank())
gp_S4b <- gp_S4b + annotation_custom(ggplotGrob(tmp), xmin = 0.6, xmax = 1.5, ymin = 0.1, ymax = 0.9)

xmap_pos <- dda_op %>% mutate(xmap = str_c(rec, "_", don)) %>% filter(c_mean > 0) %>% pull(xmap)
xmap_neg <- dda_op %>% mutate(xmap = str_c(rec, "_", don)) %>% filter(c_mean < 0) %>% pull(xmap)

gp_S4c <- sim_op %>% mutate(xmap = str_c(rec, "_", don)) %>% filter(sensitivity <= 1e+3) %>%
    mutate(sign = case_when(xmap %in% xmap_pos ~ "Beneficial", xmap %in% xmap_neg ~ "Harmful", don == "all" ~ "All", TRUE ~ "Zero")) %>%
    filter(sign != "Zero") %>%
    mutate(sign = factor(sign, levels = c("Harmful", "Beneficial", "All"))) %>%
    ggplot(aes(x = 1 - rho, y = sensitivity)) +
    facet_wrap(~sign) +
    geom_hline(yintercept = 1.0, size = 0.25, linetype = 2) +
    geom_vline(xintercept = 0.0, size = 0.25, linetype = 2) +
    geom_smooth(aes(color = sign, fill = sign), method = "glm", method.args = list(family = Gamma(link = log)), size = 0.5) +
    scale_color_brewer(palette = "Set1", guide = "none") +
    scale_fill_brewer(palette = "Set1", guide = "none") +
    xlab(expression(paste("Shuffle intensity (", 1 - rho, ")"))) +
    ylab(expression(paste("Dynamical sensitivity ", frac(CV, "CV'")))) +
    labs(tag = expression((italic(c))))
 
fig_S4 <- ({(gp_S4a | gp_S4b) / gp_S4c} + plot_layout(height = c(4, 3))) & theme_st()
ggsave("fig/fig_S4.eps", fig_S4, device = cairo_ps, fallback_resolution = 1200, family = "Helvetica", width = 16, height = 15, units = "cm")
