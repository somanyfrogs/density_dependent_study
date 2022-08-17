####
#### analysis_fish.R: Time-series analysis of marine fish monitoring data (Masuda et al. 2009)
#### for "Flexibility of species interactions and ecological stability"
####
#### ver 0.0.1: Initially written on 20190718 by K.Kawatsu
#### ver 0.0.2: Updated on 20190726
#### ver 0.0.3: Updated on 20190918
#### ver 0.1.0: Updated on 20210702
#### ver 0.1.1: Updated on 20210715
#### ver 0.1.2: Updated on 20210722
#### ver 0.2.0: Updated on 20211118
#### ver 0.3.0: Updated on 20211215
#### ver 0.4.0: Updated on 20220808

### load source code
source("R/functions.R")

### load time-series dataset of maiduru fish monitoring
### prepare 'data' directory containing maiduru_obs.csv under your current workspace
fish_orig <- read_csv("data/maiduru_obs.csv")

### make monthly averaged data
fish <- fish_orig %>%
    mutate(date = as.Date(date, origin = as.Date("2002-01-01")), year = str_sub(date, 1, 4), month = str_sub(date, 6, 7)) %>%
    select(-date) %>% group_by(year, month) %>%
    group_modify(~summarize(., across(.cols = everything(), .fns = ~mean(.x, na.rm = TRUE)))) %>% ungroup()

### summarize dataset
fish_sum <- fish %>% select(-year, -month) %>%
    pivot_longer(cols = everything(), names_to = "Sp", values_to = "Count") %>% group_by(Sp) %>%
    group_modify(~mutate(., PA = if_else(Count > 0, 1, 0))) %>%
    group_modify(~summarize(., Total = sum(Count), Obs = sum(PA))) %>%
    arrange(desc(Obs)) %>% ungroup()

fish_sum %>% write_csv("data/fish_sum.csv")

### Species selection (observed more than 1 / 3 of total data length (201))
cols <- read_csv("data/fish_sum.csv") %>% filter(Obs >= nrow(fish) / 3) %>% pull(Sp)
fish <- fish %>% mutate(across(.cols = all_of(cols), .fns = normalize))
Es <- fish %>% find_best_dim(cols = cols, range = 2:12)
Es %>% write_csv("data/fish_Es.csv")

### make surrogate data for each fish species
set.seed(123)
seeds <- sample(32768, length(cols), replace = TRUE)
system.time(surr_list <- foreach(sp = cols) %do% {
                E <- Es %>% filter(var == sp) %>% pull(E)
                emat <- fish %>% gen_emat(cols = rep(sp, E), lags = 0:-(E - 1))
                emat %>% gen_surr_ts(col = 1, time_idx = 1:nrow(emat), iter = 1000, seed = seeds[which(sp == cols)])
} %>% setNames(cols))

### Causality test
cl <- makeCluster(detectCores(), type = "PSOCK")
registerDoParallel(cl)

system.time(ccm_op <-
foreach(ref = cols, .combine = rbind) %:%
foreach(tar = cols[cols != ref], .combine = rbind, .packages = c("tidyverse", "foreach", "rpkg")) %dopar% {
    E <- Es %>% filter(var == ref) %>% pull(E)

    foreach(Tp = -3:0, .combine = rbind) %do% {
        tryCatch({
            system.time(op <- fish %>% ccm_test(E = E, Tp = Tp, ref = ref, tar = tar, surr = surr_list[[ref]]))
            op_max <- op %>% filter(data == "max")
            p_min <- op %>% filter(data == "min") %>% summarize(p = sum(op_max$rmse >= rmse) / nrow(.)) %>% pull(p)
            p_sur <- op %>% filter(data == "surr") %>% summarize(p = sum(op_max$rmse >= rmse) / nrow(.)) %>% pull(p)
            op_max %>% select(-data) %>% mutate(p_min = p_min, p_sur = p_sur)
        }, error = function(e) tibble(ref = ref, tar = tar, E = E, Tp = Tp, rho = NA, mae = NA, rmse = NA, p_min = NA, p_max = NA))
    }
})

### save CCM result of 13 fish species in 'data' directory
ccm_op %>% write_csv("data/ccm_op.csv")
stopCluster(cl)

### Jacobian-estimation analysis and sensitivity analysis
### load time-series dataset of maiduru fish monitoring
fish_orig <- read_csv("data/maiduru_obs.csv")

### make monthly averaged data
fish <- fish_orig %>%
    mutate(date = as.Date(date, origin = as.Date("2002-01-01")), year = str_sub(date, 1, 4), month = str_sub(date, 6, 7)) %>%
    select(-date) %>% group_by(year, month) %>%
    group_modify(~summarize(., across(.cols = everything(), .fns = ~mean(.x, na.rm = TRUE)))) %>% ungroup()

### Species selection (observed more than 1 / 2 total data length (201))
ccm_op <- read_csv("data/ccm_op.csv")
cols <- read_csv("data/fish_sum.csv") %>% filter(Obs >= 201 * 0.5) %>% pull(Sp)
fish <- fish %>% mutate(across(.cols = all_of(cols), .fns = normalize))

seed <- 123
rhos <- seq(0.09, 1, 0.1)
iter <- 100

seeds <- sample(32768, iter * length(rhos), replace = TRUE)

system.time(smap_op_01 <- foreach(rec = cols, .combine = rbind) %do% {
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

## save density-dependence analysis of 13 fish species in 'data' directory
foreach(i = 1:nrow(smap_op_01), .combine = rbind) %do% {smap_op_01$dda[[i]]} %>% write_csv("data/dda_op_01.csv")
foreach(i = 1:nrow(smap_op_01), .combine = rbind) %do% {smap_op_01$sim[[i]]} %>% write_csv("data/sim_op_01.csv")

### make Figure 3
dda_op <- read_csv("data/dda_op_01.csv")
sim_op <- read_csv("data/sim_op_01.csv")
cols <- read_csv("data/fish_sum.csv") %>% filter(Obs >= 201 * 0.5) %>% pull(Sp)

gp_03a <- dda_op %>% mutate(rec = factor(rec, levels = rev(cols)), don = factor(don, levels = cols)) %>%
    mutate(sign = case_when(c_mean < 0 & p_value1 <= 0.05 ~ "neg", c_mean > 0 & p_value1 <= 0.05 ~ "pos", TRUE ~ "zero")) %>%
    ggplot(aes(x = don, y = rec, color = sign)) +
    geom_point(aes(size = abs(c_mean))) +
    scale_color_manual(values = c(RColorBrewer::brewer.pal(3, "Set1")[1:2], "grey25"), guide = "none") +
    scale_size_continuous(guide = "none") +
    scale_x_discrete(labels = 1:13, position = "top") +
    scale_y_discrete(labels = 13:1) +
    xlab("Donor species' ID") + ylab("Recipient species' ID") + labs(tag = expression((italic(a))))

gp_03b <- dda_op %>% ggplot(aes(x = c_mean, y = dd_est)) +
    geom_hline(yintercept = 0.0, size = 0.25, linetype = 2) +
    geom_vline(xintercept = 0.0, size = 0.25, linetype = 2) +
    geom_errorbar(aes(xmin = c_mean - c_sd, xmax = c_mean + c_sd)) +
    geom_errorbar(aes(ymin = dd_est - dd_sd, ymax = dd_est + dd_sd)) +
    geom_point() + xlim(-0.6, 1.4) + ylim(-1.0, 0.8) +
    xlab(expression(paste("Mean per capita effect ", italic(bar(c)[ij])))) +
    ylab(expression(paste("Density dependence of ", italic(hat(c[ij]))))) +
    labs(tag = expression((italic(b))))

tmp <- dda_op %>% mutate(sign = if_else(c_mean < 0.0, "Harm", "Bene")) %>%
    filter(!is.na(sign)) %>%
    ggplot(aes(x = sign, y = dd_est)) +
    geom_hline(yintercept = 0.0, size = 0.25, linetype = 2) +
    geom_boxplot(width = 0.5) +
    ylim(-1.0, 1.0) +
    theme_st() + theme(axis.title = element_blank())
gp_03b <- gp_03b + annotation_custom(ggplotGrob(tmp), xmin = 0.6, xmax = 1.5, ymin = 0.1, ymax = 0.9)

xmap_pos <- dda_op %>% mutate(xmap = str_c(rec, "_", don)) %>% filter(c_mean > 0) %>% pull(xmap)
xmap_neg <- dda_op %>% mutate(xmap = str_c(rec, "_", don)) %>% filter(c_mean < 0) %>% pull(xmap)

gp_03c <- sim_op %>% mutate(xmap = str_c(rec, "_", don)) %>%
    mutate(sign = case_when(xmap %in% xmap_pos ~ "Beneficial", xmap %in% xmap_neg ~ "Harmful", don == "all" ~ "All", TRUE ~ "Zero")) %>%
    filter(sign != "Zero") %>%
    mutate(sign = factor(sign, levels = c("Harmful", "Beneficial", "All"))) %>%
    ggplot(aes(x = 1 - rho, y = sensitivity)) +
    facet_wrap(~sign) +
    geom_hline(yintercept = 1.0, size = 0.25, linetype = 2) +
    geom_vline(xintercept = 0.0, size = 0.25, linetype = 2) +
    geom_smooth(aes(color = sign, fill = sign), size = 0.5) +
    scale_color_brewer(palette = "Set1", guide = "none") +
    scale_fill_brewer(palette = "Set1", guide = "none") +
    xlab(expression(paste("Shuffle intensity (", 1 - rho, ")"))) +
    ylab(expression(paste("Dynamical sensitivity ", frac(CV, "CV'")))) +
    labs(tag = expression((italic(c))))

fig_03 <- ({(gp_03a | gp_03b) / gp_03c} + plot_layout(height = c(4, 3))) & theme_st()

## Prepare 'fig' directory under your current workspace
## then save figure 3 in fig
ggsave("fig/fig_03.eps", fig_03, device = cairo_ps, fallback_resolution = 600, family = "Times", width = 16, height = 15, units = "cm")

