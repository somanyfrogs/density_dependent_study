####
#### analysis_exp.R: Time-series analysis of 3-species insect experiments (Ishii & Shimada 2012, PNAS)
#### for "Temporal flexibility of species interactions and ecological stability"
####
#### ver 0.0.1: Initially written on 20190709 by K.Kawatsu
#### ver 0.0.2: Updated on 20190716
#### ver 0.0.3: Updated on 20190718
#### ver 0.1.0: Updated on 20210409
#### ver 0.1.1: Updated on 20210413
#### ver 0.2.0: Updated on 20210628
#### ver 0.2.1: Updated on 20210701
#### ver 0.3.0: Updated on 20211118
#### ver 0.4.0: Updated on 20211215
#### ver 0.4.1: Updated on 20220808
#### ver 0.4.2: Updated on 20220906

### load source code
source("R/functions.R")

### load 3-sp insect community experiments (Ishii & Shimada 2012, PNAS)
### prepare 'data' directory containing callosobruchus_3.csv under your workspace
df <- read_csv("data/callosobruchus_3.csv")

### perform density-dependence analysis of 3-species experiments
treat <- c(0.2, 0.5, 0.8)
cols <- c("Cm", "Cc", "Ac")
Tps <- c(4, 4, 2)
lgs <- c(0, 0, 3)
seed <- 123

system.time(dda_res <- foreach(br = treat, .combine = rbind) %do% {
    df_ <- df %>% filter(BR == br & !is.na(Cm) & !is.na(Cc)) %>%
        mutate(across(.cols = all_of(cols), .fns = normalize))
    lib <- df_ %>% find_knot(time_col = "Weeks")
    Es <- df_ %>% find_best_dim(cols = cols, lib = lib, range = 1:10, Tps = Tps) %>%
        select(var, E) %>% pivot_wider(names_from = var, values_from = E) %>% select(!!!syms(cols))

    foreach(i = 1:length(cols), .combine = rbind) %do% {
        ## perform Regularized S-map analysis
        df_ <- df_ %>% mutate(sv1 = !!sym(cols[i]) * lead(!!sym(cols[-i][1]), lgs[-i][1]), sv2 = !!sym(cols[i]) * lead(!!sym(cols[-i][2]), lgs[-i][2]))
        emat <- df_ %>% gen_emat(cols = c(rep(cols[i], pull(Es, cols[i]) - 2), "sv1", "sv2"), lags = c(0:-(pull(Es, cols[i]) - 3), 0, 0), knot = lib)
        op <- emat %>% smap_net(lib = lib, Tp = Tps[i], range = seq(0, 10, 1), seed = seed, alpha = 0.05)
        coef <- op %>% pull(coef) %>% .[[1]]

        ## save smap coefficients for each condition
        coef %>% write_csv(str_c("data/coef_", br, "_", cols[i], "_", cols[-i][1], "_", cols[-i][2], ".csv"))

        ## summarize estimated per-capita interaction effects of donor 1 and donor 2 on recipient species
        df_ %>% select(BR, run, Weeks, all_of(cols)) %>%
            mutate(!!cols[-i][1] := pull(coef, ncol(emat) - 1), !!cols[-i][2] := pull(coef, ncol(emat)),
                   rec_sp = cols[i], rho = op$rho, mae = op$mae, rmse = op$rmse) %>%
        pivot_longer(all_of(cols[-i]), names_to = "don_sp", values_to = "coef") %>%
        rename(rec_dens = !!sym(cols[i])) %>%
        select(BR, run, Weeks, rec_sp, don_sp, rho, mae, rmse, rec_dens, coef)
    }
})

### save result in 'data' directory
dda_res %>% write_csv("data/dda_res.csv")

## GLM analysis
dda_res <- read_csv("data/dda_res.csv")
dda_res %>% filter(rec_sp != "Ac") %>%
    ggplot(aes(x = factor(BR), y = coef)) + facet_grid(rec_sp ~ don_sp) + geom_hline(yintercept = 0, linetype = 2) + geom_violin()

dda_res %>% filter(rec_sp == "Cm" & don_sp == "Ac") %>% with(t.test(coef)) %>% tidy()
dda_res %>% filter(rec_sp == "Cm" & don_sp == "Cc") %>% with(t.test(coef)) %>% tidy()
dda_res %>% filter(rec_sp == "Cc" & don_sp == "Ac") %>% with(t.test(coef)) %>% tidy()
dda_res %>% filter(rec_sp == "Cc" & don_sp == "Cm") %>% with(t.test(coef)) %>% tidy()

a1 <- dda_res %>% filter(rec_sp == "Cm" & don_sp == "Ac") %>% with(lmer(coef ~ BR * rec_dens + (1|run)))
a1 %>% tidy()
wald.test(b = fixef(a1), Sigma = vcov(a1), Term = 4)

a2 <- dda_res %>% filter(rec_sp == "Cc" & don_sp == "Ac") %>% with(lmer(coef ~ BR * rec_dens + (1|run)))
a2 %>% tidy()
wald.test(b = fixef(a2), Sigma = vcov(a2), Term = 4)

a3 <- dda_res %>% filter(rec_sp == "Cc" & don_sp == "Cm") %>% with(lmer(coef ~ BR * rec_dens + (1|run)))
a3 %>% tidy()
wald.test(b = fixef(a3), Sigma = vcov(a3), Term = 4)

a4 <- dda_res %>% filter(rec_sp == "Cm" & don_sp == "Cc") %>% with(lmer(coef ~ BR * rec_dens + (1|run)))
a4 %>% tidy()
wald.test(b = fixef(a4), Sigma = vcov(a4), Term = 4)

### Sensitivity analysis of density-dependence alteration for the target species
treat <- c(0.2, 0.5, 0.8)
cols <- c("Cm", "Cc", "Ac")
Tps <- c(4, 4, 2)
lgs <- c(0, 0, 3)
seed <- 123
rhos <- seq(0.09, 1, 0.1)
iter <- 100

set.seed(seed)
seeds <- sample(32768, iter * length(rhos), replace = TRUE)

system.time(dds_res <- foreach(br = treat, .combine = rbind) %do% {
    df_ <- df %>% filter(BR == br & !is.na(Cm) & !is.na(Cc)) %>%
        mutate(across(.cols = all_of(cols), .fns = normalize))
    lib <- df_ %>% find_knot(time_col = "Weeks")
    Es <- df_ %>% find_best_dim(cols = cols, lib = lib, range = 1:10, Tps = Tps) %>%
        select(var, E) %>% pivot_wider(names_from = var, values_from = E) %>% select(!!!syms(cols))

    foreach(i = 1:length(cols), .combine = rbind) %do% {
        emat <- df_ %>%
            mutate(sv1 = !!sym(cols[i]) * lead(!!sym(cols[-i][1]), lgs[-i][1]), sv2 = !!sym(cols[i]) * lead(!!sym(cols[-i][2]), lgs[-i][2])) %>%
            gen_emat(cols = c(rep(cols[i], pull(Es, cols[i]) - 2), "sv1", "sv2"), lags = c(0:-(pull(Es, cols[i]) - 3), 0, 0)) %>% cbind(1)

        ### load S-map result
        coef <- read_csv(str_c("data/coef_", br, "_", cols[i], "_", cols[-i][1], "_", cols[-i][2], ".csv"))
        cnames <- names(coef)[(ncol(emat) - 2):(ncol(emat) - 1)]
        cv_o <- coef %>% as.matrix() %>% get_pred(emat) %>% get_cv()

        ## perform shuffle simulation
        foreach(rho = rhos, .combine = rbind) %:%
        foreach(j = 1:iter, .combine = rbind) %do% {
            s <- seeds[iter * (which(rhos == rho) - 1) + j]
            surr <- coef %>% mutate(across(.cols = all_of(cnames), .fns = ~{shuffle_coef(.x, dens = emat[, 1], rho = rho, seed = s)}))
            cv_1 <- coef %>% mutate(!!cnames[1] := pull(surr, cnames[1])) %>% as.matrix() %>% get_pred(emat) %>% get_cv()
            cv_2 <- coef %>% mutate(!!cnames[2] := pull(surr, cnames[2])) %>% as.matrix() %>% get_pred(emat) %>% get_cv()
            cv_b <- coef %>% mutate(!!cnames[1] := pull(surr, cnames[2]), !!cnames[2] := pull(surr, cnames[2])) %>% as.matrix() %>% get_pred(emat) %>% get_cv()

            tibble(BR = br, rho = rho, rec_sp = cols[i], !!cols[-i][1] := cv_1 / cv_o, !!cols[-i][2] := cv_2 / cv_o, both = cv_b / cv_o) %>%
                pivot_longer(cols = c(-BR, -rho, -rec_sp), names_to = "don_sp", values_to = "sensitivity")
        }
    }
})

dds_res %>% write_csv("data/dds_res.csv")

## GLM analysis
dds_res <- read_csv("data/dds_res.csv") %>% filter(rec_sp != "Ac") %>% 
    mutate(don_sp = factor(case_when(don_sp == "both" ~ "Both",
                                     don_sp == "Ac" ~ "Parasitism",
                                     TRUE ~ "Competition"), levels = c("Parasitism", "Competition", "Both")))
glm(sensitivity ~ rec_sp * rho + BR, data = filter(dds_res, don_sp == "Parasitism"), family = Gamma(link = log)) %>% tidy()
glm(sensitivity ~ rec_sp * rho + BR, data = filter(dds_res, don_sp == "Competition"), family = Gamma(link = log)) %>% tidy()
glm(sensitivity ~ rec_sp * rho + BR, data = filter(dds_res, don_sp == "Both"), family = Gamma(link = log)) %>% tidy()

## make Figure 2a
gp_02a <- read_csv("data/dda_res.csv") %>% filter(rec_sp != "Ac") %>%
    mutate(don_sp = factor(ifelse(don_sp == "Ac", "Parasitism", "Competition"), levels = c("Parasitism", "Competition"))) %>%
    ggplot(aes(x = rec_dens, y = coef, color = factor(BR))) +
    facet_grid(rec_sp ~ don_sp) +
    geom_hline(yintercept = 0.0, linetype = 2, size = 0.25) +
    geom_point(shape = 16, size = 0.25, alpha = 0.5) +
    geom_smooth(aes(fill = factor(BR)), method = "loess", span = 1, size = 0.5) +
    scale_color_brewer(palette = "Set2") + scale_fill_brewer(palette = "Set2") +
    scale_x_continuous(breaks = seq(0.0, 1.0, 0.5), limits = c(0.0, 1.0)) +
    xlab("Normalized density of recipient species") +
    ylab(expression(paste("Per capita effect  ",
                          frac(
                               paste(partialdiff, italic(x)[italic(i)](italic(t) + 1)),
                               paste(partialdiff, italic(x)[italic(ij)](italic(t))))))) +
    labs(tag = expression(bold(A)))

## make Figure 2b
gp_02b <- read_csv("data/dds_res.csv") %>%
    filter(rec_sp != "Ac") %>%
    mutate(don_sp = factor(case_when(
                                     don_sp == "both" ~ "Both",
                                     don_sp == "Ac" ~ "Parasitism",
                                     TRUE ~ "Competition"
                                     ), levels = c("Parasitism", "Competition", "Both"))) %>%
    ggplot(aes(x = 1 - rho, y = sensitivity, color = factor(BR))) +
    facet_grid(rec_sp ~ don_sp) +
    geom_hline(yintercept = 1.0, linetype = 2, size = 0.25) +
    geom_point(position = position_jitter(width = 0.02, height = 0.0), shape = 16, size = 0.25, alpha = 0.5) +
    geom_smooth(aes(fill = factor(BR)), method = "loess", span = 1, size = 0.5) +
    scale_color_brewer(palette = "Set2", guide = FALSE) + scale_fill_brewer(palette = "Set2", guide = FALSE) +
    scale_x_continuous(breaks = seq(0.0, 1.0, 0.5), limits = c(0.0, 1.0)) +
    scale_y_continuous(breaks = seq(1.0, 4.0, 1.0)) +
    xlab(expression(paste("Shuffle intensity (", 1 - rho, ")"))) +
    ylab(expression(paste("Dynamical sensitivity ", frac(CV, "CV'")))) +
    labs(tag = expression(bold(B)))

fig_02 <- (gp_02a + gp_02b + plot_layout(nrow = 1, width = c(2, 3))) & theme_st(just = c(0, 0), pos = c(0, 0))

## Prepare 'fig' directory under your current workspace
## then save figure 2 in fig
ggsave("fig/fig_02.eps", fig_02, device = cairo_ps, fallback_resolution = 1200, family = "Helvetica", width = 16, height = 8, units = "cm")

