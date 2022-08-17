####
#### analysis_dlm.R: R code for Figure 1 in "Flexibility of species interaction and ecological stability"
#### Demonstration of the concept: the analysis of density-dependent interaction
#### Time-series data are generated with a stage-structured 2-prey-1-predator model
####
#### ver 0.0.1: Initially written 20210614 by K.Kawatsu
#### ver 0.0.2: Updated on 20210623
#### ver 0.1.0: Updated on 20210701
#### ver 0.1.1: Updated on 20210716
#### ver 0.1.2: Updated on 20210718
#### ver 0.1.3: Updated on 20210722
#### ver 1.0.0: Updated on 20211117
#### ver 1.0.1: Updated on 20211122
#### ver 1.0.2: Updated on 20211215
#### ver 1.1.0: Updated on 20220808

### load source code
source("R/functions.R")

## demographic parameters
r <- c(3.65, 3.72, 3.58)
rc <- c(r[1], r[2], 1.00)
rp <- c(1.00, 1.00, r[3])

C <- matrix(c(-r[1], -0.10, 0.00,
              -0.08, -r[2], 0.00,
               0.00,  0.00, 0.00), nrow = 3, byrow = TRUE)
P <- matrix(c(0.00, 0.00, -0.20,
              0.00, 0.00, -0.10,
              0.10, 0.05, -r[3]), nrow = 3, byrow = TRUE)

beta <- 1.0

## generate time-series data
df <- gen_pp_map(rep(0.01, 3), len = 250, rc = rc, rp = rp, C = C, P = P, beta = beta)
df_ <- df %>% filter(time > 50)

## estimate per capita interaction strengths
emat <- df_ %>% gen_emat(cols = 2:4, lags = rep(0, 3))
emat[, 2:3] <- emat[, 1] * emat[, 2:3]
system.time(op1 <- smap_net(emat, range = seq(16, 18, 0.1), s = "lambda.min", alpha = 0.05))

## demonstration of interaction-effect shuffle
coef <- op1$coef[[1]]
sur1 <- coef %>% mutate(C2 = shuffle_coef(C2, lag(df_$x1), rho = 0.4, sigma = 0.01, seed = 123))
sur2 <- coef %>% mutate(C3 = shuffle_coef(C3, lag(df_$x1), rho = 0.4, sigma = 0.01, seed = 123))

## tibble of shuffled coefficients
demo_1 <- coef %>% mutate(x1 = lag(df_$x1)) %>% select(x1, C2, C3) %>% rename(C12 = C2, P13 = C3) %>%
    pivot_longer(-x1, names_to = "Coef", values_to = "estimate") %>% group_by(Coef) %>%
    group_modify(~mutate(., surrogate = shuffle_coef(estimate, x1, rho = 0.4, sigma = 0.01, seed = 123)))

## tibble of simulation demonstration
demo_2 <- df_ %>%
    mutate(original = get_pred(as.matrix(coef), cbind(emat, 1)),
           C12 = get_pred(as.matrix(sur1), cbind(emat, 1)),
           P13 = get_pred(as.matrix(sur2), cbind(emat, 1))) %>%
    select(cols = -starts_with("x")) %>%
    pivot_longer(-c(time, original), names_to = "shuffle", values_to = "prediction") %>%
    mutate(diff = case_when(
                         prediction - original < -0.01 ~ "less",
                         prediction - original > 0.01 ~ "more",
                         TRUE ~ "in"))

## tibble of the sensitivity analysis of density-dependence
set.seed(123)
cv_orig <- get_cv(get_pred(as.matrix(coef), cbind(emat, 1)))

demo_3 <-foreach(rho = seq(0.0, 1.0, 0.1), .combine = rbind) %:% foreach(iter = 1:100, .combine = rbind) %do% {
    sur1 <- coef %>% mutate(C2 = shuffle_coef(C2, lag(df_$x1), rho = rho)) %>% as.matrix()
    sur2 <- coef %>% mutate(C3 = shuffle_coef(C3, lag(df_$x1), rho = rho)) %>% as.matrix()

    tibble(rho = rho, C12 = get_cv(get_pred(sur1, cbind(emat, 1))), P13 = get_cv(get_pred(sur2, cbind(emat, 1))))
} %>% pivot_longer(-rho, names_to = "shuffle", values_to = "CV")

## make figure 1 with 4 panels
labC <- expression(italic(C)[12])
labP <- expression(italic(P)[13])
mcol <- RColorBrewer::brewer.pal(3, "Set2")[c(1, 2)]

gp_01a <- df %>% filter(between(time, 100, 200)) %>%
    pivot_longer(-time, names_to = "species", values_to = "density") %>%
    ggplot(aes(x = time)) +
    geom_line(aes(y = density, color = species), size = 0.25) +
    geom_line(data = filter(df, between(time, 100, 200)), aes(y = x1), size = 0.5) +
    scale_color_manual(values = c("black", mcol)) +
    scale_y_continuous(breaks = seq(0, 1, 0.5), limits = c(0, 1)) +
    xlab("Time") +
    ylab("Density") +
    labs(tag = expression((italic(a))))

gp_01b <- demo_1 %>% mutate(Coef = factor(Coef, levels = c("P13", "C12"))) %>%
    ggplot(aes(x = x1, color = Coef, fill = Coef)) +
    geom_smooth(aes(y = estimate), method = "loess", span = 1, size = 0.25) +
    geom_point(aes(y = surrogate), size = 0.5) +
    scale_color_manual(values = rev(mcol), labels = c(labP, labC)) +
    scale_fill_manual(values = rev(mcol), labels = c(labP, labC)) +
    xlab(expression(paste("Density (", italic(x)[1], ")"))) +
    ylab(expression(paste("Per capita effect  ",
                          frac(
                               paste(partialdiff, italic(x)[1](italic(t) + 1)),
                               paste(partialdiff, italic(x)[paste(1, italic(j))](italic(t))))))) +
    labs(tag = expression((italic(b))))

gp_01c <- demo_2 %>% filter(between(time, 150, 200)) %>%
    ggplot(aes(x = time)) + geom_line(aes(y = original), size = 0.25) +
    geom_point(aes(y = prediction, fill = shuffle, shape = diff), size = 0.75) +
    scale_fill_manual(values = mcol, guide = "none") +
    scale_shape_manual(values = c(21, 25, 24), guide = "none") +
    scale_y_continuous(breaks = seq(0, 1, 0.5), limits = c(0, 1)) +
    xlab("Time") +
    ylab(expression(paste("Prediction (", italic(x)[1], ")"))) +
    labs(tag = expression((italic(c))))

gp_01d <- demo_3 %>% mutate(shuffle = factor(shuffle, levels = c("P13", "C12"))) %>%
    ggplot(aes(x = 1 - rho, y = CV / cv_orig, color = shuffle, fill = shuffle)) +
    geom_smooth(method = "loess", span = 1, size = 0.25) +
    geom_point(position = position_jitter(width = 0.02, height = 0.0), alpha = 0.5, shape = 16, size = 0.5) +
    scale_color_manual(values = rev(mcol), labels = c(labP, labC)) +
    scale_fill_manual(values = rev(mcol), labels = c(labP, labC)) +
    scale_x_continuous(breaks = seq(0, 1, 0.5), limits = c(0, 1)) +
    scale_y_continuous(breaks = seq(0.5, 2, 0.5), limits = c(0.5, 2)) +
    xlab(expression(paste("Shuffle intensity (", 1 - rho, ")"))) +
    ylab(expression(paste("Dynamical sensitivity (", frac(CV, "CV'"), ")"))) +
    labs(tag = expression((italic(d))))

fig_01 <- (((gp_01a + plot_spacer() + plot_layout(ncol = 1, height = c(2, 1))) | gp_01b) / (gp_01c | gp_01d)) & theme_st(c(0, 0), c(0, 0))

## Prepare 'fig' directory under your current workspace
## then save figure 1 in fig
ggsave("fig/fig_01.eps", fig_01, device = cairo_ps, fallback_resolution = 600, family = "Times", width = 15, height = 12, units = "cm")

