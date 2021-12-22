#' Script for comparison of Stan implementations of random walk models
#' 

# on Markov:
#   ----   Load libraries and set workspace   ----
set_workspace <- function(markov=TRUE){
  if(markov){
    .libPaths("~/Documents/R_libraries")
    setwd("~/Documents/GitHub/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Random walk implementations/Rw1_fh_testing")
  } else {
    setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Random walk implementations/Rw1_fh_testing")
  }
}

set_workspace(markov=F)

library("rstan")
library("inlabru")
library("tidyverse")
library("ggplot2")
library("patchwork")
library("INLA")

# generate test data:
n=100
z=seq(0,6,length.out=n)
y=sin(z)+rnorm(n,mean=0,sd=0.5)
data=data.frame(y=y,z=z)

# quickly estimate in inlabru:
components = ~ -1 +
  eta(z, model = "rw1", constr = TRUE, hyper = list(prec = list(intitial = log(7500), fixed = TRUE)))

formula = y ~ eta
likelihood = like(formula = formula, family = "gaussian", data = data)

c.family <- list(hyper = list(prec = list(initial = log(4.5), fixed = TRUE)))

res.inlabru <- bru(components = components,
                   likelihood = likelihood, 
                   options = list(verbose = F, 
                                  bru_verbose = 1,
                                  control.family  = c.family))


ggplot(data.frame(eta = res.inlabru$summary.random$eta$mean, z = z, sin.z = sin(z))) + geom_point(aes(x = z, y = eta)) + geom_point(aes(x = z, y = sin.z))

###   ----   Run stan analyses   ---- 
input_stan <- list(y = data$y)

# run first stan program: rw1_stepwise.stan:

fit_stepwise <- stan(
  file = "rw1_fh_stepwise_test.stan",
  data = input_stan,
  chains = 4,
  iter = 4000,
  warmup = 400,
  refresh = 1000,
  seed = 123
)

list_of_draws_stepwise <- rstan::extract(fit_stepwise)

# save or load data
save(list_of_draws_stepwise, file = 'draws_stepwise.RData')
load('draws_stepwise.RData')

eta_df_stepwise <- data.frame(list_of_draws_stepwise$eta)
summary_stepwise <- as.data.frame(summary(fit_stepwise)$summary) %>% mutate("Implementation" = "Stepwise")

trace_stepwise <- plot(fit_stepwise, plotfun = "trace", pars = c("tau_y", "tau_eta", "eta[1]", "eta[50]", "eta[90]", "eta[100]"))
ggsave("trace_stepwise.pdf", plot = trace_stepwise, dpi="retina", device = "pdf")

#  run second stan program: rw1_by_differences.stan

fit_diff_1 <- stan(
  file = "rw1_fh_by_differences_1_test.stan",
  data = input_stan,
  chains = 4,
  iter = 4000,
  warmup = 400,
  refresh = 1000,
  seed = 123
)

trace_diff_1 <- plot(fit_diff_1, plotfun = "trace", pars = c("tau_y", "tau_eta", "eta[1]", "eta[50]", "eta[90]", "eta[100]"))
ggsave("trace_diff_1.pdf", plot = trace_diff_1, dpi="retina", device = "pdf")

list_of_draws_diff_1 <- rstan::extract(fit_diff_1)

# save or load draws:
save(list_of_draws_diff_1, file = 'draws_diff_1.RData')
load('draws_diff_1.RData')

eta_df_diff_1 <- data.frame(list_of_draws_diff_1$eta)
summary_diff_1 <- as.data.frame(summary(fit_diff_1)$summary) %>% mutate("Implementation" = "Diff.1")

#  run third stan program: rw1_by_differences_icar.stan

fit_diff_2 <- stan(
  file = "rw1_fh_by_differences_2_test.stan",
  data = input_stan,
  chains = 4,
  iter = 4000,
  warmup = 400,
  refresh = 1000,
  seed = 123
)

trace_diff_2 <- plot(fit_diff_2, plotfun = "trace", pars = c("tau_y", "tau_eta", "eta[1]", "eta[50]", "eta[90]", "eta[100]"))
ggsave("trace_diff_2.pdf", plot = trace_diff_2, dpi="retina", device = "pdf")

list_of_draws_diff_2 <- rstan::extract(fit_diff_2)

# save or load draws
save(list_of_draws_diff_2, file = 'draws_diff_2.RData')
load("draws_diff_2.RData")

eta_df_diff_2 <- data.frame(list_of_draws_diff_2$eta)
summary_diff_2 <- as.data.frame(summary(fit_diff_2)$summary) %>% mutate("Implementation" = "Diff.2")

#  run fourth stan program: rw1_by_differences_3.stan

fit_diff_3 <- stan(
  file = "rw1_fh_by_differences_3-test.stan",
  data = input_stan,
  chains = 4,
  iter = 4000,
  warmup = 400,
  refresh = 1000,
  seed = 123
)

trace_diff_3 <- plot(fit_diff_3, plotfun = "trace", pars = c("tau_y", "tau_eta", "eta[1]", "eta[50]", "eta[90]", "eta[100]"))
ggsave("trace_diff_3.pdf", plot = trace_diff_3, dpi="retina", device = "pdf")

list_of_draws_diff_3 <- rstan::extract(fit_diff_3)

# save or load draws
save(list_of_draws_diff_3, file = 'draws_diff_3.RData')
load("draws_diff_3.RData")

eta_df_diff_3 <- data.frame(list_of_draws_diff_3$eta)
summary_diff_3 <- as.data.frame(summary(fit_diff_3)$summary) %>% mutate("Implementation" = "Diff.3")

#  run fourth stan program: rw1_by_differences_4.stan

fit_diff_4 <- stan(
  file = "rw1_fh_by_differences_4_test.stan",
  data = input_stan,
  chains = 4,
  iter = 4000,
  warmup = 400,
  refresh = 1000,
  seed = 123
)

trace_diff_4 <- plot(fit_diff_4, plotfun = "trace", pars = c("tau_y", "tau_eta", "eta[1]", "eta[50]", "eta[90]", "eta[100]"))
ggsave("trace_diff_4.pdf", plot = trace_diff_4, dpi="retina", device = "pdf")

list_of_draws_diff_4 <- rstan::extract(fit_diff_4)

# save or load draws
save(list_of_draws_diff_4, file = 'draws_diff_4.RData')
load("draws_diff_4.RData")

eta_df_diff_4 <- data.frame(list_of_draws_diff_4$eta)
summary_diff_4 <- as.data.frame(summary(fit_diff_4)$summary) %>% mutate("Implementation" = "Diff.4")

###   ----   Compare results of different rw2 implementations   ----:


palette <- c('#70A4D4', '#ECC64B', '#93AD80', '#da9124', '#696B8D',
             '#3290c1', '#5d8060', '#D7B36A', '#826133', '#A85150')

eta_df_stepwise_long <- eta_df_stepwise %>% pivot_longer(cols = everything(), names_to = "Param", values_to = "Stepwise")
eta_df_diff_1_long <- eta_df_diff_1 %>% pivot_longer(cols = everything(), names_to = "Param", values_to = "Diff.1")
eta_df_diff_2_long <- eta_df_diff_2 %>% pivot_longer(cols = everything(), names_to = "Param", values_to = "Diff.2")
eta_df_diff_3_long <- eta_df_diff_3 %>% pivot_longer(cols = everything(), names_to = "Param", values_to = "Diff.3")
eta_df_diff_4_long <- eta_df_diff_4 %>% pivot_longer(cols = everything(), names_to = "Param", values_to = "Diff.4")

eta_all <- cbind(eta_df_stepwise_long, Diff.1 = eta_df_diff_1_long$Diff.1, Diff.2 = eta_df_diff_2_long$Diff.2) %>% 
  mutate(idx = parse_number(Param))

eta_all_34 <- cbind(eta_df_stepwise_long, Diff.1 = eta_df_diff_1_long$Diff.1,
                     Diff.2 = eta_df_diff_2_long$Diff.2,
                     Diff.3 = eta_df_diff_3_long$Diff.3,
                     Diff.4 = eta_df_diff_4_long$Diff.4) %>% 
  mutate(idx = parse_number(Param))

# every 10th value:

eta_10 <- eta_all %>% filter(idx %% 10 == 0)

eta.plot <- ggplot(eta_10) +
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.1, y = after_stat(density), color = "Difference 1", fill = "Difference 1"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  facet_wrap(~Param)

ggsave("eta_histogram_10.pdf", plot = eta.plot, dpi="retina", device = "pdf")

# inlcuding diff 3 and 4
eta_10_34 <- eta_all_34 %>% filter(idx %% 10 == 0)

eta.plot.34 <- ggplot(eta_10_34) +
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.1, y = after_stat(density), color = "Difference 1", fill = "Difference 1"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.3, y = after_stat(density), color = "Difference 3", fill = "Difference 3"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.4, y = after_stat(density), color = "Difference 4", fill = "Difference 4"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  facet_wrap(~Param)

ggsave("eta_histogram_10_34.pdf", plot = eta.plot.34, dpi="retina", device = "pdf")

# as density:
eta.plot.34.dens <- ggplot(eta_10_34) +
  geom_density(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), alpha = 0.2) +
  geom_density(aes(x = Diff.1, y = after_stat(density), color = "Difference 1", fill = "Difference 1"), alpha = 0.2) +
  geom_density(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), alpha = 0.2) +
  geom_density(aes(x = Diff.3, y = after_stat(density), color = "Difference 3", fill = "Difference 3"), alpha = 0.2) +
  geom_density(aes(x = Diff.4, y = after_stat(density), color = "Difference 4", fill = "Difference 4"), alpha = 0.2) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  facet_wrap(~Param)

ggsave("eta_dens_10_34.pdf", plot = eta.plot.34.dens, dpi="retina", device = "pdf", height = 5, width = 8)

# only diff.1 and stepwise:

eta.plo.wo.2 <- ggplot(eta_10) +
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.1, y = after_stat(density), color = "Difference 1", fill = "Difference 1"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  facet_wrap(~Param)

ggsave("eta_histogram_wo_2.pdf", plot = eta.plo.wo.2, dpi="retina", device = "pdf")

# only diff 2,3,4

eta.plot.234 <- ggplot(eta_10_34) +
  geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.3, y = after_stat(density), color = "Difference 3", fill = "Difference 3"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.4, y = after_stat(density), color = "Difference 4", fill = "Difference 4"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  facet_wrap(~Param)

ggsave("eta_histogram_10_234.pdf", plot = eta.plot.234, dpi="retina", device = "pdf", height = 5, width = 8)

# as density plot:
eta.plot.234.dens <- ggplot(eta_10_34) +
  geom_density(aes(x = Diff.2, color = "Difference 2", fill = "Difference 2"), alpha = 0.2) +
  geom_density(aes(x = Diff.3, color = "Difference 3", fill = "Difference 3"), alpha = 0.2) +
  geom_density(aes(x = Diff.4, color = "Difference 4", fill = "Difference 4"), alpha = 0.2) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  facet_wrap(~Param)

ggsave("eta_dens_10_234.pdf", plot = eta.plot.234.dens, dpi="retina", device = "pdf", height = 5, width = 8)

# all values of the predictor

eta.plot.full <- ggplot(eta_all) +
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.1, y = after_stat(density), color = "Difference 1", fill = "Difference 1"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  facet_wrap(~Param)

ggsave("eta_histogram_full.pdf", plot = eta.plot.full, dpi="retina", device = "pdf")

#   ---   Plot means from summary together with "true" value   ----

eta_summary_stewpwise <- summary_stepwise%>%
  rownames_to_column(var="idx") %>%
  filter(grepl("eta", idx)) %>% filter(!grepl("tau_eta", idx)) %>%
  mutate(idx = parse_number(idx))

eta_summary_diff_1 <- summary_diff_1%>%
  rownames_to_column(var="idx") %>%
  filter(grepl("eta", idx)) %>% filter(!grepl("tau_eta", idx)) %>% filter(!grepl("diff", idx)) %>%
  mutate(idx = parse_number(idx))

eta_summary_diff_2 <- summary_diff_2%>%
  rownames_to_column(var="idx") %>%
  filter(grepl("eta", idx)) %>% filter(!grepl("tau_eta", idx)) %>% filter(!grepl("scaled", idx)) %>%
  mutate(idx = parse_number(idx))

summary_eta <- data.frame(idx = eta_summary_stewpwise$idx,
                          Stepwise = eta_summary_stewpwise$mean, 
                          Diff.1 = eta_summary_diff_1$mean,
                          Diff.2 = eta_summary_diff_2$mean) %>%
  mutate(True = sin(z)[idx])

p.summary.eta <- ggplot(summary_eta, aes(x = idx)) + 
  geom_point(aes(y = Stepwise, color = "Stepwise"), shape = 1) + 
  geom_point(aes(y = Diff.1, color = "Difference 1"), shape = 1, alpha = 0.5) + 
  geom_point(aes(y = Diff.2, color = "Difference 2"), shape = 1, alpha = 0.5) + 
  geom_point(aes(y = True, color = "True"), alpha = 0.5) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs("Means of eta", x = "z", y = "")

ggsave("summary_eta.pdf", plot = p.summary.eta, dpi="retina", device = "pdf", height = 5, width = 8)

p.summary.eta.ib <- ggplot(summary_eta, aes(x = idx)) + 
  geom_point(aes(y = Stepwise, color = "Stepwise"), shape = 1) + 
  geom_point(aes(y = Diff.1, color = "Difference 1"), shape = 1, alpha = 0.5) + 
  geom_point(aes(y = Diff.2, color = "Difference 2"), shape = 1, alpha = 0.5) + 
  geom_point(data = data.frame(res.inlabru$summary.random$eta) %>% mutate(idx = summary_eta$idx), aes(x = idx, y = mean, color = "Inlabru"), alpha = 0.5, shape = 1) + 
  geom_point(aes(y = True, color = "True"), alpha = 0.5) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs("Means of eta", x = "z", y = "")

ggsave("summary_eta_ib.pdf", plot = p.summary.eta.ib, dpi="retina", device = "pdf", height = 5, width = 8)

#   ----    Plot comparison of Stan and inlabru   ----

eta_ib <- data.frame(res.inlabru$marginals.random$eta)

p.eta.1 <- ggplot(eta_all %>% filter(idx == 1)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.1), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.1, y = after_stat(density), color = "Difference 1", fill = "Difference 1"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 1", x = "", y = "")

p.eta.20 <- ggplot(eta_all %>% filter(idx == 20)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.20), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.1, y = after_stat(density), color = "Difference 1", fill = "Difference 1"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 20", x = "", y = "")

p.eta.40 <- ggplot(eta_all %>% filter(idx == 40)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.40), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.1, y = after_stat(density), color = "Difference 1", fill = "Difference 1"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 40", x = "", y = "")

p.eta.60 <- ggplot(eta_all %>% filter(idx == 60)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.60), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.1, y = after_stat(density), color = "Difference 1", fill = "Difference 1"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 60", x = "", y = "")

p.eta.80 <- ggplot(eta_all %>% filter(idx == 80)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.80), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.1, y = after_stat(density), color = "Difference 1", fill = "Difference 1"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 80", x = "", y = "")

p.eta.100 <- ggplot(eta_all %>% filter(idx == 100)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.100), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.1, y = after_stat(density), color = "Difference 1", fill = "Difference 1"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 100", x = "", y = "")

p.eta.ib <- (p.eta.1 | p.eta.20 | p.eta.40)/(p.eta.60 | p.eta.80 | p.eta.100) + 
  plot_layout(guides = "collect")

ggsave("eta_inlabru.pdf", plot = p.eta.ib, device = "pdf", height = 5, width = 8, dpi = "retina")

#   ----    Plot comparison of Diff.2 and inlabru   ----

eta_ib <- data.frame(res.inlabru$marginals.random$eta)

p.eta.1.2 <- ggplot(eta_all %>% filter(idx == 1)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.1), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 1", x = "", y = "")

p.eta.20.2 <- ggplot(eta_all %>% filter(idx == 20)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.20), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 20", x = "", y = "")

p.eta.40.2 <- ggplot(eta_all %>% filter(idx == 40)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.40), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 40", x = "", y = "")

p.eta.60.2 <- ggplot(eta_all %>% filter(idx == 60)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.60), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 60", x = "", y = "")

p.eta.80.2 <- ggplot(eta_all %>% filter(idx == 80)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.80), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 80", x = "", y = "")

p.eta.100.2 <- ggplot(eta_all %>% filter(idx == 100)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.100), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 100", x = "", y = "")

p.eta.ib.2 <- (p.eta.1.2 | p.eta.20.2 | p.eta.40.2)/(p.eta.60.2 | p.eta.80.2 | p.eta.100.2) + 
  plot_layout(guides = "collect")

ggsave("eta_inlabru_diff_2.pdf", plot = p.eta.ib.2, device = "pdf", height = 5, width = 8, dpi = "retina")

#   ----    Plot comparison of stepwise and inlabru   ----

p.eta.1.sw <- ggplot(eta_all %>% filter(idx == 1)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.1), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 1", x = "", y = "")

p.eta.20.sw <- ggplot(eta_all %>% filter(idx == 20)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.20), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 20", x = "", y = "")

p.eta.40.sw <- ggplot(eta_all %>% filter(idx == 40)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.40), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 40", x = "", y = "")

p.eta.60.sw <- ggplot(eta_all %>% filter(idx == 60)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.60), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 60", x = "", y = "")

p.eta.80.sw <- ggplot(eta_all %>% filter(idx == 80)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.80), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 80", x = "", y = "")

p.eta.100.sw <- ggplot(eta_all %>% filter(idx == 100)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.100), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 100", x = "", y = "")

p.eta.ib.sw <- (p.eta.1.sw | p.eta.20.sw | p.eta.40.sw)/(p.eta.60.sw | p.eta.80.sw | p.eta.100.sw) + 
  plot_layout(guides = "collect")

ggsave("eta_inlabru_stepwise.pdf", plot = p.eta.ib.sw, device = "pdf", height = 5, width = 8, dpi = "retina")

#   ----    Plot comparison of Difference 1 and inlabru   ----

p.eta.1.1 <- ggplot(eta_all %>% filter(idx == 1)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.1), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.1, y = after_stat(density), color = "Difference 1", fill = "Difference 1"), bins = 500, alpha = 0.5, size = 0) +
  #geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 1", x = "", y = "")

p.eta.20.1 <- ggplot(eta_all %>% filter(idx == 20)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.20), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  #geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.1, y = after_stat(density), color = "Difference 1", fill = "Difference 1"), bins = 500, alpha = 0.5, size = 0) +
  #geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 20", x = "", y = "")

p.eta.40.1 <- ggplot(eta_all %>% filter(idx == 40)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.40), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  #geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.1, y = after_stat(density), color = "Difference 1", fill = "Difference 1"), bins = 500, alpha = 0.5, size = 0) +
  #geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 40", x = "", y = "")

p.eta.60.1 <- ggplot(eta_all %>% filter(idx == 60)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.60), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  #geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.1, y = after_stat(density), color = "Difference 1", fill = "Difference 1"), bins = 500, alpha = 0.5, size = 0) +
  #geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 60", x = "", y = "")

p.eta.80.1 <- ggplot(eta_all %>% filter(idx == 80)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.80), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  #geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.1, y = after_stat(density), color = "Difference 1", fill = "Difference 1"), bins = 500, alpha = 0.5, size = 0) +
  #geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 80", x = "", y = "")

p.eta.100.1 <- ggplot(eta_all %>% filter(idx == 100)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.100), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  #geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.1, y = after_stat(density), color = "Difference 1", fill = "Difference 1"), bins = 500, alpha = 0.5, size = 0) +
  #geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 100", x = "", y = "")

p.eta.ib.1 <- (p.eta.1.1 | p.eta.20.1 | p.eta.40.1)/(p.eta.60.1 | p.eta.80.1 | p.eta.100.1) + 
  plot_layout(guides = "collect")

ggsave("eta_inlabru_diff_1.pdf", plot = p.eta.ib.1, device = "pdf", height = 5, width = 8, dpi = "retina")

#   ----    Plot comparison of Diff.3 and inlabru   ----

p.eta.1.3 <- ggplot(eta_all_34 %>% filter(idx == 1)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.1), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.3, y = after_stat(density), color = "Difference 3", fill = "Difference 3"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 1", x = "", y = "")

p.eta.20.3 <- ggplot(eta_all_34 %>% filter(idx == 20)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.20), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.3, y = after_stat(density), color = "Difference 3", fill = "Difference 3"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 20", x = "", y = "")

p.eta.40.3 <- ggplot(eta_all_34 %>% filter(idx == 40)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.40), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.3, y = after_stat(density), color = "Difference 3", fill = "Difference 3"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 40", x = "", y = "")

p.eta.60.3 <- ggplot(eta_all_34 %>% filter(idx == 60)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.60), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.3, y = after_stat(density), color = "Difference 3", fill = "Difference 3"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 60", x = "", y = "")

p.eta.80.3 <- ggplot(eta_all_34 %>% filter(idx == 80)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.80), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.3, y = after_stat(density), color = "Difference 3", fill = "Difference 3"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 80", x = "", y = "")

p.eta.100.3 <- ggplot(eta_all_34 %>% filter(idx == 100)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.100), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.3, y = after_stat(density), color = "Difference 3", fill = "Difference 3"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 100", x = "", y = "")

p.eta.ib.3 <- (p.eta.1.3 | p.eta.20.3 | p.eta.40.3)/(p.eta.60.3 | p.eta.80.3 | p.eta.100.3) + 
  plot_layout(guides = "collect")

ggsave("eta_inlabru_diff_3.pdf", plot = p.eta.ib.3, device = "pdf", height = 5, width = 8, dpi = "retina")


#   ----    Plot comparison of Diff.3 and inlabru   ----

p.eta.1.4 <- ggplot(eta_all_34 %>% filter(idx == 1)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.1), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.4, y = after_stat(density), color = "Difference 4", fill = "Difference 4"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 1", x = "", y = "")

p.eta.20.4 <- ggplot(eta_all_34 %>% filter(idx == 20)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.20), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.4, y = after_stat(density), color = "Difference 4", fill = "Difference 4"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 20", x = "", y = "")

p.eta.40.4 <- ggplot(eta_all_34 %>% filter(idx == 40)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.40), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.3, y = after_stat(density), color = "Difference 4", fill = "Difference 4"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 40", x = "", y = "")

p.eta.60.4 <- ggplot(eta_all_34 %>% filter(idx == 60)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.60), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.3, y = after_stat(density), color = "Difference 4", fill = "Difference 4"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 60", x = "", y = "")

p.eta.80.4 <- ggplot(eta_all_34 %>% filter(idx == 80)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.80), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.4, y = after_stat(density), color = "Difference 4", fill = "Difference 4"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 80", x = "", y = "")

p.eta.100.4 <- ggplot(eta_all_34 %>% filter(idx == 100)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.100), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.4, y = after_stat(density), color = "Difference 4", fill = "Difference 4"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 100", x = "", y = "")

p.eta.ib.4 <- (p.eta.1.4 | p.eta.20.4 | p.eta.40.4)/(p.eta.60.4 | p.eta.80.4 | p.eta.100.4) + 
  plot_layout(guides = "collect")

ggsave("eta_inlabru_diff_4.pdf", plot = p.eta.ib.4, device = "pdf", height = 5, width = 8, dpi = "retina")





