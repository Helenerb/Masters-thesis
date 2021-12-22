#' Script for comparison of Stan implementations of random walk models
#' 

library("rstan")
library("inlabru")
library("tidyverse")
library("ggplot2")
library("patchwork")

# set workspace to current location:
setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Random walk implementations")

# generate test data:
n=100
z=seq(0,6,length.out=n)
y=sin(z)+rnorm(n,mean=0,sd=0.5)
data=data.frame(y=y,z=z)

# quickly estimate in inlabru:
# components = ~ -1 + 
#   eta(z, model = "rw2", constr = TRUE, hyper = list(prec = list(intitial = log(1.93), fixed = TRUE)))

# components = ~ -1 +
#   eta(z, model = "rw1", constr = TRUE, hyper = list(prec = list(intitial = log(1.93), fixed = TRUE)))

components = ~ -1 +
  eta(z, model = "rw1", constr = TRUE)

formula = y ~ eta
likelihood = like(formula = formula, family = "gaussian", data = data)

res.inlabru <- bru(components = components,
                   likelihood = likelihood, 
                   options = list(verbose = F, 
                                  bru_verbose = 1))

ggplot(data.frame(eta = res.inlabru$summary.random$eta$mean, z = z, sin.z = sin(z))) + geom_point(aes(x = z, y = eta)) + geom_point(aes(x = z, y = sin.z))

###   ----   Run stan analyses   ---- 
input_stan <- list(y = data$y)

# run first stan program: rw1_stepwise.stan:

fit_stepwise <- stan(
  file = "rw1_stepwise.stan",
  data = input_stan,
  chains = 4,
  iter = 80000,
  warmup = 8000,
  refresh = 1000,
  seed = 123
)

list_of_draws_stepwise <- rstan::extract(fit_stepwise)
eta_df_stepwise <- data.frame(list_of_draws_stepwise$eta)
summary_stepwise <- as.data.frame(summary(fit_stepwise)$summary) %>% mutate("Implementation" = "Stepwise")

trace_stepwise <- plot(fit_stepwise, plotfun = "trace", pars = c("tau_y", "tau_eta", "eta[1]", "eta[50]", "eta[90]", "eta[100]"))
ggsave("trace_stepwise.pdf", plot = trace_stepwise, dpi="retina", device = "pdf")

#  run second stan program: rw1_by_differences.stan

fit_diff_1 <- stan(
  file = "rw1_by_differences.stan",
  data = input_stan,
  chains = 4,
  iter = 80000,
  warmup = 8000,
  refresh = 10000,
  seed = 123
)

trace_diff_1 <- plot(fit_diff_1, plotfun = "trace", pars = c("tau_y", "tau_eta", "eta[1]", "eta[50]", "eta[90]", "eta[100]"))
ggsave("trace_diff_1.pdf", plot = trace_diff_1, dpi="retina", device = "pdf")

list_of_draws_diff_1 <- rstan::extract(fit_diff_1)
eta_df_diff_1 <- data.frame(list_of_draws_diff_1$eta)
summary_diff_1 <- as.data.frame(summary(fit_diff_1)$summary) %>% mutate("Implementation" = "Diff.1")

#  run third stan program: rw1_by_differences_icar.stan

fit_diff_2 <- stan(
  file = "rw1_by_differences_target.stan",
  data = input_stan,
  chains = 4,
  iter = 8000,
  warmup = 800,
  refresh = 10,
  seed = 123
)

trace_diff_2 <- plot(fit_diff_2, plotfun = "trace", pars = c("tau_y", "tau_eta", "eta[1]", "eta[50]", "eta[90]", "eta[100]"))
ggsave("trace_diff_2.pdf", plot = trace_diff_2, dpi="retina", device = "pdf")

list_of_draws_diff_2 <- rstan::extract(fit_diff_2)
eta_df_diff_2 <- data.frame(list_of_draws_diff_2$eta)
summary_diff_2 <- as.data.frame(summary(fit_diff_2)$summary) %>% mutate("Implementation" = "Diff.2")

###   ----   Compare results of different rw1 implementations   ----:


palette <- c('#70A4D4', '#ECC64B', '#93AD80', '#da9124', '#696B8D',
             '#3290c1', '#5d8060', '#D7B36A', '#826133', '#A85150')

eta_df_stepwise_long <- eta_df_stepwise %>% pivot_longer(cols = everything(), names_to = "Param", values_to = "Stepwise")
eta_df_diff_1_long <- eta_df_diff_1 %>% pivot_longer(cols = everything(), names_to = "Param", values_to = "Diff.1")
eta_df_diff_2_long <- eta_df_diff_2 %>% pivot_longer(cols = everything(), names_to = "Param", values_to = "Diff.2")

eta_all <- cbind(eta_df_stepwise_long, Diff.1 = eta_df_diff_1_long$Diff.1, Diff.2 = eta_df_diff_2_long$Diff.2) %>% 
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

eta.plot.full <- ggplot(eta_all) +
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.1, y = after_stat(density), color = "Difference 1", fill = "Difference 1"), bins = 500, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 500, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  facet_wrap(~Param)

ggsave("eta_histogram_full.pdf", plot = eta.plot.full, dpi="retina", device = "pdf")

#   Plot distributions of hyperparameters:
tau_eta_df <- data.frame(Stepwise = list_of_draws_stepwise$tau_eta,
                         Diff.1 = list_of_draws_diff_1$tau_eta,
                         Diff.2 = list_of_draws_diff_2$tau_eta)

p.tau_eta <- ggplot(tau_eta_df) + 
  geom_histogram(aes(x = Stepwise, color = "Stepwise", fill = "Stepwise", y = after_stat(density)), alpha = 0.5, bins = 100, size = 0) + 
  geom_histogram(aes(x = Diff.1, color = "Difference 1", fill = "Difference 1", y = after_stat(density)), alpha = 0.5, bins = 100, size = 0) + 
  geom_histogram(aes(x = Diff.2, color = "Difference 2", fill = "Difference 2", y = after_stat(density)), alpha = 0.5, bins = 100, size = 0) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Tau eta", x = "", y = "")

tau_y_df <- data.frame(Stepwise = list_of_draws_stepwise$tau_y,
                         Diff.1 = list_of_draws_diff_1$tau_y,
                         Diff.2 = list_of_draws_diff_2$tau_y)

p.tau_y <- ggplot(tau_y_df) + 
  geom_histogram(aes(x = Stepwise, color = "Stepwise", fill = "Stepwise", y = after_stat(density)), alpha = 0.5, bins = 100, size = 0) + 
  geom_histogram(aes(x = Diff.1, color = "Difference 1", fill = "Difference 1", y = after_stat(density)), alpha = 0.5, bins = 100, size = 0) + 
  geom_histogram(aes(x = Diff.2, color = "Difference 2", fill = "Difference 2", y = after_stat(density)), alpha = 0.5, bins = 100, size = 0) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Tau y", x = "", y = "")

p.tau  <- (p.tau_eta | p.tau_y)
  
ggsave("tau.pdf", plot = p.tau, dpi="retina", device = "pdf")

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

ggsave("summary_eta.pdf", plot = p.summary.eta, dpi="retina", device = "pdf")

#   ----    Plot comparison of Stan and inlabru   ----

eta_ib <- data.frame(res.inlabru$marginals.random$eta)

p.eta.1 <- ggplot(eta_all %>% filter(Param == 1)) +
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 100, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.1, y = after_stat(density), color = "Difference 1", fill = "Difference 1"), bins = 100, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 100, alpha = 0.5, size = 0) +
  geom_ares(data = data.frame(res.inlabru$marginals.random$eta$index.1), aes(x = x, y = y, color = "Inlabru"), alpha = 0.5) + 
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 1", x = "", y = "")

p.eta.20 <- ggplot(eta_all %>% filter(Param == 20)) +
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 100, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.1, y = after_stat(density), color = "Difference 1", fill = "Difference 1"), bins = 100, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 100, alpha = 0.5, size = 0) +
  geom_ares(data = data.frame(res.inlabru$marginals.random$eta$index.20), aes(x = x, y = y, color = "Inlabru"), alpha = 0.5) + 
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 20", x = "", y = "")

p.eta.40 <- ggplot(eta_all %>% filter(Param == 40)) +
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 100, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.1, y = after_stat(density), color = "Difference 1", fill = "Difference 1"), bins = 100, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 100, alpha = 0.5, size = 0) +
  geom_ares(data = data.frame(res.inlabru$marginals.random$eta$index.40), aes(x = x, y = y, color = "Inlabru"), alpha = 0.5) + 
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 40", x = "", y = "")

p.eta.60 <- ggplot(eta_all %>% filter(Param == 60)) +
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 100, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.1, y = after_stat(density), color = "Difference 1", fill = "Difference 1"), bins = 100, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 100, alpha = 0.5, size = 0) +
  geom_ares(data = data.frame(res.inlabru$marginals.random$eta$index.60), aes(x = x, y = y, color = "Inlabru"), alpha = 0.5) + 
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 60", x = "", y = "")

p.eta.80 <- ggplot(eta_all %>% filter(Param == 80)) +
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 100, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.1, y = after_stat(density), color = "Difference 1", fill = "Difference 1"), bins = 100, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 100, alpha = 0.5, size = 0) +
  geom_ares(data = data.frame(res.inlabru$marginals.random$eta$index.80), aes(x = x, y = y, color = "Inlabru"), alpha = 0.5) + 
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 80", x = "", y = "")

p.eta.100 <- ggplot(eta_all %>% filter(Param == 100)) +
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stepwise", fill = "Stepwise"), bins = 100, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.1, y = after_stat(density), color = "Difference 1", fill = "Difference 1"), bins = 100, alpha = 0.5, size = 0) +
  geom_histogram(aes(x = Diff.2, y = after_stat(density), color = "Difference 2", fill = "Difference 2"), bins = 100, alpha = 0.5, size = 0) +
  geom_ares(data = data.frame(res.inlabru$marginals.random$eta$index.100), aes(x = x, y = y, color = "Inlabru"), alpha = 0.5) + 
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 100", x = "", y = "")

p.eta.ib <- (p.eta.1 | p.eta.20 | p.eta.40)/(p.eta.60 | p.eta.80 | p.eta.100)

ggsave("eta_inlabru.pdf", plot = p.eta.ib, device = "pdf",
       path = path, height = 5, width = 8, dpi = "retina")






