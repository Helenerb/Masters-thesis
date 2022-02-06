#' Script for comparison of Stan implementations of random walk models
#' 

# on Markov:
#   ----   Load libraries and set workspace   ----

# TODO: Change to your directory
setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Random walk implementations/Rw1_constraints")

library("rstan")
library("inlabru")
library("tidyverse")
library("ggplot2")
library("patchwork")

# generate test data:
set.seed(123)
n=100
z=seq(0,6,length.out=n)
y=sin(z)+rnorm(n,mean=0,sd=0.5)
data=data.frame(y=y,z=z)

#   ----   Analysis with gamma priors    ----

run.inlabru <- function(input.data){
  components = ~ - 1 +
    Int(1, prec.linear = 0.001, mean.linear = 0) +
    eta(z, model = "rw1", constr = TRUE, scale.model = FALSE, hyper = list(prec = list(prior ="loggamma", param = c(1, 0.00005), initial = log(1))))
  
  formula = y ~ Int + eta
  likelihood = like(formula = formula, family = "gaussian", data = input.data)
  
  results <- bru(components = components,
                     likelihood = likelihood, 
                     options = list(verbose = F, 
                                    bru_verbose = 1))
  return(results)
}

res.inlabru <- run.inlabru(input.data = data)

ggplot(data.frame(eta = res.inlabru$summary.random$eta$mean, z = z, sin.z = sin(z))) + geom_point(aes(x = z, y = eta)) + geom_point(aes(x = z, y = sin.z))

###   ----   Run stan analyses   ---- 
input_stan <- list(y = data$y)

#run stan program - stan fit object stored in fit_diff_1
fit_diff_1 <- stan(
  file = "stand_alone_rw1_stepwise.stan",
  data = input_stan,
  chains = 4,
  iter = 30000,
  warmup = 3000,
  refresh = 500,
  seed = 123
)

# save trace plot
trace_diff_1 <- plot(fit_diff_1, plotfun = "trace", pars = c("tau_y", "tau_eta", "eta[1]", "eta[50]", "eta[90]", "eta[100]"))
ggsave("trace_diff_1.pdf", plot = trace_diff_1, device = "pdf", dpi = "retina", height = 5, width = 8)

# find list of HMC draws from stan
list_of_draws_diff_1 <- rstan::extract(fit_diff_1)
save(list_of_draws_diff_1, file = 'draws_diff_1.RData')
eta_df_diff_1 <- data.frame(list_of_draws_diff_1$eta)
summary_diff_1 <- as.data.frame(summary(fit_diff_1)$summary) %>% mutate("Implementation" = "Diff.1")

#   ---- Simple plots, comparing inlabru and diff.1   ----

# color palette for plotting
palette <- c('#70A4D4', '#ECC64B', '#93AD80', '#da9124', '#696B8D',
             '#3290c1', '#5d8060', '#D7B36A', '#826133', '#A85150')

#  plot precisions
tau_eta_df <- data.frame(Diff.1 = list_of_draws_diff_1$tau_eta, Diff.1.theta = list_of_draws_diff_1$theta_eta)
tau_y_df <- data.frame(Diff.1 = list_of_draws_diff_1$tau_y)

p.theta_eta_ib <- ggplot(tau_eta_df)  + 
  geom_histogram(aes(x = Diff.1.theta, color = "Stan", fill = "Stan", y = after_stat(density)), alpha = 0.3, bins = 100) + 
  geom_area(data = data.frame(res.inlabru$internal.marginals.hyperpar$`Log precision for eta`), aes(x=x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.3) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Theta eta", x = "", y = "")

p.tau_eta_ib <- ggplot(tau_eta_df) + 
  geom_histogram(aes(x = Diff.1, color = "Stan", fill = "Stan", y = after_stat(density)), alpha = 0.3, bins = 100) + 
  geom_area(data = data.frame(res.inlabru$marginals.hyperpar$`Precision for eta`) %>% filter(x < 40), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.3) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Tau eta", x = "", y = "")

p.tau_y_ib <- ggplot(tau_y_df) + 
  geom_histogram(aes(x = Diff.1, color = "Stan", fill = "Stan", y = after_stat(density)), alpha = 0.3, bins = 100) + 
  geom_area(data = data.frame(res.inlabru$marginals.hyperpar$`Precision for the Gaussian observations`) %>% filter(x < 15), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.3) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Tau y", x = "", y = "")

p.tau_y_ib

p.tau_ib  <- (p.tau_eta_ib | p.tau_y_ib) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("tau_ib.pdf", plot = p.tau_ib, dpi="retina", device = "pdf", height = 5, width = 8)

eta_df_diff_1_long <- eta_df_diff_1 %>% pivot_longer(cols = everything(), names_to = "Param", values_to = "Diff.1")
eta_all <- data.frame(Diff.1 = eta_df_diff_1_long$Diff.1, Param = eta_df_diff_1_long$Param) %>% 
  mutate(idx = parse_number(Param))

p.eta.1.1 <- ggplot(eta_all %>% filter(idx == 1)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.1), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.1, color = "Stan", fill = "Stan", y = after_stat(density)), alpha = 0.3, bins = 100, size  = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 1", x = "", y = "")

p.eta.20.1 <- ggplot(eta_all %>% filter(idx == 20)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.20), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.1, color = "Stan", fill = "Stan", y = after_stat(density)), alpha = 0.3, bins = 100, size  = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 20", x = "", y = "")

p.eta.40.1 <- ggplot(eta_all %>% filter(idx == 40)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.40), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.1, color = "Stan", fill = "Stan", y = after_stat(density)), alpha = 0.3, bins = 100, size  = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 40", x = "", y = "")

p.eta.60.1 <- ggplot(eta_all %>% filter(idx == 60)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.60), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.1, color = "Stan", fill = "Stan", y = after_stat(density)), alpha = 0.3, bins = 100, size  = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 60", x = "", y = "")

p.eta.80.1 <- ggplot(eta_all %>% filter(idx == 80)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.80), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.1, color = "Stan", fill = "Stan", y = after_stat(density)), alpha = 0.3, bins = 100, size  = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 80", x = "", y = "")

p.eta.100.1 <- ggplot(eta_all %>% filter(idx == 100)) +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.100), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(aes(x = Diff.1, color = "Stan", fill = "Stan", y = after_stat(density)), alpha = 0.3, bins = 100, size  = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Z = 100", x = "", y = "")

p.eta.ib.1 <- (p.eta.1.1 | p.eta.20.1 | p.eta.40.1)/(p.eta.60.1 | p.eta.80.1 | p.eta.100.1) + 
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("eta_inlabru_diff_1.pdf", plot = p.eta.ib.1, device = "pdf", height = 5, width = 8, dpi = "retina")

eta_summary_diff_1 <- summary_diff_1%>%
  rownames_to_column(var="idx") %>%
  filter(grepl("eta", idx)) %>% filter(!grepl("tau_eta", idx)) %>%
  filter(!grepl("diff", idx)) %>% filter(!grepl("theta_eta", idx)) %>% filter(!grepl("theta_y", idx)) %>%
  filter(!grepl("eta_scaled", idx)) %>%
  mutate(idx = parse_number(idx))

summary_eta <- data.frame(idx = eta_summary_diff_1$idx,
                          Diff.1 = eta_summary_diff_1$mean,
                          Diff.1.025 = eta_summary_diff_1$`2.5%`,
                          Diff.1.975 = eta_summary_diff_1$`97.5%`) %>%
  mutate(True = sin(z)[idx])


p.summary.eta.ib <- ggplot(summary_eta, aes(x = idx)) + 
  geom_ribbon(aes(ymin = Diff.1.025, ymax = Diff.1.975, fill = "Stan", color = "Stan"), alpha = 0.5) + 
  geom_ribbon(data = data.frame(res.inlabru$summary.random$eta) %>% mutate(idx = summary_eta$idx), aes(x = idx, ymin = X0.025quant, ymax = X0.975quant , color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_point(aes(y = Diff.1, color = "Stan", fill = "Stan")) + 
  geom_point(data = data.frame(res.inlabru$summary.random$eta) %>% mutate(idx = summary_eta$idx), aes(x = idx, y = mean, color = "Inlabru", fill = "Inlabru")) + 
  geom_point(aes(y = True, color = "True", fill = "True")) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs("Means of eta", x = "z", y = "")

ggsave("summary_eta_ib.pdf", plot = p.summary.eta.ib, dpi="retina", device = "pdf", height = 5, width = 8)

p.intercept <- ggplot() + 
  geom_histogram(data = data.frame(Diff.1 = list_of_draws_diff_1$intercept), aes(x = Diff.1, color = "Stan", fill = "Stan", y = after_stat(density)), bins = 100, alpha = 0.5) + 
  geom_area(data = data.frame(res.inlabru$marginals.fixed$Int), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs("Intercept", x = "", y = "")

ggsave("intercept.pdf", plot = p.intercept, dpi="retina", device = "pdf", height = 5, width = 8)

#   ----   Analysis with fixed precisions    ----

# run.inlabru.fh <- function(input.data){
#   components = ~ - 1 +
#     Int(1, prec.linear = 0.001, mean.linear = 0) +
#     eta(z, model = "rw1", constr = TRUE, scale.model = FALSE, hyper = list(prec = list(initial = 2, fixed = T, prior ="loggamma", param = c(1, 0.00005))))
#   
#   formula = y ~ Int + eta
#   likelihood = like(formula = formula, family = "gaussian", data = input.data)
#   
#   results <- bru(components = components,
#                  likelihood = likelihood, 
#                  options = list(verbose = F, 
#                                 bru_verbose = 1))
#   return(results)
# }
# 
# res.inlabru.fh <- run.inlabru.fh(input.data = data)
# 
# ggplot(data.frame(eta = res.inlabru.fh$summary.random$eta$mean, z = z, sin.z = sin(z))) + geom_point(aes(x = z, y = eta)) + geom_point(aes(x = z, y = sin.z))
# 
# ###   ----   Run stan analyses   ---- 
# input_stan <- list(y = data$y)
# 
# fit_diff_1_fh <- stan(
#   file = "stand_alone_rw1_stepwise_fh.stan",
#   data = input_stan,
#   chains = 4,
#   iter = 10000,
#   warmup = 1000,
#   refresh = 500,
#   seed = 123
# )
# 
# trace_diff_1_fh <- plot(fit_diff_1_fh, plotfun = "trace", pars = c("tau_y", "tau_eta", "eta[1]", "eta[50]", "eta[90]", "eta[100]"))
# ggsave("trace_diff_1_fh.pdf", plot = trace_diff_1_fh, dpi="retina", device = "pdf")
# 
# list_of_draws_diff_1_fh <- rstan::extract(fit_diff_1_fh)
# save(list_of_draws_diff_1_fh, file = 'draws_diff_1_fh.RData')
# eta_df_diff_1_fh <- data.frame(list_of_draws_diff_1_fh$eta)
# summary_diff_1_fh <- as.data.frame(summary(fit_diff_1_fh)$summary) %>% mutate("Implementation" = "Diff.1")
# 
# #   ---- Simple plots, comparing inlabru and diff.1   ----
# 
# palette <- c('#70A4D4', '#ECC64B', '#93AD80', '#da9124', '#696B8D',
#              '#3290c1', '#5d8060', '#D7B36A', '#826133', '#A85150')
# 
# eta_df_diff_1_long_fh <- eta_df_diff_1_fh %>% pivot_longer(cols = everything(), names_to = "Param", values_to = "Diff.1")
# eta_all_fh <- data.frame(Diff.1 = eta_df_diff_1_long_fh$Diff.1, Param = eta_df_diff_1_long_fh$Param) %>% 
#   mutate(idx = parse_number(Param))
# 
# p.eta.1.1_fh <- ggplot(eta_all_fh %>% filter(idx == 1)) +
#   geom_area(data = data.frame(res.inlabru.fh$marginals.random$eta$index.1), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
#   geom_density(aes(x = Diff.1, color = "Stan", fill = "Stan"), alpha = 0.3) +
#   scale_fill_manual(name = "", values = palette) + 
#   scale_color_manual(name = "", values = palette) + 
#   theme_classic() + 
#   labs(title = "Z = 1", x = "", y = "")
# 
# p.eta.20.1_fh <- ggplot(eta_all_fh %>% filter(idx == 20)) +
#   geom_area(data = data.frame(res.inlabru.fh$marginals.random$eta$index.20), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
#   geom_density(aes(x = Diff.1, color = "Stan", fill = "Stan"), alpha = 0.3) +
#   scale_fill_manual(name = "", values = palette) + 
#   scale_color_manual(name = "", values = palette) + 
#   theme_classic() + 
#   labs(title = "Z = 20", x = "", y = "")
# 
# p.eta.40.1_fh <- ggplot(eta_all_fh %>% filter(idx == 40)) +
#   geom_area(data = data.frame(res.inlabru.fh$marginals.random$eta$index.40), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
#   geom_density(aes(x = Diff.1, color = "Stan", fill = "Stan"), alpha = 0.3) +
#   scale_fill_manual(name = "", values = palette) + 
#   scale_color_manual(name = "", values = palette) + 
#   theme_classic() + 
#   labs(title = "Z = 40", x = "", y = "")
# 
# p.eta.60.1_fh <- ggplot(eta_all_fh %>% filter(idx == 60)) +
#   geom_area(data = data.frame(res.inlabru.fh$marginals.random$eta$index.60), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
#   geom_density(aes(x = Diff.1, color = "Stan", fill = "Stan"), alpha = 0.3) +
#   scale_fill_manual(name = "", values = palette) + 
#   scale_color_manual(name = "", values = palette) + 
#   theme_classic() + 
#   labs(title = "Z = 60", x = "", y = "")
# 
# p.eta.80.1_fh <- ggplot(eta_all_fh %>% filter(idx == 80)) +
#   geom_area(data = data.frame(res.inlabru.fh$marginals.random$eta$index.80), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
#   geom_density(aes(x = Diff.1, color = "Stan", fill = "Stan"), alpha = 0.3) +
#   scale_fill_manual(name = "", values = palette) + 
#   scale_color_manual(name = "", values = palette) + 
#   theme_classic() + 
#   labs(title = "Z = 80", x = "", y = "")
# 
# p.eta.100.1_fh <- ggplot(eta_all_fh %>% filter(idx == 100)) +
#   geom_area(data = data.frame(res.inlabru.fh$marginals.random$eta$index.100), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
#   geom_density(aes(x = Diff.1, color = "Stan", fill = "Stan"), alpha = 0.3) +
#   scale_fill_manual(name = "", values = palette) + 
#   scale_color_manual(name = "", values = palette) + 
#   theme_classic() + 
#   labs(title = "Z = 100", x = "", y = "")
# 
# p.eta.ib.1_fh <- (p.eta.1.1_fh | p.eta.20.1_fh | p.eta.40.1_fh)/(p.eta.60.1_fh | p.eta.80.1_fh | p.eta.100.1_fh) + 
#   plot_layout(guides = "collect")
# 
# ggsave("eta_inlabru_diff_1_fh.pdf", plot = p.eta.ib.1_fh, device = "pdf", height = 5, width = 8, dpi = "retina")
# 
# eta_summary_diff_1_fh <- summary_diff_1_fh%>%
#   rownames_to_column(var="idx") %>%
#   filter(grepl("eta", idx)) %>% filter(!grepl("tau_eta", idx)) %>%
#   filter(!grepl("diff", idx)) %>% filter(!grepl("theta_eta", idx)) %>%
#   filter(!grepl("eta_scaled", idx)) %>%
#   mutate(idx = parse_number(idx))
# 
# summary_eta_fh <- data.frame(idx = eta_summary_diff_1_fh$idx,
#                           Diff.1 = eta_summary_diff_1_fh$mean,
#                           Diff.1.025 = eta_summary_diff_1_fh$`2.5%`,
#                           Diff.1.975 = eta_summary_diff_1_fh$`97.5%`) %>%
#   mutate(True = sin(z)[idx])
# 
# 
# p.summary.eta.ib_fh <- ggplot(summary_eta_fh, aes(x = idx)) + 
#   geom_ribbon(aes(ymin = Diff.1.025, ymax = Diff.1.975, fill = "Stan", color = "Stan"), alpha = 0.5) + 
#   geom_ribbon(data = data.frame(res.inlabru.fh$summary.random$eta) %>% mutate(idx = summary_eta_fh$idx), aes(x = idx, ymin = X0.025quant, ymax = X0.975quant , color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
#   geom_point(aes(y = Diff.1, color = "Stan", fill = "Stan")) + 
#   geom_point(data = data.frame(res.inlabru.fh$summary.random$eta) %>% mutate(idx = summary_eta_fh$idx), aes(x = idx, y = mean, color = "Inlabru", fill = "Inlabru")) + 
#   geom_point(aes(y = True, color = "True", fill = "True")) + 
#   theme_classic() + 
#   scale_color_manual(name = "", values = palette) + 
#   scale_fill_manual(name = "", values = palette) + 
#   labs("Means of eta - fixed hyper", x = "z", y = "")
# 
# ggsave("summary_eta_ib_fh.pdf", plot = p.summary.eta.ib_fh, dpi="retina", device = "pdf", height = 5, width = 8)
# 
# p.intercept_fh <- ggplot() + 
#   geom_density(data = data.frame(Diff.1 = list_of_draws_diff_1_fh$intercept), aes(x = Diff.1, color = "Stan", fill = "Stan"), alpha = 0.5) + 
#   geom_area(data = data.frame(res.inlabru.fh$marginals.fixed$Int), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
#   theme_classic() + 
#   scale_color_manual(name = "", values = palette) + 
#   scale_fill_manual(name = "", values = palette) + 
#   labs("Intercept - fixed hyper", x = "", y = "")
# 
# ggsave("intercept_fh.pdf", plot = p.intercept_fh, dpi="retina", device = "pdf", height = 5, width = 8)
