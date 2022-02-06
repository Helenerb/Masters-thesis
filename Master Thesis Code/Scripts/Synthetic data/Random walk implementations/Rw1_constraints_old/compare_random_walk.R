#' Script for comparison of Stan implementations of random walk models
#' 


# on Markov:
#   ----   Load libraries and set workspace   ----
set_workdirectory <- function(markov=TRUE){
  if(markov){
    .libPaths("~/Documents/R_libraries")
    setwd("~/Documents/GitHub/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Random walk implementations/Rw1_constraints")
  } else {
    setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Random walk implementations/Rw1_constraints")
  }
}

set_workdirectory(markov=F)

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

components = ~ -1 +
  Int(prec.linear = 0.001, mean.linear = 0) + 
  eta(z, model = "rw1", constr = TRUE, scale.model = F, hyper = list(prec = list(prior ="loggamma", param = c(1, 0.00005))))

formula = y ~ Int + eta
likelihood = like(formula = formula, family = "gaussian", data = data)

res.inlabru <- bru(components = components,
                   likelihood = likelihood, 
                   options = list(verbose = F, 
                                  bru_verbose = 1))

ggplot(data.frame(eta = res.inlabru$summary.random$eta$mean, z = z, sin.z = sin(z))) + geom_point(aes(x = z, y = eta)) + geom_point(aes(x = z, y = sin.z))

###   ----   Run stan analyses   ---- 
input_stan <- list(y = data$y)

# run first stan program: rw1_stepwise.stan:

system.time({
  fit_stepwise <- stan(
    file = "rw1_stepwise.stan",
    data = input_stan,
    chains = 4,
    iter = 10000,
    warmup = 2000,
    refresh = 2000,
    seed = 123
  )
}
)

list_of_draws_stepwise <- rstan::extract(fit_stepwise)
save(list_of_draws_stepwise, file = 'draws_stepwise.RData')
load(file = 'draws_stepwise.RData')
eta_df_stepwise <- data.frame(list_of_draws_stepwise$eta)
summary_stepwise <- as.data.frame(summary(fit_stepwise)$summary) %>% mutate("Implementation" = "Stepwise")

trace_stepwise <- plot(fit_stepwise, plotfun = "trace", pars = c("tau_y", "tau_eta", "eta[1]", "eta[50]", "theta_eta", "theta_y"))
ggsave("trace_stepwise.pdf", plot = trace_stepwise, dpi="retina", device = "pdf", height = 5, width = 8)


###   ----   Compare results of different rw1 implementations   ----:

palette <- c('#70A4D4', '#ECC64B', '#93AD80', '#da9124', '#696B8D',
             '#3290c1', '#5d8060', '#D7B36A', '#826133', '#A85150')

eta_df_stepwise_long <- eta_df_stepwise %>% pivot_longer(cols = everything(), names_to = "Param", values_to = "Stepwise")

#   ----   Plot distributions of hyperparameters:   ----
tau_eta_df <- data.frame(Stepwise = list_of_draws_stepwise$tau_eta)

tau_y_df <- data.frame(Stepwise = list_of_draws_stepwise$tau_y)

theta_eta_df <- data.frame(Stepwise = list_of_draws_stepwise$theta_eta)

theta_y_df <- data.frame(Stepwise = list_of_draws_stepwise$theta_y)

p.tau_eta_ib <- ggplot(tau_eta_df) + 
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.2, bins = 100) + 
  geom_area(data = data.frame(res.inlabru$marginals.hyperpar$`Precision for eta`) %>% filter(x < 25), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Tau eta", x = "", y = "")

p.tau_y_ib <- ggplot(tau_y_df) + 
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.2, bins = 100) + 
  geom_area(data = data.frame(res.inlabru$marginals.hyperpar$`Precision for the Gaussian observations`) %>% filter(x < 15), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Tau y", x = "", y = "")

p.theta_eta_ib <- ggplot(theta_eta_df) + 
  geom_histogram(aes(x = Stepwise, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.2, bins = 100) + 
  geom_area(data = data.frame(res.inlabru$internal.marginals.hyperpar$`Log precision for eta`), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Theta eta", x = "", y = "")

p.theta_y_ib <- ggplot(theta_y_df) + 
  geom_histogram(aes(x = Stepwise, color = "Stan", fill = "Stan", y = after_stat(density)), alpha = 0.2, bins = 100) + 
  geom_area(data = data.frame(res.inlabru$internal.marginals.hyperpar$`Log precision for the Gaussian observations`), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Theta y", x = "", y = "")

p.tau_ib  <- (p.tau_eta_ib | p.tau_y_ib) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
p.theta_ib <- (p.theta_eta_ib | p.theta_y_ib) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("tau_ib.pdf", plot = p.tau_ib, dpi="retina", device = "pdf", height = 5, width = 8)
ggsave("theta_ib.pdf", plot = p.theta_ib, dpi="retina", device = "pdf", height = 5, width = 8)

#   ----    Plot comparison of the precisions for each of the methods   ----

# stepwise
p.tau_eta_sw  <- ggplot(tau_eta_df) + 
  geom_histogram(aes(x = Stepwise, color = "Stan", fill = "Stan", y = after_stat(density)), alpha = 0.2, bins = 100) + 
  geom_area(data = data.frame(res.inlabru$marginals.hyperpar$`Precision for eta`) %>% filter(x < 25), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  theme_classic() + 
  scale_color_manual(name = "", values = c(palette[2], palette[1])) + 
  scale_fill_manual(name = "", values = c(palette[2], palette[1])) + 
  labs(title = "Tau eta", x = "", y = "")

p.tau_y_sw <- ggplot(tau_y_df) + 
  geom_histogram(aes(x = Stepwise, color = "Stan", fill = "Stan", y = after_stat(density)), alpha = 0.2, bins = 100) + 
  geom_area(data = data.frame(res.inlabru$marginals.hyperpar$`Precision for the Gaussian observations`) %>% filter(x < 15), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  theme_classic() + 
  scale_color_manual(name = "", values = c(palette[2], palette[1])) + 
  scale_fill_manual(name = "", values = c(palette[2], palette[1])) + 
  labs(title = "Tau y", x = "", y = "")

p.theta_eta_sw  <- ggplot(theta_eta_df) + 
  geom_histogram(aes(x = Stepwise, color = "Stan", fill = "Stan", y = after_stat(density)), alpha = 0.2, bins = 100) + 
  geom_area(data = data.frame(res.inlabru$internal.marginals.hyperpar$`Log precision for eta`) %>% filter(x < 25), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  theme_classic() + 
  scale_color_manual(name = "", values = c(palette[2], palette[1])) + 
  scale_fill_manual(name = "", values = c(palette[2], palette[1])) + 
  labs(title = "Theta eta", x = "", y = "")

p.theta_y_sw <- ggplot(theta_y_df) + 
  geom_histogram(aes(x = Stepwise, color = "Stan", fill = "Stan", y = after_stat(density)), alpha = 0.2, bins = 100) + 
  geom_area(data = data.frame(res.inlabru$internal.marginals.hyperpar$`Log precision for the Gaussian observations`) %>% filter(x < 15), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  theme_classic() + 
  scale_color_manual(name = "", values = c(palette[2], palette[1])) + 
  scale_fill_manual(name = "", values = c(palette[2], palette[1])) + 
  labs(title = "Theta y", x = "", y = "")

p.tau_sw <- (p.tau_eta_sw | p.tau_y_sw) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
p.theta_sw <- (p.theta_eta_sw | p.theta_y_sw) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("tau_ib_sw.pdf", plot = p.tau_sw, dpi="retina", device = "pdf", height = 5, width = 8)
ggsave("theta_ib_sw.pdf", plot = p.theta_sw, dpi="retina", device = "pdf", height = 5, width = 8)


#   ----    Plot comparison of Stan and inlabru   ----

eta_ib <- data.frame(res.inlabru$marginals.random$eta)

p.eta.1 <- ggplot() +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.1), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(data = eta_df_stepwise, aes(x = X1, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 100, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Eta[1]", x = "", y = "")

p.eta.20 <- ggplot() +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.20), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(data = eta_df_stepwise, aes(x = X20, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 100, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Eta[20]", x = "", y = "")

p.eta.40 <- ggplot() +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.40), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(data = eta_df_stepwise, aes(x = X40, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 100, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Eta[40]", x = "", y = "")

p.eta.60 <- ggplot() +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.60), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(data = eta_df_stepwise, aes(x = X60, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 100, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Eta[60]", x = "", y = "")

p.eta.80 <- ggplot() +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.80), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(data = eta_df_stepwise, aes(x = X80, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 100, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Eta[80]", x = "", y = "")

p.eta.100 <- ggplot() +
  geom_area(data = data.frame(res.inlabru$marginals.random$eta$index.100), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(data = eta_df_stepwise, aes(x = X100, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 100, alpha = 0.5, size = 0) +
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Eta[100]", x = "", y = "")

p.eta.ib <- (p.eta.1 | p.eta.20 | p.eta.40)/(p.eta.60 | p.eta.80 | p.eta.100) + 
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("eta_inlabru.pdf", plot = p.eta.ib, device = "pdf", height = 5, width = 8, dpi = "retina")

#   ----   Plot intercept    ----
p.intercept <- ggplot() + 
  #geom_area(data = data.frame(res.inlabru$marginals.fixed$Int), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(data = data.frame(int = list_of_draws_stepwise$intercept), aes(x = int, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins = 100, size  = 0) + 
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  theme_classic() + 
  labs(title = "Intercept", x = "", y = "")

ggsave("intercept.pdf", plot = p.intercept, device = "pdf", height = 4, width = 6.4, dpi = "retina")
