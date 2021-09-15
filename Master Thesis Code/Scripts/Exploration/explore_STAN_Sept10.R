# script used for exploring the STAN library for modelling with
# Hamiltonian MC

library("rstan")
library(tidyverse)
library(ggplot2)

# setting up data for testing:
T = 10  # number of time steps
X = 10  # number of age groups

obs <- data.frame(expand.grid(x=1:X, t=1:T))

# underlying functions for age and period effects:
true_alpha = cos(1:X/2 + 2) - 3
plot(true_alpha)

true_beta = rnorm(10, mean=0, sd=0.05)
true_beta = true_beta/sum(true_beta)
plot(true_beta)

true_kappa = (1:T)^2/100
true_kappa = true_kappa - sum(true_kappa)/length(true_kappa)
plot(true_kappa)

true_data = obs %>%
  mutate(alpha = true_alpha[x]) %>%
  mutate(beta = true_beta[x]) %>%
  mutate(kappa = true_kappa[t]) %>%
  mutate(epsilon = rnorm(T*X, mean=0, sd=0.01)) %>%
  mutate(eta = alpha + beta*kappa + epsilon) %>%
  mutate(E = rep(1000, 100)) %>%
  mutate(rate = E*exp(eta)) %>%
  mutate(Y = rpois(100,rate))
  
# Generate input data in the form of a named list:

exploration_data <- list(
  X=X,
  T=T,
  x=true_data$x,
  t=true_data$t,
  E = true_data$E,
  Y = true_data$Y
)

intit_f <- function() list(sigma_alpha = 0.1,
                           sigma_beta = 0.05,
                           sigma_kappa = 0.1,
                           sigma_epsilon = 0.01,
                           alpha = true_alpha,
                           beta = true_beta,
                           kappa = true_kappa,
                           eta = true_data$eta)

# run stan
fit <- stan(
  file="explore_STAN.stan",
  data = exploration_data,
  chains=4,
  warmup = 5000,
  iter = 20000,
  refresh = 500,
  seed=123#,
  #init = intit_f
)

pairs(fit, pars=c("alpha[1]", "beta[1]", "kappa[1]", "phi", "eta[1]"), condition = "accept_stat__")

print(fit)
plot(fit)
fit_df <- as.data.frame(fit)
fit_summary <- summary(fit)
print(fit_summary$summary)
fir_summary_df <- as.data.frame(fit_summary$summary)

# look at alpha
summary_alpha <- fir_summary_df %>%
  rownames_to_column("parameter") %>%
  filter(grepl('alpha', parameter)) %>%
  filter(!grepl('sigma_alpha', parameter)) %>%
  extract(parameter, into="index", regex="([0-9]+)") %>%
  mutate(index = parse_number(index)) %>%
  mutate(true_alpha = true_alpha[index])

plot_alpa <- ggplot(data=summary_alpha) +
  geom_point(aes(x=index, y=mean, color="estimated")) + 
  geom_line(aes(x=index, y=`2.5%`, color="estimated")) + 
  geom_line(aes(x=index, y=`97.5%`, color="estimated")) +
  geom_point(aes(x=index, y=true_alpha, color="true value")) +
  ggtitle("Alpha")
plot_alpa  

# look at beta
summary_beta <- fir_summary_df %>%
  rownames_to_column("parameter") %>%
  filter(grepl('beta', parameter)) %>%
  filter(!grepl('sigma_beta', parameter)) %>%
  extract(parameter, into="index", regex="([0-9]+)") %>%
  mutate(index = parse_number(index)) %>%
  mutate(true_beta = true_beta[index])

plot_beta <- ggplot(data=summary_beta) +
  geom_point(aes(x=index, y=mean, color="estimated")) + 
  geom_line(aes(x=index, y=`2.5%`, color="estimated")) + 
  geom_line(aes(x=index, y=`97.5%`, color="estimated")) +
  geom_point(aes(x=index, y=true_beta, color="true value")) +
  ggtitle("Beta")
plot_beta

# look at kappa
summary_kappa <- fir_summary_df %>%
  rownames_to_column("parameter") %>%
  filter(grepl('kappa', parameter)) %>%
  filter(!grepl('sigma_kappa', parameter)) %>%
  extract(parameter, into="index", regex="([0-9]+)") %>%
  mutate(index = parse_number(index)) %>%
  mutate(true_kappa = true_kappa[index])

plot_kappa <- ggplot(data=summary_kappa) +
  geom_point(aes(x=index, y=mean, color="estimated")) + 
  geom_line(aes(x=index, y=`2.5%`, color="estimated")) + 
  geom_line(aes(x=index, y=`97.5%`, color="estimated")) +
  geom_point(aes(x=index, y=true_kappa, color="true value")) +
  ggtitle("Kappa")
plot_kappa

# look at eta
summary_eta <- fir_summary_df %>%
  rownames_to_column("parameter") %>%
  filter(grepl('eta', parameter)) %>%
  filter(!grepl('beta', parameter)) %>%
  extract(parameter, into="index", regex="([0-9]+)") %>%
  mutate(index = parse_number(index)) %>%
  mutate(true_eta = true_data$eta)

plot_eta <- ggplot(data=summary_eta) +
  geom_point(aes(x=index, y=mean, color="estimated")) + 
  geom_line(aes(x=index, y=`2.5%`, color="estimated")) + 
  geom_line(aes(x=index, y=`97.5%`, color="estimated")) +
  geom_point(aes(x=index, y=true_eta, color="true value")) +
  ggtitle("Eta")
plot_eta
