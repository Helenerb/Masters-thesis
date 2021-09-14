# script used for exploring the STAN library for modelling with
# Hamiltonian MC

library("rstan")
library(tidyverse)

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
  warmup = 2000,
  iter = 10000,
  refresh = 500,
  seed=123,
  init = intit_f
)

print(fit)
plot(fit)
fit_df <- as.data.frame(fit)
fit_summary <- summary(fit)
print(fit_summary$summary)
fir_summary_df <- as.data.frame(fit_summary$summary)

summary_alpha <- fir_summary_df %>%
  rownames_to_column("parameter") %>%
  filter(grepl('alpha', parameter)) %>%
  filter(!grepl('sigma_alpha', parameter)) %>%
  extract(parameter, into="index", regex="([0-9]+)")
  