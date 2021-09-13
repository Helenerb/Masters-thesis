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
  mutate(Y = E*exp(eta))
  
# Generate input data in the form of a named list:

exploration_data <- list(
  E = true_data$E,
  Y = true_data$Y
)


