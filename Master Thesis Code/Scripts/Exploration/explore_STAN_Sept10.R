# script used for exploring the STAN library for modelling with
# Hamiltonian MC

# Note: The configuration in this script has convergence issues, 
# and displays identifiability problems. Return to figure out why. 

library("rstan")
library(tidyverse)
library(ggplot2)

library(INLA)
library(inlabru)

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
  mutate(Y = rpois(100,rate)) %>%
  mutate(x.c = x) %>%
  mutate(t.c = t) %>%
  mutate(xts = seq_along(x))

#   ---- For comparison - run the analysis with inlabru   ---- 
pc.prior <- list(prec = list(prior = "pc.prec", param = c(1,0.05)))

A.mat = matrix(1, nrow = 1, ncol = 10) 
e.vec = 1

#in general: unsure about scale.model...
# also unsure about the ~-1
comp.inlabru = ~ -1 +
  alpha(x, model = "rw1", values = 1:X, constr = FALSE, hyper = pc.prior, scale.model = TRUE) + # ikke helt sikker på om du skal ha med scale.model...
  beta(x.c, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior) +
  phi(t, model = "linear", prec.linear = 1) +
  kappa(t.c, model = "rw1", values = 1:T, constr = TRUE, hyper = pc.prior, scale.model = TRUE) +
  epsilon(xts, model = "iid", hyper = pc.prior)

formula = Y ~ -1 + alpha + beta*phi + beta*kappa + epsilon
likelihood = like(formula = formula, family = "poisson", data = true_data, E = true_data$E)

c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)

results.inlabru = bru(components = comp.inlabru,
               likelihood,
               options = list(verbose = F,
                              bru_verbose = 1, 
                              num.threads = "1:1",
                              control.compute = c.c,
                              control.predictor = list(link = 1),  # The log link function
                              bru_max_iter = 20
               ))

# plot the results:

palette.basis <- c('#70A4D4', '#ECC64B', '#93AD80', '#da9124', '#696B8D',
                   '#3290c1',
                   '#5d8060', '#D7B36A', '#826133', '#A85150')

data.alpha = cbind(results.inlabru$summary.random$alpha, alpha.true = true_alpha[results.inlabru$summary.random$alpha$ID])
p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
  geom_point(aes(y = alpha.true, color = "True value", fill = "True value")) + 
  geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
  scale_color_manual(name = "",
                     values = palette.basis ) +
  scale_fill_manual(name = "",
                    values = palette.basis ) +
  labs(title="Alpha", x = "x", y='')

p.alpha

data.beta = cbind(results.inlabru$summary.random$beta, beta.true = true_beta[results.inlabru$summary.random$beta$ID])
p.beta <- ggplot(data = data.beta, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
  geom_point(aes(y = beta.true, color = "True value", fill = "True value")) + 
  geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
  scale_color_manual(name = "",
                     values = palette.basis ) +
  scale_fill_manual(name = "",
                    values = palette.basis ) +
  labs(x = "x", y = "beta", title = "Beta")

p.beta

data.kappa = cbind(results.inlabru$summary.random$kappa, kappa.true = true_kappa[results.inlabru$summary.random$kappa$ID]) %>%
  mutate(kappa.drifted = mean + ID*results.inlabru$summary.fixed$mean[1])

p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
  geom_point(aes(y = kappa.true, color = "True value", fill = "True value")) + 
  geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
  geom_point(aes(y = kappa.drifted, color="Estimated drift")) + 
  scale_color_manual(name = "",
                     values = palette.basis ) +
  scale_fill_manual(name = "",
                    values = palette.basis ) +
  labs(x = "t", y = "kappa", title = "Kappa")

p.kappa

p.phi <- ggplot(data.frame(results.inlabru$marginals.fixed)) + 
  geom_area(aes(x = phi.x, y = phi.y, fill = "Estimated"), alpha = 0.4) + 
  geom_vline(data = results.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Estimated", fill = "Estimated")) + 
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  labs(x = "Value of phi", y = " ", title = "Phi")
p.phi

data.eta <- data.frame({eta.sim = results.inlabru$summary.linear.predictor$mean[1:100]}) %>%
  mutate(true.eta = true_data$eta)
p.eta <- ggplot(data = data.eta) +
  geom_point(aes(x = eta.sim, y = true.eta), color = palette.basis[1]) + 
  labs(x="Estimated eta", y="True value for eta", title = "Eta")
p.eta

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
