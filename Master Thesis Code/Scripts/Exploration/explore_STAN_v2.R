# second attempt at successful model fit with STAN

library("rstan")
library(tidyverse)
library(ggplot2)

library(INLA)
library(inlabru)

seed = 324
set.seed(seed)

nx = 10
nt = 10

N = 500

at.risk = 1000

#obs <- data.frame(expand.grid(x=1:nx, t=1:nt))
obs = data.frame(x = sample(1:nx, N, replace = TRUE), t = sample(1:nx, N, replace = TRUE))

#  Standard deviation of 0.1 for the iid effects (beta)
tau.beta = 1/0.1**2

#  Standard deviation of 0.01 for the noise 
tau.epsilon = 1/0.01**2

kappa_true = 0.3*cos((1:nt)*pi/5)
kappa_true = kappa_true - mean(kappa_true)  #  sum to zero

# change this into an effect of x
alpha_true = cos(((1:nx - 3)* pi)/6)
alpha_true = alpha_true - mean(alpha_true)  # sum to zero

phi_true = -0.5  

#  sample synthetic data:
beta_true = rnorm(nx, 0, sqrt(1/tau.beta))
beta_true = 1/nx + beta_true - mean(beta_true)   # sum to one

data <- obs %>%
  mutate(alpha = alpha_true[x]) %>%
  mutate(beta = beta_true[x]) %>%
  mutate(phi = phi_true) %>%
  mutate(phi.t = phi*t) %>%
  mutate(kappa = kappa_true[t]) %>%
  mutate(epsilon = rnorm(N, mean=0, sd = sqrt(1/tau.epsilon))) %>%
  mutate(x.c = x, t.c = t, xt = seq_along(x)) %>%
  mutate(E = at.risk) %>%
  mutate(eta = alpha + beta*phi.t + beta*kappa + epsilon) %>%
  mutate(Y = rpois(N, E*exp(eta)))

#  helper values for constraining of beta:
A.mat = matrix(1, nrow = 1, ncol = 10)
e.vec = 1

# common uninformative prior distribution:
pc.prior <- list(prec = list(prior = "pc.prec", param = c(1,0.5)))

comp = ~ -1 + 
  Int(1) + 
  alpha(x, model = "rw1", values=1:nx, constr = TRUE, hyper = pc.prior) + 
  phi(t, model = "linear", prec.linear = 1) +
  beta(x.c, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior) + 
  kappa(t.c, model = "rw1", values = 1:nt, constr = TRUE, hyper = pc.prior) +
  epsilon(xt, model = "iid", hyper = pc.prior)

formula = Y ~ -1 + Int + alpha + beta*phi + beta*kappa + epsilon

likelihood = like(formula = formula, family = "poisson", data = data, E = data$E)

c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)  # control.compute

res.inlabru = bru(components = comp,
              likelihood, 
              options = list(verbose = F,
                             bru_verbose = 1, 
                             num.threads = "1:1",
                             control.compute = c.c,
                             bru_max_iter=30
              ))

# plot results from inlabru:

palette.basis <- c('#70A4D4', '#ECC64B', '#93AD80', '#da9124', '#696B8D',
                   '#3290c1',
                   '#5d8060', '#D7B36A', '#826133', '#A85150')

data.alpha = cbind(res.inlabru$summary.random$alpha, alpha.true = alpha[res.inlabru$summary.random$alpha$ID])
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

data.beta = cbind(res.inlabru$summary.random$beta, beta.true = beta[res.inlabru$summary.random$beta$ID])
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

data.kappa = cbind(res.inlabru$summary.random$kappa, kappa.true = kappa[res.inlabru$summary.random$kappa$ID])

p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
  geom_point(aes(y = kappa.true, color = "True value", fill = "True value")) + 
  geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
  scale_color_manual(name = "",
                     values = palette.basis ) +
  scale_fill_manual(name = "",
                    values = palette.basis ) +
  labs(x = "t", y = "kappa", title = "Kappa")

p.kappa

p.phi <- ggplot(data.frame(res.inlabru$marginals.fixed)) + 
  geom_area(aes(x = phi.x, y = phi.y, fill = "Estimated"), alpha = 0.4) + 
  geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[2], color = "Estimated", fill = "Estimated")) + 
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  labs(x = "Value of phi", y = " ", title = "Phi")
p.phi

data.eta <- data.frame({eta.sim = res.inlabru$summary.linear.predictor$mean[1:N]}) %>%
  mutate(true.eta = data$eta)
p.eta <- ggplot(data = data.eta) +
  geom_point(aes(x = eta.sim, y = true.eta), color = palette.basis[1]) + 
  labs(x="Estimated eta", y="True value for eta", title = "Eta")
p.eta
