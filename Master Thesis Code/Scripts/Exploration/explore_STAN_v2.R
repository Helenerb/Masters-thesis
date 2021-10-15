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

#N = 500
N = nx*nt

at.risk = 1000

obs <- data.frame(expand.grid(x=1:nx, t=1:nt))
#obs = data.frame(x = sample(1:nx, N, replace = TRUE), t = sample(1:nx, N, replace = TRUE))

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
pc.prior <- list(prec = list(prior = "pc.prec", param = c(1,0.05)))

comp = ~ -1 + 
  Int(1) + 
  alpha(x, model = "rw1", values=1:nx, constr = TRUE, hyper = pc.prior, scale.model=TRUE) + 
  phi(t, model = "linear", prec.linear = 1) +
  beta(x.c, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior) + 
  kappa(t.c, model = "rw1", values = 1:nt, constr = TRUE, hyper = pc.prior, scale.model=TRUE) +
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

p.int <- ggplot(data.frame(res.inlabru$marginals.fixed)) + 
  geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
  geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Estimated")) + 
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  labs(x = "Value of intercept", y = " ", title = "Intercept - inlabru")
p.int

data.alpha = cbind(res.inlabru$summary.random$alpha, alpha.true = alpha_true[res.inlabru$summary.random$alpha$ID])
p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
  geom_point(aes(y = alpha.true, color = "True value", fill = "True value")) + 
  geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
  scale_color_manual(name = "",
                     values = palette.basis ) +
  scale_fill_manual(name = "",
                    values = palette.basis ) +
  labs(title="Alpha - inlabru", x = "x", y='')

p.alpha

data.beta = cbind(res.inlabru$summary.random$beta, beta.true = beta_true[res.inlabru$summary.random$beta$ID])
p.beta <- ggplot(data = data.beta, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
  geom_point(aes(y = beta.true, color = "True value", fill = "True value")) + 
  geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
  scale_color_manual(name = "",
                     values = palette.basis ) +
  scale_fill_manual(name = "",
                    values = palette.basis ) +
  labs(x = "x", y = "beta", title = "Beta - inlabru")

p.beta

data.kappa = cbind(res.inlabru$summary.random$kappa, kappa.true = kappa_true[res.inlabru$summary.random$kappa$ID])

p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
  geom_point(aes(y = kappa.true, color = "True value", fill = "True value")) + 
  geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
  scale_color_manual(name = "",
                     values = palette.basis ) +
  scale_fill_manual(name = "",
                    values = palette.basis ) +
  labs(x = "t", y = "kappa", title = "Kappa - inlabru")

p.kappa

p.phi <- ggplot(data.frame(res.inlabru$marginals.fixed)) + 
  geom_area(aes(x = phi.x, y = phi.y, fill = "Estimated"), alpha = 0.4) + 
  geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[2], color = "Estimated", fill="Estimated")) + 
  geom_vline(aes(xintercept = phi_true, color="True", fill="True")) +
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  labs(x = "Value of phi", y = " ", title = "Phi - inlabru")
p.phi

data.eta <- data.frame({eta.sim = res.inlabru$summary.linear.predictor$mean[1:N]}) %>%
  mutate(true.eta = data$eta)
p.eta <- ggplot(data = data.eta) +
  geom_point(aes(x = eta.sim, y = true.eta), color = palette.basis[1]) + 
  labs(x="Estimated eta", y="True value for eta", title = "Eta")
p.eta

p.eta.2 <- ggplot(data = data.eta) +
  geom_point(aes(x=data$xt, y = eta.sim, color="Estimated")) +
  geom_point(data=data, aes(x=xt, y = eta, color="True")) +
  labs(x=" ", y="Eta", title="Eta- inlabru")
p.eta.2

#   ----   Implement the model in STAN   ---- 

exploration_data <- list(
  X=nx,
  T=nt,
  x=data$x,
  t=data$t,
  E = data$E,
  Y = data$Y
)

fit <- stan(
  file="explore_STAN_v2.stan",
  data = exploration_data,
  chains=4,
  warmup = 3000,
  iter = 30000,
  refresh = 1000,
  seed=123
)

traceplot(fit, pars=c("eta[1]", "eta[2]", "eta[3]", "eta[4]", "eta[5]",
                      "eta[6]", "eta[7]", "eta[8]", "eta[9]", "eta[10]" ))

traceplot(fit, pars=c("kappa[1]", "kappa[2]", "kappa[3]", "kappa[4]",
                      "kappa[5]", "kappa[6]", "kappa[7]", "kappa[8]",
                      "kappa[9]", "kappa[10]" ))

traceplot(fit)

#   ----   View results from model fitting using STAN   ----
pairs(fit, pars=c("sigma_alpha", "sigma_beta", "sigma_kappa", "sigma_epsilon", "phi", "intercept", "alpha[1]", "beta[1]", "kappa[1]", "phi", "eta[1]"))

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
  filter(!grepl('alpha_raw', parameter)) %>%
  mutate(index = parse_number(parameter)) %>%
  mutate(true_alpha = alpha_true[index])

plot_alpa <- ggplot(data=summary_alpha) +
  geom_point(aes(x=index, y=mean, color="estimated")) + 
  geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
  geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha=0.5) +
  geom_point(aes(x=index, y=true_alpha, color="true value")) +
  ggtitle("Alpha")
plot_alpa  

# look at beta
summary_beta <- fir_summary_df %>%
  rownames_to_column("parameter") %>%
  filter(grepl('beta', parameter)) %>%
  filter(!grepl('sigma_beta', parameter)) %>%
  mutate(index = parse_number(parameter)) %>%
  mutate(true_beta = beta_true[index])

plot_beta <- ggplot(data=summary_beta) +
  geom_point(aes(x=index, y=mean, color="estimated")) + 
  geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
  geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha = 0.5) +
  geom_point(aes(x=index, y=true_beta, color="true value")) +
  ggtitle("Beta -STAN")
plot_beta

# look at kappa
summary_kappa <- fir_summary_df %>%
  rownames_to_column("parameter") %>%
  filter(grepl('kappa', parameter)) %>%
  filter(!grepl('sigma_kappa', parameter)) %>%
  mutate(index = parse_number(parameter)) %>%
  mutate(true_kappa = kappa_true[index]) %>%
  mutate(kappa_drifted = true_kappa + phi_true*index)

plot_kappa <- ggplot(data=summary_kappa) +
  geom_point(aes(x=index, y=mean, color="estimated")) + 
  geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
  geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha=0.5) +
  #geom_point(aes(x=index, y=true_kappa, color="true value")) +
  geom_point(aes(x=index, y=kappa_drifted, color="true with drift"))+
  ggtitle("Kappa; time effect with drift - STAN")
plot_kappa

# look at eta
summary_eta <- fir_summary_df %>%
  rownames_to_column("parameter") %>%
  filter(grepl('eta', parameter)) %>%
  filter(!grepl('beta', parameter)) %>%
  mutate(index = parse_number(parameter)) %>%
  mutate(true_eta = data$eta)

plot_eta <- ggplot(data=summary_eta) +
  geom_point(aes(x=index, y=mean, color="estimated")) + 
  geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
  geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha=0.5) +
  geom_point(aes(x=index, y=true_eta, color="true value")) +
  ggtitle("Eta - STAN")
plot_eta


#   ----   estimated hyperparameters by inlabru and STAN compared   ----

p.prec.alpha <- ggplot(data.frame(res.inlabru$marginals.hyperpar) %>%
                         filter(Precision.for.alpha.x < 35)) + 
  geom_area(aes(x = Precision.for.alpha.x, y = Precision.for.alpha.y),fill = palette.basis[1], alpha = 0.4) + 
  geom_vline(data = res.inlabru$summary.hyperpar, aes(xintercept = mean[1]), color = palette.basis[1]) + 
  geom_vline(aes(xintercept = 1/sqrt(fir_summary_df$mean[1])), color=palette.basis[2]) + 
  labs(x = "Value of precision", y = " ", title = "Precision for alpha - comparison")
p.prec.alpha
# two values above 200

p.prec.beta <- ggplot(data.frame(res.inlabru$marginals.hyperpar) %>%
                        filter(Precision.for.beta.x < 200)) + 
  geom_area(aes(x = Precision.for.beta.x, y = Precision.for.beta.y, fill = "Inlabru"), alpha = 0.4) + 
  geom_vline(data = res.inlabru$summary.hyperpar, aes(xintercept = mean[2], color = "Inlabru", fill = "Inlabru")) + 
  geom_vline(aes(xintercept = tau.beta, color = "True value", fill = "True value")) + 
  geom_vline(aes(xintercept = 1/sqrt(fir_summary_df$mean[2]), color="STAN")) + 
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  labs(x = "Value of precision", y = " ", title = "Precision for beta - comparsion")
p.prec.beta
# two values above 1000

p.prec.kappa <- ggplot(data.frame(res.inlabru$marginals.hyperpar) %>%
                         filter(Precision.for.kappa.x < 175)) + 
  geom_area(aes(x = Precision.for.kappa.x, y = Precision.for.kappa.y), alpha = 0.4, fill = palette.basis[1]) + 
  geom_vline(data = res.inlabru$summary.hyperpar, aes(xintercept = mean[3]), color = palette.basis[1]) + 
  geom_vline(aes(xintercept = 1/sqrt(fir_summary_df$mean[3])), color=palette.basis[2]) + 
  labs(x = "Value of precision", y = " ", title = "Precision for kappa - comparison")
p.prec.kappa
# two values above 1000

p.prec.epsilon <- ggplot(data.frame(res.inlabru$marginals.hyperpar) %>%
                           filter(Precision.for.epsilon.x < 150000)) + 
  geom_area(aes(x = Precision.for.epsilon.x, y = Precision.for.epsilon.y, fill = "Inlabru"), alpha = 0.4) + 
  geom_vline(data = res.inlabru$summary.hyperpar, aes(xintercept = mean[4], color = "Inlabru", fill = "Inlabru")) + 
  geom_vline(aes(xintercept = 1/sqrt(fir_summary_df$mean[4]), color="STAN")) +
  geom_vline(aes(xintercept = tau.epsilon, color = "True value", fill = "True value")) + 
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  labs(x = "Value of precision", y = " ", title = "Precision for epsilon - comparison")
p.prec.epsilon

#   ----   plot results together   ---- 

p.int.c <- ggplot(data.frame(res.inlabru$marginals.fixed)) + 
  geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
  geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Inlabru")) + 
  geom_vline(data = fir_summary_df, aes(xintercept = mean[6], color = "STAN")) +
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  labs(x = "Value of intercept", y = " ", title = "Intercept - comparison")
p.int.c

p.alpha.c <- ggplot(data = data.alpha, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
  geom_point(aes(y = mean, color = "Inlabru", fill = "Inlabru")) + 
  geom_point(data = summary_alpha, aes(x = index, y = mean, color = "STAN", fill = "STAN")) +
  geom_line(data = summary_alpha, aes(x = index, y = `2.5%`, color = "STAN", fill = "STAN"), alpha = 0.6) +
  geom_line(data = summary_alpha, aes(x = index, y = `97.5%`, color = "STAN", fill = "STAN"), alpha = 0.6) +
  geom_point(aes(y = alpha.true, color = "True value", fill = "True value")) + 
  scale_color_manual(name = "",
                     values = palette.basis ) +
  scale_fill_manual(name = "",
                    values = palette.basis ) +
  labs(title="Alpha - comparison", x = "x", y='')

p.alpha.c

p.beta.c <- ggplot(data = data.beta, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
  geom_point(data = summary_beta, aes(x = index, y = mean, color = "STAN", fill = "STAN")) +
  geom_line(data = summary_beta, aes(x = index, y = `2.5%`, color = "STAN", fill = "STAN"), alpha = 0.6) +
  geom_line(data = summary_beta, aes(x = index, y = `97.5%`, color = "STAN", fill = "STAN"), alpha = 0.6) +
  geom_point(aes(y = beta.true, color = "True value", fill = "True value")) + 
  geom_point(aes(y = mean, color = "Inlabru", fill = "Inlabru")) + 
  scale_color_manual(name = "",
                     values = palette.basis ) +
  scale_fill_manual(name = "",
                    values = palette.basis ) +
  labs(x = "x", y = "beta", title = "Beta - comparison")
p.beta.c

p.kappa.c <- ggplot(data = data.kappa, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
  geom_point(aes(y = kappa.true, color = "True undrifted", fill = "True undrifted")) + 
  geom_point(aes(y = mean, color = "Inlabru", fill = "Inlabru")) + 
  geom_point(data = summary_kappa, aes(x = index, y = mean, color = "STAN", fill = "STAN")) +
  geom_line(data = summary_kappa, aes(x = index, y = `2.5%`, color = "STAN", fill = "STAN"), alpha = 0.6) +
  geom_line(data = summary_kappa, aes(x = index, y = `97.5%`, color = "STAN", fill = "STAN"), alpha = 0.6) +
  geom_point(data = summary_kappa, aes(x = index, y = kappa_drifted, color = "True drifted", fill = "True drifted")) +
  scale_color_manual(name = "",
                     values = palette.basis ) +
  scale_fill_manual(name = "",
                    values = palette.basis ) +
  labs(x = "t", y = "kappa", title = "Kappa - comparison")

p.kappa.c

p.phi.c <- ggplot(data.frame(res.inlabru$marginals.fixed)) + 
  geom_area(aes(x = phi.x, y = phi.y, fill = "Inlabru"), alpha = 0.4) + 
  geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[2], color = "Inlabru", fill="Inlabru")) + 
  geom_vline(aes(xintercept = phi_true, color="True", fill="True")) +
  geom_vline(data=fir_summary_df, aes(xintercept=mean[5], color="STAN", fill="STAN")) + 
  geom_vline(data=fir_summary_df, aes(xintercept=`2.5%`[5], color="STAN", fill="STAN"), alpha=0.5) +
  geom_vline(data=fir_summary_df, aes(xintercept=`97.5%`[5], color="STAN", fill="STAN"), alpha=0.5) +
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) +
  labs(x = "Value of phi", y = " ", title = "Phi - comparison")
p.phi.c

p.eta.c <- ggplot(data = data.eta) +
  geom_point(aes(x=data$xt, y = eta.sim, color="Inlabru")) +
  geom_point(data=data, aes(x=xt, y = eta, color="True")) +
  geom_point(data=summary_eta, aes(x=index, y=mean, color="STAN")) +
  scale_color_manual(name = " ", values = palette.basis) + 
  labs(x=" ", y="Eta", title="Eta- comparison")
p.eta.c


