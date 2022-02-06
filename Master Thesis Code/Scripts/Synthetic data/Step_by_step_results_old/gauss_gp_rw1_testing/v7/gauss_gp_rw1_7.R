#' This program investigates inference with:
#' - The traditional_lc_log_prec version - Gaussian lc-model. 
#' All random effects are modelled as iid, in inlabru and in stan
#' Constraints are imposed as usual in inlabru, and with soft constraints in stan
#' They hyperparameters are fixed. 
#' 

#   ----   Load libraries and set workspace   ----
library("tidyverse")
library("inlabru")
library("ggplot2")
library("INLA")
library("patchwork")
library("rstan")

setwd("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master\ Thesis\ Code")

investigation.name <- "gauss_gp_rw1_testing"
investigation.path <- file.path(investigation.name, "v7")

#   ----    Retrieve the data   ----

synthetic.male.lung.v7 <- function(){
  obs <- read.csv("Data/synthetic_male_lung_7.csv")
  obs <- obs %>% mutate(x.old = x, x = x - 9, x.c = x) %>%
    select(-X) %>%
    mutate(x = x + 1, t = t+1, x.c = x.c + 1)
  
  obs.trad <- obs %>% 
    select(c(x, t, xt, age.int, year, x.c, alpha, beta, kappa, intercept, epsilon,
             eta, tau.alpha, tau.beta, tau.kappa, tau.epsilon, E)) %>%
    #mutate(eta = eta) %>%
    mutate(eta.no.error = intercept + alpha + beta*kappa) %>%
    mutate(mr_gaussian = exp(eta)) %>%
    mutate(Y_gaussian  = mr_gaussian * E)
  
  underlying.effects <- list(obs = obs.trad, nx = 18, nt = 18,
                             alpha.true = {obs %>% filter(t == 0)}$alpha,
                             beta.true = {obs %>% filter(t == 0)}$beta,
                             kappa.true = {obs %>% filter(x == 0)}$kappa,
                             intercept = unique(obs$intercept),
                             age.intercept.true = unique(obs$intercept),
                             tau.alpha.true = unique(obs$tau.alpha),
                             tau.beta.true = unique(obs$tau.beta),
                             tau.kappa.true = unique(obs$tau.kappa),
                             tau.epsilon.true = unique(obs$tau.epsilon))
  
  return(list(obs = obs, underlying.effects = underlying.effects))
}

# We use this data for both inlabru and stan
config.data <- synthetic.male.lung.v7()
obs <- config.data$obs
underlying.effects <- config.data$underlying.effects

#   ----   Run STAN analysis   ----
# Running traditional lc version of stan, 
# implemented with log-precisions, random effects as iid, and no constraints

stan.output  <- file.path("Scripts/Synthetic data/Step_by_step_results", investigation.path)
source("Scripts/Synthetic\ data/run_stan_functions.R")

run_stan <- function(stan_program, obs, chains, warmup, iter, output.path, config.name, markov=TRUE){
  
  input_stan <- list(
    X=length(unique(obs$x)),
    T=length(unique(obs$t)),
    x=(obs$x),
    t=(obs$t),
    log_mr = obs$eta,
    nx = length(unique(obs$x)),
    nt = length(unique(obs$t))
  )
  
  print(input_stan)
  
  stan_fit <- stan(
    file = stan_program,
    data = input_stan,
    chains = chains,
    iter = iter,
    warmup = warmup,
    refresh = iter%/%10,
    seed = 123
  )
  
  store_stan_results_traditional(
    fit=stan_fit, output.path=output.path, config=config.name,
    chains=chains, warmup=warmup, iter=iter, stan_program=stan_program)
  
  return(stan_fit)
}

stan_fit <- run_stan(
  stan_program=file.path(stan.output, "stan_gauss_gp_rw1.stan"),
  obs = obs, chains=4, warmup = 400, iter = 4000, output.path = stan.output,
  config.name = investigation.name, markov=F)

inlabru.gaus.gp.rw1 <- function(obs, max_iter=30){
  #'Implements inlabru analysis for lc model, fixing the precisions and modelling all random effects as iid
  #'
  #'@param obs: Contains the observed data and the real underlying random effects
  #'@param max_iter (int): maximum number of iterations in inlabru
  
  nx = length(unique(obs$x))
  nt = length(unique(obs$t))
  
  # constraints for the age effect beta
  A.beta = matrix(1, nrow = 1, ncol = nx)
  e.beta = 1
  
  # fixed.theta.alpha <- list(prec = list(initial = log(1.96), fixed = T))
  # fixed.theta.beta <- list(prec = list(initial = log(64), fixed = T))
  # fixed.theta.kappa <- list(prec = list(initial = log(336), fixed = T))
  # fixed.theta.epsilon <- list(prec = list(initial = log(420), fixed = T))
  
  loggamma.prior <- list(prec = list(prior = 'loggamma', param = c(1,0.00005), initial = log(1)))
  loggamma.prior.high.variance <- list(prec = list(prior = 'loggamma', param = c(1,0.005), initial = log(1)))
  
  comp = ~ -1 +
    Int(1, prec.linear = 0.001, mean.linear = 0) +
    alpha(x, model = "rw1", hyper = loggamma.prior, constr = TRUE) +
    beta(x.c, model = "iid", hyper = loggamma.prior, extraconstr = list(A = A.beta, e = e.beta)) +
    kappa(t, model = "rw1", hyper = loggamma.prior.high.variance, constr = TRUE)
  
  formula = eta ~ Int + alpha + beta*kappa
  
  likelihood = like(formula = formula, family = "gaussian", data = obs)
  
  c.compute <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE, return.marginals.predictor = TRUE)  # control.compute
  c.family <- list(hyper = loggamma.prior)
  
  res.inlabru = bru(components = comp,
                    likelihood, 
                    options = list(verbose = F,
                                   bru_verbose = 4, 
                                   num.threads = "1:1",
                                   control.compute = c.compute,
                                   bru_max_iter=max_iter,
                                   #control.predictor = list(compute = TRUE),
                                   control.family  = c.family
                    ))
  return(res.inlabru)
}

res.inlabru <- inlabru.gaus.gp.rw1(obs, max_iter = 100)

source("Scripts/Functions/plotters.R")
source("Scripts/Synthetic data/plot_inlabru_vs_underlying.R")
source("Scripts/Synthetic data/plot_inlabru_stan_compared.R")
source("Scripts/Synthetic data/plot_stan_vs_underlying.R")

output.path <- stan.output

plots.summaries.inlabru <- plot.inlabru.vs.underlying.traditional.lc(
  res.inlabru,
  underlying.effects,
  path.to.storage = output.path,
  save=F)

load(file.path(stan.output, paste("stan_", investigation.name, ".Rda", sep = "")))

load(file=file.path(stan.output, "draws_intercept.RData"))
load(file=file.path(stan.output, "draws_tau_epsilon.RData"))
load(file.path(stan.output, "draws_tau_alpha.RData"))
load(file.path(stan.output, "draws_tau_beta.RData"))
load(file.path(stan.output, "draws_tau_kappa.RData"))
load(file.path(stan.output, "draws_alpha.RData"))
load(file.path(stan.output, "draws_beta.RData"))
load(file.path(stan.output, "draws_kappa.RData"))
load(file.path(stan.output, "draws_eta_100.RData"))
load(file.path(stan.output, "draws_eta.RData"))
load(file.path(stan.output, "draws_eta_reduced.RData"))

stan.marginals <- list(intercept_draws = intercept_draws,
                       tau_epsilon_draws = tau_epsilon_draws,
                       tau_alpha_draws = tau_alpha_draws,
                       tau_beta_draws = tau_beta_draws,
                       tau_kappa_draws = tau_kappa_draws,
                       alpha_draws = alpha_draws,
                       beta_draws = beta_draws,
                       kappa_draws = kappa_draws,
                       eta_draws = eta_draws)

stan.res <- produce.stan.plots(stan_df=stan_lc_df,
                               underlying.effects=underlying.effects,
                               plot.func=plot.stan.vs.underlying.synthetic.cancer,
                               save.func=function(...) {save.stan.plots.lc.rw2(..., save=F)},
                               path.to.storage=output.path,
                               summaries.func=produce.summaries.stan.traditional)

plots_compared <- produce.compared.plots(
  stan.summaries = stan.res$summaries,
  stan.marginals = stan.marginals,
  inlabru.summaries = plots.summaries.inlabru$summaries,
  res.inlabru = res.inlabru,
  underlying.effects = underlying.effects,
  plot.func = function(...) {plot.inlabru.stan.traditional.lc(..., cohort=FALSE, tau.beta.cutoff = 400, tau.kappa.cutoff = 150, tau.alpha.cutoff = 100, a45=F)},
  #plot.func = function(...) {plot.inlabru.stan.traditional.lc.fixed.hypers(..., cohort=FALSE, tau.beta.cutoff = 700, tau.kappa.cutoff = 500, tau.alpha.cutoff = 10, a45=F)},
  #plot.func = function(...) {plot.inlabru.stan.traditional.lc.no.beta(..., cohort=FALSE, tau.beta.cutoff = 5000, tau.kappa.cutoff = 5000, tau.alpha.cutoff = 100, a45=F)},
  #plot.func = function(...) {plot.inlabru.stan.traditional.lc.fixed.hypers.no.beta(..., cohort=FALSE, tau.beta.cutoff = 5000, tau.kappa.cutoff = 5000, tau.alpha.cutoff = 100, a45=F)},
  save.func = function(...) {save.compared.rw2(..., cohort=FALSE)},
  path.to.storage=output.path)

#   ----   Sample predictor   ----

stan.predictor.df <- data.frame(eta_draws)
plot.predictor.inlabru.stan.compared(res.inlabru, stan.predictor.df, path.to.storage = output.path, a45=T)

p.predictor <- ggplot() + 
  geom_area(data = data.frame(res.inlabru$marginals.linear.predictor$APredictor.054), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(data = stan.predictor.df, aes(x = X54, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins = 100) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Predictor[54]", x = "", y = "")
p.predictor

ggsave("predictor_54.pdf", p.predictor, path=output.path, device = "pdf", dpi="retina", height = 5, width = 8)

p.alpha <- ggplot() + 
  geom_area(data = data.frame(res.inlabru$marginals.random$alpha$index.9), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_histogram(data = stan.alpha.df, aes(x = X9, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins = 100) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Alpha[9]", x = "", y = "")
p.alpha

#   ----   Plot marginals of random effects   ----

stan.beta.df <- data.frame(beta_draws)
plot.beta.inlabru.stan.compared(res.inlabru, stan.beta.df, path.to.storage = output.path, a45=T)

stan.kappa.df <- data.frame(kappa_draws)
plot.kappa.inlabru.stan.compared(res.inlabru, stan.kappa.df, path.to.storage = output.path)

list_of_draws <- rstan::extract(stan_fit)

stan.alpha.df <- data.frame(alpha_draws)
plot.alpha.inlabru.stan.compared(res.inlabru, stan.alpha.df, path.to.storage = output.path, a45=T)

p.intercept <- ggplot() + 
  geom_area(data = data.frame(res.inlabru$marginals.fixed$Int), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  geom_density(data = data.frame("x" = list_of_draws$intercept), aes(x = x, color = "Stan", fill = "Stan"), alpha = 0.2) +
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Intercept", x = "", y = "")

ggsave("intercept.pdf", p.intercept, path=output.path, device = "pdf", dpi="retina", height = 5, width = 8)

#   ----   Plot densities of hyperparameters   ----

p.tau.alpha <- ggplot() + 
  geom_density(data=data.frame("x" = list_of_draws$tau_alpha), aes(x = x, color = "Stan", fill = "Stan"), alpha = 0.2) + 
  geom_area(data = data.frame(res.inlabru$marginals.hyperpar$`Precision for alpha`) %>% filter(x < 70), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Tau alpha", x = "", y = "")

p.theta.alpha <- ggplot() + 
  geom_density(data=data.frame("x" = list_of_draws$theta_alpha), aes(x = x, color = "Stan", fill = "Stan"), alpha = 0.2) + 
  geom_area(data = data.frame(res.inlabru$internal.marginals.hyperpar$`Log precision for alpha`), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Theta alpha", x = "", y = "")

p.tau.kappa <- ggplot() + 
  geom_density(data=data.frame("x" = list_of_draws$tau_kappa), aes(x = x, color = "Stan", fill = "Stan"), alpha = 0.2) + 
  geom_area(data = data.frame(res.inlabru$marginals.hyperpar$`Precision for kappa`) %>% filter(x < 3000), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Tau kappa", x = "", y = "")

p.theta.kappa <- ggplot() + 
  geom_density(data=data.frame("x" = list_of_draws$theta_kappa), aes(x = x, color = "Stan", fill = "Stan"), alpha = 0.2) + 
  geom_area(data = data.frame(res.inlabru$internal.marginals.hyperpar$`Log precision for kappa`), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Theta kappa", x = "", y = "")

p.tau.eta <- ggplot() + 
  geom_density(data=data.frame("x" = list_of_draws$tau_epsilon), aes(x = x, color = "Stan", fill = "Stan"), alpha = 0.2) + 
  geom_area(data = data.frame(res.inlabru$marginals.hyperpar$`Precision for the Gaussian observations`) %>% filter(x < 350), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Tau eta", x = "", y = "")

p.theta.eta <- ggplot() + 
  geom_density(data=data.frame("x" = list_of_draws$theta_epsilon), aes(x = x, color = "Stan", fill = "Stan"), alpha = 0.2) + 
  geom_area(data = data.frame(res.inlabru$internal.marginals.hyperpar$`Log precision for the Gaussian observations`), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Theta eta", x = "", y = "")

p.tau.beta <- ggplot() + 
  geom_density(data=data.frame("x" = list_of_draws$tau_beta), aes(x = x, color = "Stan", fill = "Stan"), alpha = 0.2) + 
  geom_area(data = data.frame(res.inlabru$marginals.hyperpar$`Precision for beta`) %>% filter(x < 350), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Tau beta", x = "", y = "")

p.theta.beta <- ggplot() + 
  geom_density(data=data.frame("x" = list_of_draws$theta_beta), aes(x = x, color = "Stan", fill = "Stan"), alpha = 0.2) + 
  geom_area(data = data.frame(res.inlabru$internal.marginals.hyperpar$`Log precision for beta`), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Theta beta", x = "", y = "")

p.tau <- (p.tau.alpha | p.tau.beta)/(p.tau.kappa |p.tau.eta) + plot_layout(guides = "collect")
p.theta <- (p.theta.alpha | p.theta.beta)/(p.theta.kappa | p.theta.eta) + plot_layout(guides = "collect")

ggsave("marginals_tau.pdf", p.tau, path=output.path, device = "pdf", dpi="retina", height = 5, width = 8)
ggsave("marginals_theta.pdf", p.theta, path=output.path, device = "pdf", dpi="retina", height = 5, width = 8)

#   ----   Plot histograns of hyperparameters   ----

p.tau.alpha.h <- ggplot() + 
  geom_histogram(data=data.frame("x" = list_of_draws$tau_alpha), aes(x = x, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.3, bins = 100, size=0.2) + 
  geom_area(data = data.frame(res.inlabru$marginals.hyperpar$`Precision for alpha`) %>% filter(x < 70), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Tau alpha", x = "", y = "")

p.theta.alpha.h <- ggplot() + 
  geom_histogram(data=data.frame("x" = list_of_draws$theta_alpha), aes(x = x, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.3, bins = 100, size=0.2) + 
  geom_area(data = data.frame(res.inlabru$internal.marginals.hyperpar$`Log precision for alpha`), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Theta alpha", x = "", y = "")

p.tau.kappa.h <- ggplot() + 
  geom_histogram(data=data.frame("x" = list_of_draws$tau_kappa), aes(x = x, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.3, bins = 100, size=0.2) + 
  geom_area(data = data.frame(res.inlabru$marginals.hyperpar$`Precision for kappa`) %>% filter(x < 3000), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Tau kappa", x = "", y = "")

p.theta.kappa.h <- ggplot() + 
  geom_histogram(data=data.frame("x" = list_of_draws$theta_kappa), aes(x = x, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.3, bins = 100, size=0.2) + 
  geom_area(data = data.frame(res.inlabru$internal.marginals.hyperpar$`Log precision for kappa`), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Theta kappa", x = "", y = "")

p.tau.eta.h <- ggplot() + 
  geom_histogram(data=data.frame("x" = list_of_draws$tau_epsilon), aes(x = x, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.3, bins = 100, size=0.2) + 
  geom_area(data = data.frame(res.inlabru$marginals.hyperpar$`Precision for the Gaussian observations`) %>% filter(x < 350), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Tau eta", x = "", y = "")

p.theta.eta.h <- ggplot() + 
  geom_histogram(data=data.frame("x" = list_of_draws$theta_epsilon), aes(x = x, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.3, bins = 100, size=0.2) + 
  geom_area(data = data.frame(res.inlabru$internal.marginals.hyperpar$`Log precision for the Gaussian observations`), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Theta eta", x = "", y = "")

p.tau.beta.h <- ggplot() + 
  geom_histogram(data=data.frame("x" = list_of_draws$tau_beta), aes(x = x, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.3, bins = 100, size=0.2) + 
  geom_area(data = data.frame(res.inlabru$marginals.hyperpar$`Precision for beta`) %>% filter(x < 350), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Tau beta", x = "", y = "")

p.theta.beta.h <- ggplot() + 
  geom_histogram(data=data.frame("x" = list_of_draws$theta_beta), aes(x = x, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.3, bins = 100, size=0.2) + 
  geom_area(data = data.frame(res.inlabru$internal.marginals.hyperpar$`Log precision for beta`), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Theta beta", x = "", y = "")

p.tau.h <- (p.tau.alpha.h | p.tau.beta.h)/(p.tau.kappa.h |p.tau.eta.h) + plot_layout(guides = "collect")
p.theta.h <- (p.theta.alpha.h | p.theta.beta.h)/(p.theta.kappa.h | p.theta.eta.h) + plot_layout(guides = "collect")

ggsave("marginals_tau_histogram.pdf", p.tau.h, path=output.path, device = "pdf", dpi="retina", height = 5, width = 8)
ggsave("marginals_theta_histogram.pdf", p.theta.h, path=output.path, device = "pdf", dpi="retina", height = 5, width = 8)



