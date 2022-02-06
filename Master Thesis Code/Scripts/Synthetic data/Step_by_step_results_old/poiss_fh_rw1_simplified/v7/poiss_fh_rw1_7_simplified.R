#' This program investigates inference with:
#' - The traditional_lc_log_prec version - Gaussian lc-model. 
#' All random effects are modelled as iid, in inlabru and in stan
#' Constraints are imposed as usual in inlabru, and with soft constraints in stan
#' They hyperparameters are fixed. 
#' 
#' 
#' NOTE: For this exact simplified axample, the model should experience overdispersion
#' This is because we have created the data with an error term that we do not include
#' in the simpler version of the model. 

#   ----   Load libraries and set workspace   ----
#   ----   Load libraries and set workspace   ----
set_workdirectory <- function(markov=TRUE){
  if(markov){
    .libPaths("~/Documents/R_libraries")
    setwd("~/Documents/GitHub/Masteroppgave/Masters-thesis/Master Thesis Code")
  } else {
    setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code")
  }
}

set_workdirectory(markov=F)

library("tidyverse")
library("inlabru")
library("ggplot2")
library("INLA")
library("patchwork")
library("rstan")

#setwd("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master\ Thesis\ Code")

#investigation.name <- "poiss_fh_rw1_simplified"
output.path <- "Scripts/Synthetic data/Step_by_step_results/poiss_fh_rw1_simplified/v7"

#   ----    Retrieve the data   ----

synthetic.male.lung.v7 <- function(){
  obs <- read.csv("Data/synthetic_male_lung_7.csv")
  obs <- obs %>% mutate(x.old = x, x = x - 9, x.c = x) %>%
    select(-X) %>%
    mutate(x = x + 1, t = t + 1, x.c = x.c + 1)
  
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

stan.output  <- output.path
#source("Scripts/Synthetic\ data/run_stan_functions.R")

run_stan <- function(stan_program, obs, chains, warmup, iter, output.path, config.name, markov=TRUE){
  
  input_stan <- list(
    X=length(unique(obs$x)),
    T=length(unique(obs$t)),
    x=(obs$x),
    t=(obs$t),
    Y = obs$Y,
    nx = length(unique(obs$x)),
    nt = length(unique(obs$t)),
    E = obs$E
  )
  
  stan_fit <- stan(
    file = stan_program,
    data = input_stan,
    chains = chains,
    iter = iter,
    warmup = warmup,
    refresh = iter%/%10,
    seed = 123
  )
  
  # store_stan_results(
  #   fit=stan_fit, output.path=output.path, config=config.name,
  #   chains=chains, warmup=warmup, iter=iter, stan_program=stan_program)
  
  return(stan_fit)
}

stan_fit <- run_stan(
  stan_program=file.path(stan.output, "stan_poiss_fh_rw1_7_simplified.stan"),
  obs = obs, chains=3, warmup = 400, iter = 4000, output.path = stan.output,
  config.name = investigation.name, markov=F)

stan_fit_epsilon <- run_stan(
  stan_program=file.path(stan.output, "stan_poiss_fh_rw1_7.stan"),
  obs = obs, chains=3, warmup = 400, iter = 4000, output.path = stan.output,
  config.name = investigation.name, markov=F)

summary_stan <- data.frame(summary(stan_fit)) %>% select(1:10) %>% rename_with( ~ gsub("summary.", "", .x))
save(summary_stan, file = file.path(output.path, "summary_stan_simplified.Rda"))
list_of_draws <- rstan::extract(stan_fit)
save(list_of_draws, file = file.path(output.path, "draws_simplified.RData"))

summary_stan_epsilon <- data.frame(summary(stan_fit_epsilon)) %>% select(1:10) %>% rename_with( ~ gsub("summary.", "", .x))
save(summary_stan_epsilon, file = file.path(output.path, "summary_stan_epsilon.Rda"))
list_of_draws_epsilon <- rstan::extract(stan_fit_epsilon)
save(list_of_draws_epsilon, file = file.path(output.path, "draws_epsilon.RData"))

inlabru.pois.fh.rw1.simplified <- function(obs, max_iter=30){
  #'Implements inlabru analysis for lc model, fixing the precisions and modelling all random effects as iid
  #'
  #'@param obs: Contains the observed data and the real underlying random effects
  #'@param max_iter (int): maximum number of iterations in inlabru
  
  nx = length(unique(obs$x))
  nt = length(unique(obs$t))
  
  # constraints for the age effect beta
  A.beta = matrix(1, nrow = 1, ncol = nx)
  e.beta = 1
  
  fixed.theta.alpha <- list(prec = list(initial = log(1.96), fixed = T))
  fixed.theta.beta <- list(prec = list(initial = log(64), fixed = T))
  fixed.theta.kappa <- list(prec = list(initial = log(336), fixed = T))
  fixed.theta.epsilon <- list(prec = list(initial = log(420), fixed = T))
  
  comp = ~ -1 +
    alpha(x, model = "rw1", hyper = fixed.theta.alpha, constr = FALSE) +
    beta(x.c, model = "iid", hyper = fixed.theta.beta, extraconstr = list(A = A.beta, e = e.beta)) +
    kappa(t, model = "rw1", hyper = fixed.theta.kappa, constr = TRUE)
    #epsilon(xt, model = "iid", hyper = fixed.theta.epsilon, constr = FALSE)
  
  #formula = Y ~ Int + alpha + beta*kappa + epsilon
  formula = Y ~ alpha + beta*kappa
  
  likelihood = like(formula = formula, family = "poisson", data = obs, E = obs$E)
  
  c.compute <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE, return.marginals.predictor = TRUE)  # control.compute
  
  res.inlabru = bru(components = comp,
                    likelihood, 
                    options = list(verbose = F,
                                   bru_verbose = 4, 
                                   num.threads = "1:1",
                                   control.compute = c.compute,
                                   bru_max_iter=max_iter,
                                   control.predictor = list(link = 1)
                    ))
  
  return(res.inlabru)
}

inlabru.pois.fh.rw1.epsilon <- function(obs, max_iter=30){
  #'Implements inlabru analysis for lc model, fixing the precisions and modelling all random effects as iid
  #'
  #'@param obs: Contains the observed data and the real underlying random effects
  #'@param max_iter (int): maximum number of iterations in inlabru
  
  nx = length(unique(obs$x))
  nt = length(unique(obs$t))
  
  # constraints for the age effect beta
  A.beta = matrix(1, nrow = 1, ncol = nx)
  e.beta = 1
  
  fixed.theta.alpha <- list(prec = list(initial = log(1.96), fixed = T))
  fixed.theta.beta <- list(prec = list(initial = log(64), fixed = T))
  fixed.theta.kappa <- list(prec = list(initial = log(336), fixed = T))
  fixed.theta.epsilon <- list(prec = list(initial = log(420), fixed = T))
  
  comp = ~ -1 +
    alpha(x, model = "rw1", hyper = fixed.theta.alpha, constr = FALSE) +
    beta(x.c, model = "iid", hyper = fixed.theta.beta, extraconstr = list(A = A.beta, e = e.beta)) +
    kappa(t, model = "rw1", hyper = fixed.theta.kappa, constr = TRUE) + 
    epsilon(xt, model = "iid", hyper = fixed.theta.epsilon, constr = FALSE)
  
  formula = Y ~ alpha + beta*kappa + epsilon
  
  likelihood = like(formula = formula, family = "poisson", data = obs, E = obs$E)
  
  c.compute <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE, return.marginals.predictor = TRUE)  # control.compute
  
  res.inlabru = bru(components = comp,
                    likelihood, 
                    options = list(verbose = F,
                                   bru_verbose = 4, 
                                   num.threads = "1:1",
                                   control.compute = c.compute,
                                   bru_max_iter=max_iter,
                                   control.predictor = list(link = 1)
                    ))
  
  return(res.inlabru)
}

res.inlabru.simplified <- inlabru.pois.fh.rw1.simplified(obs, max_iter = 100)
res.inlabru.epsilon <- inlabru.pois.fh.rw1.epsilon(obs, max_iter = 100)

save(res.inlabru.simplified, file = file.path(output.path, "res_inlabru.RData"))
save(res.inlabru.epsilon, file = file.path(output.path, "res_inlabru.RData"))

source("Scripts/Functions/plotters.R")
#source("Scripts/Synthetic data/plot_inlabru_vs_underlying.R")
#source("Scripts/Synthetic data/plot_inlabru_stan_compared.R")
#source("Scripts/Synthetic data/plot_stan_vs_underlying.R")

# plots.summaries.inlabru <- plot.inlabru.vs.underlying.synthetic.cancer.fixed.effects(
#   res.inlabru,
#   underlying.effects,
#   path.to.storage = output.path,
#   save=F)
# 
# load(file.path(stan.output, paste("stan_", investigation.name, ".Rda", sep = "")))
# 
# load(file=file.path(stan.output, "draws_intercept.RData"))
# load(file=file.path(stan.output, "draws_tau_epsilon.RData"))
# load(file.path(stan.output, "draws_tau_alpha.RData"))
# load(file.path(stan.output, "draws_tau_beta.RData"))
# load(file.path(stan.output, "draws_tau_kappa.RData"))
# load(file.path(stan.output, "draws_alpha.RData"))
# load(file.path(stan.output, "draws_beta.RData"))
# load(file.path(stan.output, "draws_kappa.RData"))
# load(file.path(stan.output, "draws_eta_100.RData"))
# load(file.path(stan.output, "draws_eta.RData"))
# load(file.path(stan.output, "draws_eta_reduced.RData"))
# 
# stan.marginals <- list(intercept_draws = intercept_draws,
#                        tau_epsilon_draws = tau_epsilon_draws,
#                        tau_alpha_draws = tau_alpha_draws,
#                        tau_beta_draws = tau_beta_draws,
#                        tau_kappa_draws = tau_kappa_draws,
#                        alpha_draws = alpha_draws,
#                        beta_draws = beta_draws,
#                        kappa_draws = kappa_draws,
#                        eta_draws = eta_draws)
# 
stan.res <- produce.stan.plots(stan_df=stan_lc_df,
                               underlying.effects=underlying.effects,
                               plot.func=plot.stan.vs.underlying.synthetic.cancer,
                               save.func=function(...) {save.stan.plots.lc.rw2(..., save=F)},
                               path.to.storage=output.path,
                               summaries.func=produce.summaries.stan.poiss.lc)
# 
# plots_compared <- produce.compared.plots(
#   stan.summaries = stan.res$summaries,
#   stan.marginals = stan.marginals,
#   inlabru.summaries = plots.summaries.inlabru$summaries,
#   res.inlabru = res.inlabru,
#   underlying.effects = underlying.effects,
#   plot.func = function(...) {plot.inlabru.stan.compared.poisson.lc.fixed.hypers(..., tau.beta.cutoff = 700, tau.kappa.cutoff = 500, tau.alpha.cutoff = 100)},
#   #plot.func = function(...) {plot.inlabru.stan.traditional.lc(..., cohort=FALSE, tau.beta.cutoff = 700, tau.kappa.cutoff = 500, tau.alpha.cutoff = 10, a45=F)},
#   #plot.func = function(...) {plot.inlabru.stan.traditional.lc.fixed.hypers(..., cohort=FALSE, tau.beta.cutoff = 700, tau.kappa.cutoff = 500, tau.alpha.cutoff = 10, a45=F)},
#   #plot.func = function(...) {plot.inlabru.stan.traditional.lc.no.beta(..., cohort=FALSE, tau.beta.cutoff = 5000, tau.kappa.cutoff = 5000, tau.alpha.cutoff = 100, a45=F)},
#   #plot.func = function(...) {plot.inlabru.stan.traditional.lc.fixed.hypers.no.beta(..., cohort=FALSE, tau.beta.cutoff = 5000, tau.kappa.cutoff = 5000, tau.alpha.cutoff = 100, a45=F)},
#   save.func = function(...) {save.compared.rw2(..., cohort=FALSE, png=F)},
#   path.to.storage=output.path)

#   ----   Plot resutls - simplified   ---
#   ----   Plot marginals of predictor   ----

stan.predictor.df.simplified <- data.frame(list_of_draws$eta)
plot.predictor.inlabru.stan.compared(res.inlabru.simplified, stan.predictor.df.simplified, path.to.storage = output.path, a45=T, filename = "predictor_simplified")

#   ----   Plot marginals of random effects   ----

stan.beta.df.simp <- data.frame(list_of_draws$beta)
plot.beta.inlabru.stan.compared(res.inlabru.simplified, stan.beta.df.simp, path.to.storage = output.path, a45=T, filename = "beta_simplified")

stan.kappa.df.simp <- data.frame(list_of_draws$kappa)
plot.kappa.inlabru.stan.compared(res.inlabru.simplified, stan.kappa.df.simp, path.to.storage = output.path, filename = "kappa_simplified")

stan.alpha.df.simp <- data.frame(list_of_draws$alpha)
plot.alpha.inlabru.stan.compared(res.inlabru.simplified, stan.alpha.df.simp, path.to.storage = output.path, filename = "alpha_simplified", a45=T)

#   ----   Quickly - plot comparison of inlabru estimates and true random effects   ----
p.alpha.simplified <- ggplot() + 
  geom_point(data = obs, aes(x = x, y = alpha + intercept, color = "True", fill = "True")) + 
  geom_ribbon(data = data.frame(res.inlabru.simplified$summary.random$alpha), aes(x = ID, ymin = X0.025quant, ymax = X0.975quant, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) +
  geom_point(data = data.frame(res.inlabru.simplified$summary.random$alpha), aes(x = ID, y = mean, fill = "Inlabru", color = "Inlabru")) + 
  theme_classic() + 
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title = "Alpha", x = "x", y = "")

p.beta.simplified <- ggplot() + 
  geom_point(data = obs, aes(x = x, y = beta, color = "True", fill = "True")) + 
  geom_errorbar(data = data.frame(res.inlabru.simplified$summary.random$beta), aes(x = ID, ymin = X0.025quant, ymax = X0.975quant, color = "Inlabru", fill = "Inlabru"), alpha = 0.7) +
  geom_point(data = data.frame(res.inlabru.simplified$summary.random$beta), aes(x = ID, y = mean, fill = "Inlabru", color = "Inlabru")) + 
  theme_classic() + 
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title = "Beta", x = "x", y = "")

p.kappa.simplified <- ggplot() + 
  geom_point(data = obs, aes(x = t, y = kappa, color = "True", fill = "True")) + 
  geom_ribbon(data = data.frame(res.inlabru.simplified$summary.random$kappa), aes(x = ID, ymin = X0.025quant, ymax = X0.975quant, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) +
  geom_point(data = data.frame(res.inlabru.simplified$summary.random$kappa), aes(x = ID, y = mean, fill = "Inlabru", color = "Inlabru")) + 
  theme_classic() + 
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title = "Kappa", x = "t", y = "")

p.random.simplified <- (p.alpha.simplified | p.beta.simplified)/(p.kappa.simplified) + plot_layout(guides = "collect")

ggsave("random_effects_simplified.pdf", p.random.simplified, path = output.path, dpi = "retina", width = 8, height = 5)

#   ----   Plot resutls - with epsilon  ---
#   ----   Plot marginals of predictor   ----

stan.predictor.df.eps<- data.frame(list_of_draws_epsilon$eta)
plot.predictor.inlabru.stan.compared(res.inlabru.epsilon, stan.predictor.df.eps, path.to.storage = output.path, a45=T, filename = "predictor_eps")

#   ----   Plot marginals of random effects   ----

stan.beta.df.eps <- data.frame(list_of_draws_epsilon$beta)
plot.beta.inlabru.stan.compared(res.inlabru.epsilon, stan.beta.df.eps, path.to.storage = output.path, a45=T, filename = "beta_eps")

stan.kappa.df.eps <- data.frame(list_of_draws_epsilon$kappa)
plot.kappa.inlabru.stan.compared(res.inlabru.epsilon, stan.kappa.df.eps, path.to.storage = output.path, filename = "kappa_eps")

stan.alpha.df.eps <- data.frame(list_of_draws_epsilon$alpha)
plot.alpha.inlabru.stan.compared(res.inlabru.epsilon, stan.alpha.df.eps, path.to.storage = output.path, filename = "alpha_eps", a45=T)

#   ----   Quickly - plot comparison of inlabru estimates and true random effects   ----
p.alpha.eps <- ggplot() + 
  geom_point(data = obs, aes(x = x, y = alpha + intercept, color = "True", fill = "True")) + 
  geom_ribbon(data = data.frame(res.inlabru.epsilon$summary.random$alpha), aes(x = ID, ymin = X0.025quant, ymax = X0.975quant, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) +
  geom_point(data = data.frame(res.inlabru.epsilon$summary.random$alpha), aes(x = ID, y = mean, fill = "Inlabru", color = "Inlabru")) + 
  theme_classic() + 
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title = "Alpha", x = "x", y = "")

p.beta.eps <- ggplot() + 
  geom_point(data = obs, aes(x = x, y = beta, color = "True", fill = "True")) + 
  geom_errorbar(data = data.frame(res.inlabru.epsilon$summary.random$beta), aes(x = ID, ymin = X0.025quant, ymax = X0.975quant, color = "Inlabru", fill = "Inlabru"), alpha = 0.7) +
  geom_point(data = data.frame(res.inlabru.epsilon$summary.random$beta), aes(x = ID, y = mean, fill = "Inlabru", color = "Inlabru")) + 
  theme_classic() + 
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title = "Beta", x = "x", y = "")

p.kappa.eps <- ggplot() + 
  geom_point(data = obs, aes(x = t, y = kappa, color = "True", fill = "True")) + 
  geom_ribbon(data = data.frame(res.inlabru.epsilon$summary.random$kappa), aes(x = ID, ymin = X0.025quant, ymax = X0.975quant, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) +
  geom_point(data = data.frame(res.inlabru.epsilon$summary.random$kappa), aes(x = ID, y = mean, fill = "Inlabru", color = "Inlabru")) + 
  theme_classic() + 
  scale_fill_manual(name = "", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title = "Kappa", x = "t", y = "")

p.random.eps <- (p.alpha.eps | p.beta.eps)/(p.kappa.eps) + plot_layout(guides = "collect")
ggsave("random_effects_eps.pdf", p.random.eps, path = output.path, dpi = "retina", width = 8, height = 5)

#   ----   Plot results - epsilon   ----

stan.epsilon.df.eps <- data.frame(list_of_draws_epsilon$epsilon)
plot.epsilon.inlabru.stan.compared(res.inlabru.epsilon, stan.epsilon.df.eps, path.to.storage = output.path, filename = "epsilon_eps", a45=T)

#   ----   Plot convergence diagnostics for stan   ----

p.neff.simplified <- ggplot(data = summary_stan) +
  geom_histogram(aes(x = n_eff), color = palette[1], fill = palette[1], alpha = 0.7, bins = 100) + 
  theme_classic() + 
  labs(title = "n_eff, simplified", x = "n_eff", y = "")

p.rhat.simplified <- ggplot(data = summary_stan) + 
  geom_histogram(aes(x = Rhat), color = palette[1], fill = palette[1], alpha = 0.5, bins = 100) + 
  theme_classic() + 
  labs(title = "Rhat, simplified", x = "Rhat", y = "")

p.convercence.simp <- (p.neff.simplified | p.rhat.simplified) + plot_layout(guides = "collect")

ggsave("convergence_simp.pdf", p.convercence.simp, path = output.path, dpi = "retina", height = 5, width = 8)

p.neff.epsilon <- ggplot(data = summary_stan_epsilon) +
  geom_histogram(aes(x = n_eff), color = palette[1], fill = palette[1], alpha = 0.7, bins = 100) + 
  theme_classic() + 
  labs(title = "n_eff, epsilon", x = "n_eff", y = "")

p.rhat.epsilon <- ggplot(data = summary_stan_epsilon) + 
  geom_histogram(aes(x = Rhat), color = palette[1], fill = palette[1], alpha = 0.5, bins = 100) + 
  theme_classic() + 
  labs(title = "Rhat, epsilon", x = "Rhat", y = "")

p.convercence.eps <- (p.neff.epsilon | p.rhat.epsilon) + plot_layout(guides = "collect")

ggsave("convergence_eps.pdf", p.convercence.eps, path = output.path, dpi = "retina", height = 5, width = 8)

