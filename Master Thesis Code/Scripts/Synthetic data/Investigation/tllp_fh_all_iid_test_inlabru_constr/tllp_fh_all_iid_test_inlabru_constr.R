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

# TODO: Change to where you want to run from
setwd("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master\ Thesis\ Code")

investigation.name <- "tllp_fh_all_iid_test_inlabru_constr"

#   ----    Retrieve the data   ----

synthetic.male.lung.v7 <- function(){
  # TODO: Change to where you have stored the csv file
  obs <- read.csv("Data/synthetic_male_lung_7.csv")
  obs <- obs %>% mutate(x.old = x, x = x - 9, x.c = x) %>%
    select(-X)
  
  obs.trad <- obs %>% 
    select(c(x, t, xt, age.int, year, x.c, alpha, beta, kappa, intercept, epsilon,
             eta, tau.alpha, tau.beta, tau.kappa, tau.epsilon, E, Y)) %>%
    #mutate(eta = eta) %>%
    mutate(eta.no.error = intercept + alpha + beta*kappa) %>%
    mutate(mr_gaussian = exp(eta)) %>%
    mutate(Y_gaussian  = mr_gaussian * E)
  
  underlying.effects <- list(obs = obs.trad, nx = 9, nt = 18,
                             alpha.true = {obs %>% filter(t == 0)}$alpha,
                             beta.true = {obs %>% filter(t == 0)}$beta,
                             kappa.true = {obs %>% filter(x == 0)}$kappa,
                             intercept = unique(obs$intercept),
                             age.intercept.true = unique(obs$intercept),
                             tau.alpha.true = unique(obs$tau.alpha),
                             tau.beta.true = unique(obs$tau.beta),
                             tau.kappa.true = unique(obs$tau.kappa),
                             tau.epsilon.true = unique(obs$tau.epsilon))
  
  return(list(obs = obs.trad, underlying.effects = underlying.effects))
}

# We use this data for both inlabru and stan
config.data <- synthetic.male.lung.v7()
obs <- config.data$obs
underlying.effects <- config.data$underlying.effects

# Gaussian models

inlabru.traditional.lc.fixed.hypers.all.iid.no.constr <- function(obs, max_iter=30){
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
  fixed.theta.beta <- list(prec = list(initial = log(100), fixed = T))
  fixed.theta.kappa <- list(prec = list(initial = log(70), fixed = T))
  fixed.theta.epsilon <- list(prec = list(initial = log(400), fixed = T))
  
  comp = ~ -1 +
    Int(1, prec.linear = 0.001, mean.linear = 0) +
    alpha(x, model = "iid", hyper = fixed.theta.alpha, constr = FALSE) +
    beta(x.c, model = "iid", hyper = fixed.theta.beta, constr = FALSE) +
    kappa(t, model = "iid", hyper = fixed.theta.kappa, constr = FALSE)
  
  formula = eta ~ Int + alpha + beta*kappa
  
  likelihood = like(formula = formula, family = "gaussian", data = obs)
  
  c.compute <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)  # control.compute
  c.family <- list(hyper = fixed.theta.epsilon)
  
  res.inlabru = bru(components = comp,
                    likelihood, 
                    options = list(verbose = F,
                                   bru_verbose = 1, 
                                   num.threads = "1:1",
                                   control.compute = c.compute,
                                   bru_max_iter=max_iter,
                                   #control.predictor = list(compute = TRUE),
                                   control.family  = c.family
                    ))
  return(res.inlabru)
}

inlabru.traditional.lc.fixed.hypers.all.iid <- function(obs, max_iter=30){
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
  fixed.theta.beta <- list(prec = list(initial = log(100), fixed = T))
  fixed.theta.kappa <- list(prec = list(initial = log(70), fixed = T))
  fixed.theta.epsilon <- list(prec = list(initial = log(400), fixed = T))
  
  comp = ~ -1 +
    Int(1, prec.linear = 0.001, mean.linear = 0) +
    alpha(x, model = "iid", hyper = fixed.theta.alpha, constr = TRUE) +
    beta(x.c, model = "iid", hyper = fixed.theta.beta, extraconstr = list(A = A.beta, e = e.beta)) +
    kappa(t, model = "iid", hyper = fixed.theta.kappa, constr = TRUE)
  
  formula = eta ~ Int + alpha + beta*kappa
  
  likelihood = like(formula = formula, family = "gaussian", data = obs)
  
  c.compute <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)  # control.compute
  c.family <- list(hyper = fixed.theta.epsilon)
  
  res.inlabru = bru(components = comp,
                    likelihood, 
                    options = list(verbose = F,
                                   bru_verbose = 1, 
                                   num.threads = "1:1",
                                   control.compute = c.compute,
                                   bru_max_iter=max_iter,
                                   control.predictor = list(compute = TRUE),
                                   control.family  = c.family
                    ))
  return(res.inlabru)
}

inlabru.traditional.lc.fixed.hypers.all.iid.only.beta.constr <- function(obs, max_iter=30){
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
  fixed.theta.beta <- list(prec = list(initial = log(100), fixed = T))
  fixed.theta.kappa <- list(prec = list(initial = log(70), fixed = T))
  fixed.theta.epsilon <- list(prec = list(initial = log(400), fixed = T))
  
  comp = ~ -1 +
    Int(1, prec.linear = 0.001, mean.linear = 0) +
    alpha(x, model = "iid", hyper = fixed.theta.alpha, constr = FALSE) +
    beta(x.c, model = "iid", hyper = fixed.theta.beta, extraconstr = list(A = A.beta, e = e.beta)) +
    kappa(t, model = "iid", hyper = fixed.theta.kappa, constr = FALSE)
  
  formula = eta ~ Int + alpha + beta*kappa
  
  likelihood = like(formula = formula, family = "gaussian", data = obs)
  
  c.compute <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)  # control.compute
  c.family <- list(hyper = fixed.theta.epsilon)
  
  res.inlabru = bru(components = comp,
                    likelihood, 
                    options = list(verbose = F,
                                   bru_verbose = 1, 
                                   num.threads = "1:1",
                                   control.compute = c.compute,
                                   bru_max_iter=max_iter,
                                   control.predictor = list(compute = TRUE),
                                   control.family  = c.family
                    ))
  return(res.inlabru)
}

inlabru.traditional.lc.fixed.hypers.all.iid.alpha.kappa.constr <- function(obs, max_iter=30){
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
  fixed.theta.beta <- list(prec = list(initial = log(100), fixed = T))
  fixed.theta.kappa <- list(prec = list(initial = log(70), fixed = T))
  fixed.theta.epsilon <- list(prec = list(initial = log(400), fixed = T))
  
  comp = ~ -1 +
    Int(1, prec.linear = 0.001, mean.linear = 0) +
    alpha(x, model = "iid", hyper = fixed.theta.alpha, constr = TRUE) +
    beta(x.c, model = "iid", hyper = fixed.theta.beta, constr = FALSE) +
    kappa(t, model = "iid", hyper = fixed.theta.kappa, constr = TRUE)
  
  formula = eta ~ Int + alpha + beta*kappa
  
  likelihood = like(formula = formula, family = "gaussian", data = obs)
  
  c.compute <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE, return.marginals.predictor = TRUE)  # control.compute
  c.family <- list(hyper = fixed.theta.epsilon)
  
  res.inlabru = bru(components = comp,
                    likelihood, 
                    options = list(verbose = F,
                                   bru_verbose = 1, 
                                   num.threads = "1:1",
                                   control.compute = c.compute,
                                   bru_max_iter=max_iter,
                                   control.predictor = list(compute = TRUE),
                                   control.family  = c.family
                    ))
  return(res.inlabru)
}

inlabru.traditional.lc.fixed.hypers.some.rw1.no.constr <- function(obs, max_iter=30){
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
  fixed.theta.beta <- list(prec = list(initial = log(100), fixed = T))
  fixed.theta.kappa <- list(prec = list(initial = log(70), fixed = T))
  fixed.theta.epsilon <- list(prec = list(initial = log(400), fixed = T))
  
  comp = ~ -1 +
    Int(1, prec.linear = 0.001, mean.linear = 0) +
    alpha(x, model = "rw1", hyper = fixed.theta.alpha, constr = FALSE) +
    beta(x.c, model = "iid", hyper = fixed.theta.beta, constr = FALSE) +
    kappa(t, model = "rw1", hyper = fixed.theta.kappa, constr = FALSE)
  
  formula = eta ~ Int + alpha + beta*kappa
  
  likelihood = like(formula = formula, family = "gaussian", data = obs)
  
  c.compute <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE, return.marginals.predictor = TRUE)  # control.compute
  c.family <- list(hyper = fixed.theta.epsilon)
  
  res.inlabru = bru(components = comp,
                    likelihood, 
                    options = list(verbose = F,
                                   bru_verbose = 1, 
                                   num.threads = "1:1",
                                   control.compute = c.compute,
                                   bru_max_iter=max_iter,
                                   control.predictor = list(compute = TRUE),
                                   control.family  = c.family
                    ))
  return(res.inlabru)
}

inlabru.traditional.lc.fixed.hypers.all.iid.no.constr.extra.epsilon <- function(obs, max_iter=30){
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
  fixed.theta.beta <- list(prec = list(initial = log(100), fixed = T))
  fixed.theta.kappa <- list(prec = list(initial = log(70), fixed = T))
  fixed.theta.epsilon <- list(prec = list(initial = log(400), fixed = T))
  
  comp = ~ -1 +
    Int(1, prec.linear = 0.001, mean.linear = 0) +
    alpha(x, model = "iid", hyper = fixed.theta.alpha, constr = FALSE) +
    beta(x.c, model = "iid", hyper = fixed.theta.beta, constr = FALSE) +
    kappa(t, model = "iid", hyper = fixed.theta.kappa, constr = FALSE) + 
    epsilon(xt, model="iid", hyper = fixed.theta.epsilon, constr = FALSE)
  
  formula = eta ~ Int + alpha + beta*kappa + epsilon
  
  likelihood = like(formula = formula, family = "gaussian", data = obs)
  
  c.compute <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)  # control.compute
  c.family <- list(hyper = fixed.theta.epsilon)
  
  res.inlabru = bru(components = comp,
                    likelihood, 
                    options = list(verbose = F,
                                   bru_verbose = 1, 
                                   num.threads = "1:1",
                                   control.compute = c.compute,
                                   bru_max_iter=max_iter,
                                   #control.predictor = list(compute = TRUE),
                                   control.family  = c.family
                    ))
  return(res.inlabru)
}

# Poisson models

inlabru.lc.fh.iid.no.constr <- function(obs, max_iter=30){
  #'Implements inlabru analysis for lc model using an ar1c to model the period effect
  #'
  #'@param obs: Contains the observed data and the real underlying random effects
  #'@param max_iter (int): maximum number of iterations in inlabru
  
  #nx = length(unique(obs$x))
  #nt = length(unique(obs$t))
  
  # constraints for the age effect beta
  #A.beta = matrix(1, nrow = 1, ncol = nx)  
  #e.beta = 1  
  
  fixed.theta.alpha <- list(prec = list(initial = log(1.96), fixed = T))
  fixed.theta.beta <- list(prec = list(initial = log(100), fixed = T))
  fixed.theta.kappa <- list(prec = list(initial = log(70), fixed = T))
  fixed.theta.epsilon <- list(prec = list(initial = log(400), fixed = T))
  
  comp = ~ -1 +
    Int(1, prec.linear = 0.001, mean.linear = 0) +
    alpha(x, model = "iid", hyper = fixed.theta.alpha, constr = FALSE) +
    beta(x.c, model = "iid", hyper = fixed.theta.beta, constr = FALSE) +
    kappa(t, model = "iid", hyper = fixed.theta.kappa, constr = FALSE) +
    epsilon(xt, model = "iid", hyper = fixed.theta.epsilon, constr = FALSE)
  
  formula = Y ~ Int + alpha + beta*kappa + epsilon
  
  likelihood = like(formula = formula, family = "poisson", data = obs, E = obs$E)
  
  c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)  # control.compute
  
  res.inlabru = bru(components = comp,
                    likelihood, 
                    options = list(verbose = F,
                                   bru_verbose = 1, 
                                   num.threads = "1:1",
                                   control.compute = c.c,
                                   bru_max_iter=max_iter,
                                   control.predictor = list(link = 1)
                    ))
  return(res.inlabru)
}

inlabru.lc.fh.iid.beta.constr <- function(obs, max_iter=30){
  #'Implements inlabru analysis for lc model using an ar1c to model the period effect
  #'
  #'@param obs: Contains the observed data and the real underlying random effects
  #'@param max_iter (int): maximum number of iterations in inlabru
  
  nx = length(unique(obs$x))
  nt = length(unique(obs$t))
  
  # constraints for the age effect beta
  A.beta = matrix(1, nrow = 1, ncol = nx)
  e.beta = 1
  
  fixed.theta.alpha <- list(prec = list(initial = log(1.96), fixed = T))
  fixed.theta.beta <- list(prec = list(initial = log(100), fixed = T))
  fixed.theta.kappa <- list(prec = list(initial = log(70), fixed = T))
  fixed.theta.epsilon <- list(prec = list(initial = log(400), fixed = T))
  
  comp = ~ -1 +
    Int(1, prec.linear = 0.001, mean.linear = 0) +
    alpha(x, model = "iid", hyper = fixed.theta.alpha, constr = FALSE) +
    beta(x.c, model = "iid", hyper = fixed.theta.beta, extraconstr = list(A = A.beta, e = e.beta)) +
    kappa(t, model = "iid", hyper = fixed.theta.kappa, constr = FALSE) +
    epsilon(xt, model = "iid", hyper = fixed.theta.epsilon, constr = FALSE)
  
  formula = Y ~ Int + alpha + beta*kappa + epsilon
  
  likelihood = like(formula = formula, family = "poisson", data = obs, E = obs$E)
  
  c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)  # control.compute
  
  res.inlabru = bru(components = comp,
                    likelihood, 
                    options = list(verbose = F,
                                   bru_verbose = 1, 
                                   num.threads = "1:1",
                                   control.compute = c.c,
                                   bru_max_iter=max_iter,
                                   control.predictor = list(link = 1)
                    ))
  return(res.inlabru)
}

inlabru.lc.fh.iid.all.constr <- function(obs, max_iter=30){
  #'Implements inlabru analysis for lc model using an ar1c to model the period effect
  #'
  #'@param obs: Contains the observed data and the real underlying random effects
  #'@param max_iter (int): maximum number of iterations in inlabru
  
  nx = length(unique(obs$x))
  nt = length(unique(obs$t))
  
  # constraints for the age effect beta
  A.beta = matrix(1, nrow = 1, ncol = nx)
  e.beta = 1
  
  fixed.theta.alpha <- list(prec = list(initial = log(1.96), fixed = T))
  fixed.theta.beta <- list(prec = list(initial = log(100), fixed = T))
  fixed.theta.kappa <- list(prec = list(initial = log(70), fixed = T))
  fixed.theta.epsilon <- list(prec = list(initial = log(400), fixed = T))
  
  comp = ~ -1 +
    Int(1, prec.linear = 0.001, mean.linear = 0) +
    alpha(x, model = "iid", hyper = fixed.theta.alpha, constr = TRUE) +
    beta(x.c, model = "iid", hyper = fixed.theta.beta, extraconstr = list(A = A.beta, e = e.beta)) +
    kappa(t, model = "iid", hyper = fixed.theta.kappa, constr = TRUE) +
    epsilon(xt, model = "iid", hyper = fixed.theta.epsilon, constr = FALSE)
  
  formula = Y ~ Int + alpha + beta*kappa + epsilon
  
  likelihood = like(formula = formula, family = "poisson", data = obs, E = obs$E)
  
  c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)  # control.compute
  
  res.inlabru = bru(components = comp,
                    likelihood, 
                    options = list(verbose = F,
                                   bru_verbose = 1, 
                                   num.threads = "1:1",
                                   control.compute = c.c,
                                   bru_max_iter=max_iter,
                                   control.predictor = list(link = 1)
                    ))
  return(res.inlabru)
}

inlabru.lc.fh.iid.no.constr.no.epsilon <- function(obs, max_iter=30){
  #'Implements inlabru analysis for lc model using an ar1c to model the period effect
  #'
  #'@param obs: Contains the observed data and the real underlying random effects
  #'@param max_iter (int): maximum number of iterations in inlabru
  
  #nx = length(unique(obs$x))
  #nt = length(unique(obs$t))
  
  # constraints for the age effect beta
  #A.beta = matrix(1, nrow = 1, ncol = nx)  
  #e.beta = 1  
  
  fixed.theta.alpha <- list(prec = list(initial = log(1.96), fixed = T))
  fixed.theta.beta <- list(prec = list(initial = log(100), fixed = T))
  fixed.theta.kappa <- list(prec = list(initial = log(70), fixed = T))
  #fixed.theta.epsilon <- list(prec = list(initial = log(400), fixed = T))
  
  comp = ~ -1 +
    Int(1, prec.linear = 0.001, mean.linear = 0) +
    alpha(x, model = "iid", hyper = fixed.theta.alpha, constr = FALSE) +
    beta(x.c, model = "iid", hyper = fixed.theta.beta, constr = FALSE) +
    kappa(t, model = "iid", hyper = fixed.theta.kappa, constr = FALSE)
    #epsilon(xt, model = "iid", hyper = fixed.theta.epsilon, constr = FALSE)
  
  formula = Y ~ Int + alpha + beta*kappa
  
  likelihood = like(formula = formula, family = "poisson", data = obs, E = obs$E)
  
  c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)  # control.compute
  
  res.inlabru = bru(components = comp,
                    likelihood, 
                    options = list(verbose = F,
                                   bru_verbose = 1, 
                                   num.threads = "1:1",
                                   control.compute = c.c,
                                   bru_max_iter=max_iter,
                                   control.predictor = list(link = 1)
                    ))
  return(res.inlabru)
}

# Inlabru results - gaussian

res.inlabru.no.constr <- inlabru.traditional.lc.fixed.hypers.all.iid.no.constr(obs, max_iter = 100)

res.inlabru.constrained <- inlabru.traditional.lc.fixed.hypers.all.iid(obs, max_iter = 100)

res.inlabru.some.rw1.no.constr <- inlabru.traditional.lc.fixed.hypers.some.rw1.no.constr(obs, max_iter = 100)

res.inlabru.all.iid.beta.constr <- inlabru.traditional.lc.fixed.hypers.all.iid.only.beta.constr(obs, max_iter = 100)

res.inlabru.all.iid.alpha.kappa.constr <- inlabru.traditional.lc.fixed.hypers.all.iid.alpha.kappa.constr(obs, max_iter = 100)

res.inlabru.all.iid.no.constr.extra.epsilon <- inlabru.traditional.lc.fixed.hypers.all.iid.no.constr.extra.epsilon(obs, max_iter = 100)

# Inlabru results - poisson

res.inlabru.poisson.no.constr <- inlabru.lc.fh.iid.no.constr(obs, max_iter = 30)

res.inlabru.poisson.beta.constr <- inlabru.lc.fh.iid.beta.constr(obs, max_iter = 30)

res.inlabru.poisson.all.constr <- inlabru.lc.fh.iid.all.constr(obs, max_iter = 30)

res.inlabru.poisson.no.constr.no.epsilon <- inlabru.lc.fh.iid.no.constr.no.epsilon(obs, max_iter = 30)


#   -----    For plotting  - requires cloned GitHub repo and correct output paths -----


source("Scripts/Functions/plotters.R")
source("Scripts/Synthetic data/plot_inlabru_vs_underlying.R")
source("Scripts/Synthetic data/plot_inlabru_stan_compared.R")
source("Scripts/Synthetic data/plot_stan_vs_underlying.R")

output.path <- file.path("Scripts/Synthetic data/Investigation", investigation.name)

plots.summaries.inlabru.no.constr <- plot.inlabru.vs.underlying.traditional.lc.fixed.effects(
  res.inlabru.no.constr,
  underlying.effects,
  path.to.storage = file.path(output.path, "no\ constr"),
  save=T, png = F)

plots.summaries.inlabru.constrained <- plot.inlabru.vs.underlying.traditional.lc.fixed.effects(
  res.inlabru.constrained,
  underlying.effects,
  path.to.storage = file.path(output.path, "constr"),
  save=T, png= F)

plots.summaries.inlabru.some.rw1.no.constr <- plot.inlabru.vs.underlying.traditional.lc.fixed.effects(
  res.inlabru.some.rw1.no.constr,
  underlying.effects,
  path.to.storage = file.path(output.path, "no\ constr\ rw1"),
  save=T, png = F)

plots.summaries.inlabru.all.iid.beta.constr <- plot.inlabru.vs.underlying.traditional.lc.fixed.effects(
  res.inlabru.all.iid.beta.constr,
  underlying.effects,
  path.to.storage = file.path(output.path, "only\ beta\ constr"),
  save=T, pdf = T, png = F)

plots.summaries.inlabru.all.iid.alpha.kappa.constr <- plot.inlabru.vs.underlying.traditional.lc.fixed.effects(
  res.inlabru.all.iid.alpha.kappa.constr,
  underlying.effects,
  path.to.storage = file.path(output.path, "alpha\ kappa\ constr"),
  save=T, pdf = T, png = F)

plots.summaries.inlabru.poisson.no.constraints <- plot.inlabru.vs.underlying.traditional.lc.fixed.effects(
  res.inlabru.poisson.no.constr,
  underlying.effects,
  path.to.storage = file.path(output.path, "poisson\ no\ constr"),
  save=T, pdf = T, png = F)

plots.summaries.inlabru.poisson.beta.constraints <- plot.inlabru.vs.underlying.traditional.lc.fixed.effects(
  res.inlabru.poisson.beta.constr,
  underlying.effects,
  path.to.storage = file.path(output.path, "poisson\ beta\ constr"),
  save=T, pdf = T, png = F)

plots.summaries.inlabru.poisson.no.constraints.no.epsilon <- plot.inlabru.vs.underlying.traditional.lc.fixed.effects(
  res.inlabru.poisson.no.constr.no.epsilon,
  underlying.effects,
  path.to.storage = file.path(output.path, "poisson\ no\ constr\ no\ epsilon"),
  save=T, pdf = T, png = F)

plots.summaries.inlabru.all.iid.no.constr.extra.epsilon <- plot.inlabru.vs.underlying.traditional.lc.fixed.effects(
  res.inlabru.all.iid.no.constr.extra.epsilon,
  underlying.effects,
  path.to.storage = file.path(output.path, "no\ constr\ extra\ epsilon"),
  save=T, pdf = T, png = F)

plots.summaries.inlabru.all.iid.no.constr.extra.epsilon <- plot.inlabru.vs.underlying.traditional.lc.fixed.effects(
  res.inlabru.poisson.all.constr,
  underlying.effects,
  path.to.storage = file.path(output.path, "poisson\ all\ constr"),
  save=T, pdf = T, png = F)


#   ----    Sample eta.no.error from Poisson distributions  ----

# Poisson, no constraints
samps.Poisson.eta.no.error <- generate(
  res.inlabru.poisson.no.constr,
  data = data.frame(x = obs$x, t = obs$t, x.c = obs$x.c, xt = obs$xt),
  formula = ~ Int + alpha + beta*kappa,
  n.sample = 10000)

samps.Poisson.eta.no.error.df <- data.frame(
  xt = obs$xt, 
  mean = apply(samps.Poisson.eta.no.error, 1, mean),
  X0.025quant = apply(samps.Poisson.eta.no.error, 1, quantile, 0.025),
  X0.975quant = apply(samps.Poisson.eta.no.error, 1, quantile, 0.975)) %>%
  left_join(obs, by = c("xt" = "xt"), suffix = c("", ".obs"))

p.eta.no.error <- ggplot(samps.Poisson.eta.no.error.df) + 
  geom_point(aes(x = xt, y = eta, color = "Eta with error", fill = "Eta with error"), shape = 4) + 
  geom_point(aes(x = xt, y = mean, color = "Inlabru", fill = "Inlabru")) + 
  geom_ribbon(aes(x = xt, ymin = X0.025quant, ymax = X0.975quant, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_point(aes(x = xt, y = eta.no.error, color = "True", fill = "True")) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Estimated and true Int + alpha + beta*kappa (eta no error)")

p.eta.no.error

save.figure(p.eta.no.error, name = "eta_no_error", path = file.path(output.path, "poisson\ no\ constr"), png=F)

# Poisson, with constraints on beta:

samps.Poisson.beta.constr.eta.no.error <- generate(
  res.inlabru.poisson.beta.constr,
  data = data.frame(x = obs$x, t = obs$t, x.c = obs$x.c, xt = obs$xt),
  formula = ~ Int + alpha + beta*kappa,
  n.sample = 10000)

samps.Poisson.beta.constr.eta.no.error.df <- data.frame(
  xt = obs$xt, 
  mean = apply(samps.Poisson.beta.constr.eta.no.error, 1, mean),
  X0.025quant = apply(samps.Poisson.beta.constr.eta.no.error, 1, quantile, 0.025),
  X0.975quant = apply(samps.Poisson.beta.constr.eta.no.error, 1, quantile, 0.975)) %>%
  left_join(obs, by = c("xt" = "xt"), suffix = c("", ".obs"))

p.eta.no.error.beta.constr <- ggplot(samps.Poisson.beta.constr.eta.no.error.df) + 
  geom_point(aes(x = xt, y = mean, color = "Inlabru", fill = "Inlabru")) + 
  geom_ribbon(aes(x = xt, ymin = X0.025quant, ymax = X0.975quant, fill = "Inlabru"), alpha = 0.4) + 
  geom_point(aes(x = xt, y = eta.no.error, color = "True", fill = "True")) + 
  geom_point(aes(x = xt, y = eta, color = "True, with error", fill = "True, with error"), shape = 4) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Estimated and true Int + alpha + beta*kappa (eta no error)")

p.eta.no.error.beta.constr

save.figure(p.eta.no.error.beta.constr, name = "eta_no_error", path = file.path(output.path, "poisson\ beta\ constr"), png=F)

# Poisson, with all constraints:

samps.Poisson.all.constr.eta.no.error <- generate(
  res.inlabru.poisson.all.constr,
  data = data.frame(x = obs$x, t = obs$t, x.c = obs$x.c, xt = obs$xt),
  formula = ~ Int + alpha + beta*kappa,
  n.sample = 10000)

samps.Poisson.all.constr.eta.no.error.df <- data.frame(
  xt = obs$xt, 
  mean = apply(samps.Poisson.all.constr.eta.no.error, 1, mean),
  X0.025quant = apply(samps.Poisson.all.constr.eta.no.error, 1, quantile, 0.025),
  X0.975quant = apply(samps.Poisson.all.constr.eta.no.error, 1, quantile, 0.975)) %>%
  left_join(obs, by = c("xt" = "xt"), suffix = c("", ".obs"))

p.eta.no.error.all.constr <- ggplot(samps.Poisson.all.constr.eta.no.error.df) + 
  geom_point(aes(x = xt, y = mean, color = "Inlabru", fill = "Inlabru")) + 
  geom_ribbon(aes(x = xt, ymin = X0.025quant, ymax = X0.975quant, fill = "Inlabru"), alpha = 0.4) + 
  geom_point(aes(x = xt, y = eta.no.error, color = "True", fill = "True")) + 
  geom_point(aes(x = xt, y = eta, color = "True, with error", fill = "True, with error"), shape = 4) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Estimated and true Int + alpha + beta*kappa (eta no error)")

p.eta.no.error.all.constr

save.figure(p.eta.no.error.all.constr, name = "eta_no_error", path = file.path(output.path, "poisson\ all\ constr"), png=F)




