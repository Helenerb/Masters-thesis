#' This program investigates inference with:
#' - The traditional_lc_log_prec version - Gaussian lc-model. 
#' All random effects are modelled as iid, in inlabru and in stan
#' Constraints are imposed as usual in inlabru, and with soft constraints in stan
#' The hyperparameters are fixed. 
#' 

#   ----   Load libraries and set workspace   ----
library("tidyverse")
library("inlabru")
library("ggplot2")
library("INLA")
library("patchwork")
library("rstan")

setwd("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master\ Thesis\ Code")

investigation.name <- "sml_fh_all_iid_no_constr_7"

#   ----    Retrieve the data   ----

synthetic.male.lung.v7 <- function(){
  obs <- read.csv("Data/synthetic_male_lung_7.csv")
  obs <- obs %>% mutate(x.old = x, x = x - 9, x.c = x)
  
  underlying.effects <- list(obs = obs, nx = 9, nt = 18,
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

stan.output  <- file.path("Scripts/Synthetic data/Investigation", investigation.name)
source("Scripts/Synthetic\ data/run_stan_functions.R")

run_stan <- function(stan_program, obs, chains, warmup, iter, output.path, config.name, markov=TRUE){
  
  stan_fit <- run_stan_program_lc(
    list(obs = obs), chains=chains,warmup=warmup,
    iter=iter, stan_program=stan_program)
  
  store_stan_results_traditional(
    fit=stan_fit, output.path=output.path, config=config.name,
    chains=chains, warmup=warmup, iter=iter, stan_program=stan_program,
    cohort=FALSE)
}

run_stan(
  stan_program="Scripts/Synthetic\ data/Stan\ analyses/stan_programs/stan_lc_fh_iid_no_constr.stan",
  obs = obs, chains=4, warmup = 4000, iter = 40000, output.path = stan.output,
  config.name = investigation.name, markov=F)

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
  
  c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE, return.marginals.predictor = TRUE)  # control.compute
  
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

res.inlabru <- inlabru.lc.fh.iid.no.constr(obs, max_iter = 100)

source("Scripts/Functions/plotters.R")
source("Scripts/Synthetic data/plot_inlabru_vs_underlying.R")
source("Scripts/Synthetic data/plot_inlabru_stan_compared.R")
source("Scripts/Synthetic data/plot_stan_vs_underlying.R")

output.path <- stan.output

plots.summaries.inlabru <- plot.inlabru.vs.underlying.synthetic.cancer.fixed.effects(
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
                               save.func=save.stan.plots.lc.rw2,
                               path.to.storage=output.path,
                               summaries.func=produce.summaries.stan.lc.rw2)


plots_compared <- produce.compared.plots(
  stan.summaries = stan.res$summaries,
  stan.marginals = stan.marginals,
  inlabru.summaries = plots.summaries.inlabru$summaries,
  res.inlabru = res.inlabru,
  underlying.effects = underlying.effects,
  #plot.func = function(...) {plot.inlabru.stan.traditional.lc(..., cohort=FALSE, tau.beta.cutoff = 700, tau.kappa.cutoff = 500, tau.alpha.cutoff = 10, a45=F)},
  #plot.func = function(...) {plot.inlabru.stan.traditional.lc.no.beta(..., cohort=FALSE, tau.beta.cutoff = 5000, tau.kappa.cutoff = 5000, tau.alpha.cutoff = 100, a45=F)},
  plot.func = function(...) {plot.inlabru.stan.compared.rw2(..., cohort=FALSE, tau.beta.cutoff = 5000, tau.kappa.cutoff = 5000, tau.alpha.cutoff = 100, a45=F)},
  save.func = function(...) {save.compared.rw2(..., cohort=FALSE)},
  path.to.storage=output.path)


#   ----   Sample predictor   ----

inlabru.samps.predictor <- generate(
  res.inlabru,
  data = data.frame(x = obs$x, t = obs$t, x.c = obs$x.c, xt = obs$xt),
  formula = ~ Int + alpha + beta*kappa + epsilon,
  n.sample = 10000)

inlabru.predictor.df <- data.frame(t(inlabru.samps.predictor))

stan.samps.predictor <- eta_draws[sample(nrow(eta_draws), size = 10000, replace = F),]

stan.predictor.df <- data.frame(eta_draws)

plot.predictor.inlabru.stan.compared(inlabru.predictor.df, stan.predictor.df, path.to.storage = output.path, a45=T)

