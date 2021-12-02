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

investigation.name <- "tllp_fh_rw1_extra_error_term_stan_7"

#   ----    Retrieve the data   ----

synthetic.male.lung.v7 <- function(){
  obs <- read.csv("Data/synthetic_male_lung_7.csv")
  obs <- obs %>% mutate(x.old = x, x = x - 9, x.c = x) %>%
    select(-X)
  
  obs.trad <- obs %>% 
    select(c(x, t, xt, age.int, year, x.c, alpha, beta, kappa, intercept, epsilon,
             eta, tau.alpha, tau.beta, tau.kappa, tau.epsilon, E)) %>%
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

#   ----   Run STAN analysis   ----
# Running traditional lc version of stan, 
# implemented with log-precisions. 

stan.output  <- file.path("Scripts/Synthetic data/Investigation", investigation.name)
source("Scripts/Synthetic\ data/run_stan_functions.R")

run_stan <- function(stan_program, obs, chains, warmup, iter, output.path, config.name, markov=TRUE){
  
  stan_fit <- run_stan_program_traditional_lc(
    list(obs = obs), chains=chains,warmup=warmup,
    iter=iter, stan_program=stan_program)
  
  store_stan_results_traditional.extra.error.term(
    fit=stan_fit, output.path=output.path, config=config.name,
    chains=chains, warmup=warmup, iter=iter, stan_program=stan_program,
    cohort=FALSE)
}

run_stan(
  stan_program="Scripts/Synthetic\ data/Stan\ analyses/stan_programs/stan_tllp_sc_fh_rw1_extra_error_term.stan",
  obs = obs, chains=4, warmup = 4000, iter = 40000, output.path = stan.output,
  config.name = investigation.name, markov=F)

inlabru.traditional.lc.fixed.hypers.rw1 <- function(obs, max_iter=30){
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
    alpha(x, model = "rw1", hyper = fixed.theta.alpha, constr = TRUE) +
    beta(x.c, model = "iid", extraconstr = list(A = A.beta, e = e.beta), hyper = fixed.theta.beta) +
    kappa(t, model = "rw1", constr = TRUE, hyper = fixed.theta.kappa)
  
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
                                   control.family  = c.family
                    ))
  return(res.inlabru)
}

res.inlabru <- inlabru.traditional.lc.fixed.hypers.rw1(obs, max_iter = 100)

source("Scripts/Functions/plotters.R")
source("Scripts/Synthetic data/plot_inlabru_vs_underlying.R")
source("Scripts/Synthetic data/plot_inlabru_stan_compared.R")
source("Scripts/Synthetic data/plot_stan_vs_underlying.R")

output.path <- stan.output

plots.summaries.inlabru <- plot.inlabru.vs.underlying.traditional.lc.fixed.effects(
  res.inlabru,
  underlying.effects,
  path.to.storage = output.path,
  save=T)

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
load(file.path(stan.output, "draws_error_term.RData"))

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
                               summaries.func=produce.summaries.stan.traditional)

plots_compared <- produce.compared.plots(
  stan.summaries = stan.res$summaries,
  stan.marginals = stan.marginals,
  inlabru.summaries = plots.summaries.inlabru$summaries,
  res.inlabru = res.inlabru,
  underlying.effects = underlying.effects,
  #plot.func = function(...) {plot.inlabru.stan.traditional.lc(..., cohort=FALSE, tau.beta.cutoff = 700, tau.kappa.cutoff = 500, tau.alpha.cutoff = 10, a45=F)},
  #plot.func = function(...) {plot.inlabru.stan.traditional.lc.no.beta(..., cohort=FALSE, tau.beta.cutoff = 5000, tau.kappa.cutoff = 5000, tau.alpha.cutoff = 100, a45=F)},
  plot.func = function(...) {plot.inlabru.stan.traditional.lc.fixed.hypers(..., cohort=FALSE, tau.beta.cutoff = 5000, tau.kappa.cutoff = 5000, tau.alpha.cutoff = 100, a45=F)},
  save.func = function(...) {save.compared.rw2(..., cohort=FALSE)},
  path.to.storage=output.path)

#   ----   Sample predictor   ----

stan.predictor.df <- data.frame(eta_draws)

plot.predictor.inlabru.stan.compared(res.inlabru, stan.predictor.df, path.to.storage = output.path, a45=T)

# ----   Plot offset.predictor   ----

offset.df <- data.frame(offset = res.inlabru$offset.linear.predictor[1:162], xt = obs$xt)

p.offset <- ggplot(offset.df) + geom_point(aes(x= xt, y= offset), color = palette[1]) + 
  labs("Offset of inlabru")
p.offset

save.figure(p.offset, name = "offset_inlabru", path = output.path, png = F)


#   ----   Plot difference between true eta and eta from stan   ----

p.eta.diff.stan <- ggplot(stan.res$summaries$summary_eta %>% mutate(diff = mean - true_eta)) +
  geom_point(aes(x = index, y = diff)) + 
  labs(title = "Difference between stan estimate of eta and true eta")
p.eta.diff.stan

save.figure(p.eta.diff.stan, name = "eta_diff_stan", path = output.path, png = F)

#   ----   Plot estimated error term   ----

error.term.stan.df <- stan_lc_df %>%
  rownames_to_column("parameter") %>%
  filter(grepl("error_term", parameter)) %>%
  filter(!grepl("tau_error_term", parameter)) %>%
  mutate(index = parse_number(parameter))

g.error.term <- ggplot(error.term.stan.df) + 
  geom_point(aes(x = index, y = mean), color = palette[1]) + 
  geom_errorbar(aes(x = index, ymin = `2.5%`, ymax = `97.5%`), color = palette[1]) + 
  theme_classic() + 
  labs(title = "Error terms, as estimated by stan")

g.error.term

save.figure(g.error.term, name = "error_term_stan", path = output.path, png = F)
