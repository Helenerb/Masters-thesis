# Running traditional lc version of stan, 
# implemented with log-precisions. 

set_workspace <- function(config, markov=TRUE){
  if(markov){
    .libPaths("~/Documents/R_libraries")
    setwd("~/Documents/GitHub/Masteroppgave/Masters-thesis/Master Thesis Code")
    output.path <- file.path("~/Documents/GitHub/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses", config)
  } else {
    setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code")
    output.path <- file.path('~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses', config)
  }
  return(output.path)
}

run_stan_tllp_v4_unconstrained <- function(stan_program, chains=4, warmup=1000, iter=10000, markov=TRUE){
  output.path <- set_workspace(config="traditional_lc_log_prec/v4_unconstrained", markov)
  library("tidyverse")
  source("Scripts/Synthetic\ data/run_stan_functions.R")
  
  source("Scripts/Synthetic\ data/config_synthetic_male_lung_v4.R")
  
  data = synthetic.male.lung.v4()
  
  # reformat
  obs.lc.poiss <- data$underlying.effects$obs
  
  print(head(obs.lc.poiss))
  
  obs.trad <- obs.lc.poiss %>% 
    select(c(x, t, xt, age.int, year, x.c, alpha, beta, kappa, intercept, epsilon,
             eta, tau.alpha, tau.beta, tau.kappa, tau.epsilon, E)) %>%
    mutate(mr_gaussian = exp(eta)) %>%
    mutate(eta.no.error = intercept + alpha + beta*kappa) %>%
    mutate(Y_gaussian  = mr_gaussian * E)
  
  print(head(obs.trad))
  
  stan_fit <- run_stan_program_traditional_lc(list(obs = obs.trad), chains=chains, warmup=warmup, iter=iter, stan_program=stan_program)
  
  store_stan_results_traditional(fit=stan_fit, output.path=output.path, config="traditional_lc_log_prec", chains=chains, warmup=warmup, iter=iter, stan_program=stan_program, cohort=FALSE)
}

run_stan_tllp_v4_unconstrained(stan_program="Scripts/Synthetic\ data/Stan\ analyses/stan_programs/stan_traditional_lc_no_constraints.stan", chains=2, warmup = 20000, iter = 200000, markov=T)
