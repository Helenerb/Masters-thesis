# Running configuraton v10.3 with stan on Markov

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

run_stan_sml_rw1_7 <- function(stan_program, chains=4, warmup=1000, iter=10000, markov=TRUE){
  output.path <- set_workspace(config="synthetic_male_lung_7/rw1", markov)
  library("tidyverse")
  source("Scripts/Synthetic\ data/run_stan_functions.R")
  
  source("Scripts/Synthetic\ data/config_synthetic_male_lung_v7.R")
  
  data = synthetic.male.lung.a45.v7()
  
  stan_fit <- run_stan_program_lc(data$underlying.effects, chains=chains, warmup=warmup, iter=iter, stan_program=stan_program)
  
  store_stan_results(fit=stan_fit, output.path=output.path, config="synthetic_male_lung_7", chains=chains, warmup=warmup, iter=iter, stan_program=stan_program, cohort=FALSE)
}

run_stan_sml_rw1_7(stan_program="Scripts/Synthetic\ data/Stan\ analyses/stan_programs/stan_lc_rw1.stan", chains=1, warmup = 10000, iter = 100000, markov=T)
