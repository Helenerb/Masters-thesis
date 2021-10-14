# Running configuraton v17 with stan on Markov

# base script containing functions running STAN analysis

set_workspace <- function(config, markov=TRUE){
  if(markov){
    .libPaths("~/Documents/R_libraries")
    setwd("~/Documents/GitHub/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data")
    output.path <- file.path("~/Documents/GitHub/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses", config)
  } else {
    setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data")
    output.path <- file.path('~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses', config)
  }
  return(output.path)
}

run_stan_v17 <- function(stan_program, chains= 4, warmup = 1000, iter = 10000, markov = TRUE){
  output.path <- set_workspace("v17", markov)
  source("run_stan_functions.R")
  
  source("configurations_synthetic_data.R")
  
  data = configuration.v17()
  
  stan_fit <- run_stan_program_lcc(data, chains=chains, warmup=warmup, iter=iter, stan_program=stan_program)
  
  store_stan_results(fit=stan_fit, output.path=output.path, config="v17")
}

run_stan_v17("stan_analysis_lcc_phi_kappa.stan", chains=4, warmup = 1000, iter = 10000, markov=TRUE)
