# running the v10 configuration with the stan program where
# the period effect is modeled as a sum of a linear term and a driftless random
# walk

# Running configuraton v10 with stan on Markov

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

run_stan_v10 <- function(stan_program, chains=4, warmup=1000, iter=10000, markov=TRUE){
  output.path <- set_workspace(config="v10", markov)
  source("run_stan_functions.R")
  
  source("configurations_synthetic_data.R")
  
  data = configuration.v10.2()
  
  stan_fit <- run_stan_program_lc(data, chains=chains, warmup=warmup, iter=iter, stan_program=stan_program)
  
  store_stan_results(fit=stan_fit, output.path=output.path, config="v10")
}

run_stan_v10(stan_program="Stan\ analyses/stan_programs/stan_analysis_lc_v4.stan", chains=4, warmup = 4000, iter = 40000, markov=TRUE)
