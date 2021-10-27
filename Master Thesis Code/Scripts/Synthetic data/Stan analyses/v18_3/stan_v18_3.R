# running the v18_3 configuration with the stan program where
# the period effect is modelled as a random walk of order two, summed to zero

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

run_stan_v18_3 <- function(stan_program, chains=4, warmup=1000, iter=10000, markov=TRUE){
  output.path <- set_workspace(config="v18_3", markov)
  source("run_stan_functions.R")
  
  source("configurations_synthetic_data.R")
  
  data = configuration.v18.3()
  
  stan_fit <- run_stan_program_lcc(data, chains=chains, warmup=warmup, iter=iter, stan_program=stan_program)
  
  store_stan_results(fit=stan_fit, output.path=output.path, config="v18_3")
}

run_stan_v18_3(stan_program="Stan\ analyses/stan_programs/stan_analysis_cohort_rw2.stan", chains=4, warmup = 50, iter = 50, markov=TRUE)
