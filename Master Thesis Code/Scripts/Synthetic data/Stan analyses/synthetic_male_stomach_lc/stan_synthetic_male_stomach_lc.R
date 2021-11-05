# Running configuraton v10.3 with stan on Markov

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

run_stan_male_stomach_lc <- function(stan_program, chains=4, warmup=1000, iter=10000, markov=TRUE){
  output.path <- set_workspace(config="synthetic_male_stomach_lc", markov)
  source("run_stan_functions.R")
  
  source("../Real\ data/synthetic_male_stomach_lc.R")
  
  data = synthetic.male.stomach.lc()
  
  stan_fit <- run_stan_program_lc(data, chains=chains, warmup=warmup, iter=iter, stan_program=stan_program)
  
  store_stan_results(fit=stan_fit, output.path=output.path, config="synthetic_male_stomach_lc", chains=chains, warmup=warmup, iter=iter, stan_program=stan_program, cohort=FALSE)
}

run_stan_male_stomach_lc(stan_program="Stan\ analyses/stan_programs/stan_analysis_lc_rw2.stan", chains=4, warmup = 1000, iter = 10000, markov=TRUE)
