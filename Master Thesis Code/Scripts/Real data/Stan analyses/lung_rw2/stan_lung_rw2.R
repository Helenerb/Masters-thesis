# runs stan analysis for lung cancer data, adapted to be runable
# on Markov. 

set_workspace <- function(config, markov=TRUE){
  if(markov){
    .libPaths("~/Documents/R_libraries")
    setwd("~/Documents/GitHub/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Real data")
    output.path <- file.path("~/Documents/GitHub/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Real data/Stan analyses", config)
  } else {
    setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Real data")
    output.path <- file.path('~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Real data/Stan analyses', config)
  }
  return(output.path)
}

run_stan_v18_3 <- function(stan_program, chains=4, warmup=1000, iter=10000, markov=TRUE){
  config = "lung_rw2"
  output.path <- set_workspace(config=config, markov)
  source("../Synthetic\ data/run_stan_functions.R")
  
  source("../Functions/formatters.R")
  
  population <- format_population_data("../../Data/population-germany.xlsx",
                                       save=FALSE)
  cancer.data =  format_cancer_data("../../Data/lungCancer-germany.xls", 
                             population=population, save=FALSE)
  
  data = list(obs = cancer.data,
              nx = max(cancer.data$x) + 1,
              config_name = config)
  
  print("nx: ")
  print(data$nx)
  
  stan_fit <- run_stan_program_lcc(data, chains=chains, warmup=warmup, iter=iter, stan_program=stan_program)
  
  store_stan_results(fit=stan_fit, output.path=output.path, config=config)
}

run_stan_v18_3(stan_program="../Synthetic\ data/Stan\ analyses/stan_programs/stan_analysis_cohort_rw2.stan", chains=4, warmup = 50, iter = 100, markov=TRUE)
