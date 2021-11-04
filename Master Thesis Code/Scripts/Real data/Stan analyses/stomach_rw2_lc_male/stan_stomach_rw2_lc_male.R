# runs stan analysis for stomach cancer data, adapted to be runable
# on Markov. 

# this version runs the lc-model 

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

run_stan_stomach_rw2_lc_male <- function(stan_program, chains=4, warmup=1000, iter=10000, markov=TRUE){
  config = "stomach_rw2_lc_male"
  output.path <- set_workspace(config=config, markov)
  source("../Synthetic\ data/run_stan_functions.R")
  
  source("../Functions/formatters.R")
  
  population <- format_population_data("../../Data/population-germany.xlsx",
                                       save=FALSE)
  cancer.data =  format_cancer_data("../../Data/stomachCancer-germany.xls", 
                                    population=population, save=FALSE)
  
  cancer.data <- cancer.data %>%
    mutate(female.mr = female/female.t, male.mr = male/male.t)
  
  cancer.male <- cancer.data %>%
    select(x, x.c, t, xt, cohort, c, age, age.int, year, birth.year, male.t, male, male.mr) %>%
    mutate(E = male.t, Y = male, `mortality rate` = male.mr)
  
  # Hack: set nx=1, since the real-data cohort indices are not negative. 
  data = list(obs = cancer.male,
              nx = 1,
              config_name = config)
  
  stan_fit <- run_stan_program_lc(data, chains=chains, warmup=warmup, iter=iter, stan_program=stan_program)
  
  store_stan_results(fit=stan_fit, output.path=output.path, config=config,
                     chains=chains, warmup=warmup, iter=iter,
                     stan_program=stan_program, cohort=FALSE)
  
}

run_stan_stomach_rw2_lc_male(stan_program="../Synthetic\ data/Stan\ analyses/stan_programs/stan_analysis_lc_rw2.stan", chains=4, warmup = 50, iter = 100, markov=TRUE)
