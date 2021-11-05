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

run_stan_stomach_rw2_lc_female <- function(stan_program, chains=4, warmup=1000, iter=10000, markov=TRUE){
  config = "stomach_rw2_lc_female"
  output.path <- set_workspace(config=config, markov)
  source("../Synthetic\ data/run_stan_functions.R")
  
  source("../Functions/formatters.R")
  
  population <- format_population_data("../../Data/population-germany.xlsx",
                                       save=FALSE)
  cancer.data =  format_cancer_data("../../Data/stomachCancer-germany.xls", 
                                    population=population, save=FALSE)
  
  cancer.data <- cancer.data %>%
    mutate(female.mr = female/female.t, male.mr = male/male.t)
  
  cancer.female <- cancer.data %>%
    select(x, x.c, t, xt, cohort, c, age, age.int, year, birth.year, female.t, female, female.mr) %>%
    mutate(E = female.t, Y = female, `mortality rate` = female.mr)
  
  # Hack: set nx=1, since the real-data cohort indices are not negative. 
  data = list(obs = cancer.female,
              nx = 1,
              config_name = config)
  
  stan_fit <- run_stan_program_lc(data, chains=chains, warmup=warmup, iter=iter, stan_program=stan_program)
  
  store_stan_results(fit=stan_fit, output.path=output.path, config=config,
                     chains=chains, warmup=warmup, iter=iter,
                     stan_program=stan_program, cohort=FALSE)
  
}

run_stan_stomach_rw2_lc_female(stan_program="../Synthetic\ data/Stan\ analyses/stan_programs/stan_analysis_lc_rw2.stan", chains=6, warmup = 20000, iter = 200000, markov=TRUE)
