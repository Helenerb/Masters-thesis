# base script containing functions running STAN analysis

save.figure <- function(plot, name, path){
  #'plot: <gg object>
  #'name: <string>  on the format '<name>.png'
  #'
  ggsave(paste(name, '.png', sep=""),
         plot = plot,
         device = "png",
         path = path,
         height = 5, width = 8, 
         dpi = "retina"
  )
  ggsave(paste(name, '.pdf', sep=""),
         plot = plot,
         device = "pdf",
         path = path,
         height = 5, width = 8, 
         dpi = "retina"
  )
}

run_stan_analysis <- function(stan_program, chains = 4, warmup = 1000,
                              iter = 10000, markov = TRUE, cohort = FALSE){
  #'Runs STAN to perform Bayesian inference
  #'
  #'A collective function for setting up the correct environment and running the specified
  #'STAN program. 
  #'
  #'@param stan_program The name of the file containing a full STAN program. On the format '<name of file>.stan'
  #'@param data Named list: The input data for the stan program and the underlying effects used to produce the data
  #'@param chains Number of chains in stan program
  #'@param warmup Number of warm up iterations
  #'@param iter Number of iterations
  #'@param markov If program will be run on Markov server
  #'@param cohort IF a cohort effect is included in the model
  #'
  if(markov){
    .libPaths("~/Documents/R_libraries")
    setwd("~/Documents/GitHub/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data")
    output.path <- '~/Documents/GitHub/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Output'
  } else {
    setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data")
    output.path <- '~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Output'
  }
  
  source("configurations_synthetic_data.R")
  
  data = configuration.v10()
  
  obs <- data$obs
  
  library("rstan")
  
  input_stan.lc <- list(
    X=length(unique(obs$x)),
    T=length(unique(obs$t)),
    x=(obs$x + 1),
    t=(obs$t + 1),
    ts = obs$t,
    E = obs$E,
    Y = obs$Y,
    nx = length(unique(obs$x)),
    nt = length(unique(obs$t))
  )
  
  run_info <- list(
    stan.file = "stan_analysis_lc_v4.stan",
    chains = chains,
    warmup = warmup,
    iter = iter,
    configuration = data$config_name,
    type_of_prior = "loggamma, less informative on kappa"
  )
  
  fit <- stan(
    file="stan_analysis_lc_v4.stan",
    data = input_stan.lc,
    chains=chains,
    warmup = warmup,
    iter = iter,
    refresh = 1000,
    seed=123
  )
  
  fit_summary.stan.lc <- summary(fit)
  stan_lc_df <- as.data.frame(fit_summary.stan.lc$summary)  # TODO: Save!
  
  results.path <- file.path(output.path, 'Data')
  
  save(stan_lc_df, file=file.path(results.path, paste('stan_', data$config_name, '_.Rda', sep = "")))
  
  figures.path <- file.path(output.path, 'Figures',data$config_name)
  
  # Saving trace plots etc:
  trace.eta <- traceplot(fit, pars=c("eta[1]", "eta[2]", "eta[3]", "eta[4]", "eta[5]",
                                     "eta[6]", "eta[7]", "eta[8]", "eta[9]", "eta[10]" ))
  save.figure(trace.eta,'trace_eta', figures.path)
  
  save.figure(traceplot(fit, pars=c("kappa[1]", "kappa[2]", "kappa[3]", "kappa[4]",
                                              "kappa[5]", "kappa[6]", "kappa[7]", "kappa[8]",
                                              "kappa[9]", "kappa[10]" )),
                        'trace_kappa', figures.path)
  
  save.figure(traceplot(fit, pars=c("beta[1]", "beta[2]", "beta[3]", "beta[4]",
                                              "beta[5]", "beta[6]", "beta[7]", "beta[8]",
                                              "beta[9]", "beta[10]" )),
                        'trace_beta', figures.path)
  
  
  save.figure(traceplot(fit, pars=c("alpha[1]", "alpha[3]", "alpha[11]", "alpha[13]",
                                              "alpha[5]", "alpha[15]", "alpha[7]", "alpha[8]",
                                              "alpha[12]", "alpha[10]" )),
                        'trace_alpha', figures.path)
  
  save.figure(traceplot(fit, pars=c("alpha_raw[1]", "alpha_raw[3]", "alpha_raw[11]", "alpha_raw[13]",
                                              "alpha_raw[5]", "alpha_raw[15]", "alpha_raw[7]", "alpha_raw[8]",
                                              "alpha_raw[12]", "alpha_raw[10]" )),
                        'trace_alpha', figures.path)
  
  save.figure(traceplot(fit, pars=c("beta_raw[1]", "beta_raw[3]", "beta_raw[11]", "beta_raw[13]",
                                              "beta_raw[5]", "beta_raw[15]", "beta_raw[7]", "beta_raw[8]",
                                              "beta_raw[12]", "beta_raw[10]" )),
                        'trace_beta', figures.path)
  
  save.figure(traceplot(fit), 'general_trace', figures.path)
  
}

run_stan_analysis(stan_program="stan_analysis_lc_v4.stan", chains=2, warmup=100, iter= 200, markov = FALSE)
