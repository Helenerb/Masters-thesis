# functions related to saving and running stan programs

library("ggplot2")

save.figure <- function(plot, name, path){
  #'plot: <gg object>
  #'name: <string>  on the format '<name>.png'
  #'
  
  # ggsave(paste(name, '.png', sep=""),
  #        plot = plot,
  #        device = "png",
  #        path = path,
  #        height = 5, width = 8, 
  #        dpi = "retina"
  # )
  ggsave(paste(name, '.pdf', sep=""),
         plot = plot,
         device = "pdf",
         path = path,
         height = 5, width = 8, 
         dpi = "retina"
  )
}

run_stan_program_lc <- function(data, chains, warmup, iter, stan_program="stan_analysis_lc_v4.stan"){
  
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
  
  print("We input the following to Stan: ")
  print(input_stan.lc)
  
  fit <- stan(
    file=stan_program,
    data = input_stan.lc,
    chains=chains,
    warmup = warmup,
    iter = iter,
    refresh = iter%/%10,
    seed=123
  )
  
  return(fit)
}

run_stan_program_lcc <- function(data, chains, warmup, iter, stan_program="stan_analysis_lcc_phi_kappa.stan"){
  
  obs <- data$obs
  
  library("rstan")
  
  input_stan.lc <- list(
    X=length(unique(obs$x)),
    T=length(unique(obs$t)),
    C=length(unique(obs$c)),
    x=(obs$x + 1),
    t=(obs$t + 1),
    c=(obs$c + data$nx),
    ts = obs$t,
    E = obs$E,
    Y = obs$Y,
    nx = length(unique(obs$x)),
    nt = length(unique(obs$t))
  )
  
  run_info <- list(
    stan.file = stan_program,
    chains = chains,
    warmup = warmup,
    iter = iter,
    configuration = data$config_name
  )
  
  fit <- stan(
    file=stan_program,
    data = input_stan.lc,
    chains=chains,
    warmup = warmup,
    iter = iter,
    refresh = iter%/%10,
    seed=123
  )
  
  return(fit)
}

store_stan_results <- function(fit, output.path, config, stan_program = "", chains = "", warmup = "", iter = "", cohort=TRUE){
  
  # save full stanfit object
  # tryCatch({
  #   saveRDS(fit, file.path(output.path, 'stan_fit.rds'))
  # },
  # error = function(cond){
  #   message("Could not save full stan fit: ")
  #   message(cond)
  # },
  # warning = function(cond){
  #   message("Warning: ")
  #   message(cond)
  # })
  
  list_of_draws <- rstan::extract(fit)
  
  figures.path <- output.path
  
  tau_alpha_draws <- list_of_draws$tau_alpha
  save.figure(qplot(data.frame(x = tau_alpha_draws)$x, geom="histogram"), "tau_alpha_draws", figures.path)
  
  tau_beta_draws <- list_of_draws$tau_beta
  save.figure(qplot(data.frame(x = tau_beta_draws)$x, geom="histogram"), "tau_beta_draws", figures.path)
  
  tau_kappa_draws <- list_of_draws$tau_kappa
  save.figure(qplot(data.frame(x = tau_kappa_draws)$x, geom="histogram"), "tau_kappa_draws", figures.path)
  
  if(cohort){
    tau_gamma_draws <- list_of_draws$tau_gamma
    save.figure(qplot(data.frame(x = tau_gamma_draws)$x, geom="histogram"), "tau_gamma_draws", figures.path)
  }
  
  tau_epsilon_draws <- list_of_draws$tau_epsilon
  save.figure(qplot(data.frame(x = tau_epsilon_draws)$x, geom="histogram"), "tau_epsilon_draws", figures.path)
  
  intercept_draws <- list_of_draws$intercept
  save.figure(qplot(data.frame(x = intercept_draws)$x, geom="histogram"), "intercept_draws", figures.path)
  
  alpha_draws <- list_of_draws$alpha
  beta_draws <- list_of_draws$beta
  kappa_draws <- list_of_draws$kappa
  if(cohort){
    gamma_draws <- list_of_draws$gamma
  }
  eta_draws <- list_of_draws$eta
  
  # save info about run to txt file
  run_info <- list(
    stan.file = stan_program,
    chains = chains,
    warmup = warmup,
    iter = iter,
    configuration = config
  )
  lapply(run_info, write, file.path(output.path, 'run_info.txt'), append = TRUE, ncolumns = 1000)
  
  elapsed_time <- get_elapsed_time(fit)
  write.table(elapsed_time, file=file.path(output.path, 'elapsed_time.txt'))
  
  fit_summary.stan.lc <- summary(fit)
  stan_lc_df <- as.data.frame(fit_summary.stan.lc$summary)
  
  results.path <- output.path
  
  save(stan_lc_df, file=file.path(results.path, paste('stan_', config, '.Rda', sep = "")))
  
  tryCatch({
    save(tau_alpha_draws, file = file.path(output.path, 'draws_tau_alpha.RData'))
    save(alpha_draws, file = file.path(output.path, 'draws_alpha.RData'))
    save(tau_beta_draws, file = file.path(output.path, 'draws_tau_beta.RData'))
    save(beta_draws, file = file.path(output.path, 'draws_beta.RData'))
    save(tau_kappa_draws, file = file.path(output.path, 'draws_tau_kappa.RData'))
    save(kappa_draws, file = file.path(output.path, 'draws_kappa.RData'))
    if(cohort){
      save(tau_gamma_draws, file = file.path(output.path, 'draws_tau_gamma.RData'))
      save(gamma_draws, file = file.path(output.path, 'draws_gamma.RData'))
    }
    save(tau_epsilon_draws, file = file.path(output.path, 'draws_tau_epsilon.RData'))
    save(intercept_draws, file = file.path(output.path, 'draws_intercept.RData'))
    save(eta_draws, file = file.path(output.path, 'draws_eta.RData'))
  },
  error = function(cond){
    message("Could not save lists of marginals \n")
    message(cond)
    message(" ")
  },
  warning = function(cond){
    message("Trying to save list of marginals gives the following warning \n")
    message(cond)
    message(" ")
  })
  
  # Saving trace plots etc:
  tryCatch({
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
                'trace_alpha_raw', figures.path)
    
    save.figure(traceplot(fit, pars=c("beta_raw[1]", "beta_raw[3]", "beta_raw[11]", "beta_raw[13]",
                                      "beta_raw[5]", "beta_raw[15]", "beta_raw[7]", "beta_raw[8]",
                                      "beta_raw[12]", "beta_raw[10]" )),
                'trace_beta_raw', figures.path)
    if(cohort){
      save.figure(traceplot(fit,pars = c("tau_alpha", "tau_beta", "tau_kappa",
                                         "tau_gamma", "tau_epsilon",
                                         "intercept")), 'trace_hyperpars',
                  figures.path)
    } else {
      save.figure(traceplot(fit,pars = c("tau_alpha", "tau_beta", "tau_kappa",
                                         "tau_epsilon",
                                         "intercept")), 'trace_hyperpars',
                  figures.path)
    }
  },
  error = function(cond){
    message("Error in saving trace plot: ")
    message(cond)
    message(" ")
  },
  warning = function(cond){
    message("Warning when saving trace plot: ")
    message(cond)
    message(" ")
  })
  
}