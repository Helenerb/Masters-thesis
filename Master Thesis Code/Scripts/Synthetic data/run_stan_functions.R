# functions related to saving and running stan programs

library("ggplot2")

save.figure <- function(plot, name, path, pdf=TRUE, png=TRUE){
  #'@param plot: <gg object>
  #'@param name: <string>  on the format '<name>'
  #'@param path (string): path to where figures should be stored
  
  if(png){
    ggsave(paste(name, '.png', sep=""),
           plot = plot,
           device = "png",
           path = path,
           height = 5, width = 8,
           dpi = "retina"
    )
  }
  
  if(pdf){
    ggsave(paste(name, '.pdf', sep=""),
           plot = plot,
           device = "pdf",
           path = path,
           height = 5, width = 8, 
           dpi = "retina"
    )
  }
}

run_stan_program_lc <- function(data, chains, warmup, iter, stan_program){
  
  obs <- data$obs
  
  library("rstan")
  
  input_stan.lc <- list(
    X=length(unique(obs$x)),
    T=length(unique(obs$t)),
    x=(obs$x + 1),
    t=(obs$t + 1),
    #ts = obs$t,
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

run_stan_program_traditional_lc <- function(data, chains, warmup, iter, stan_program="stan_lc_traditional.stan"){
  # Deprecate! 
  obs.trad <- data$obs
  
  library("rstan")
  
  input_stan.lc <- list(
    X=length(unique(obs.trad$x)),
    T=length(unique(obs.trad$t)),
    x=(obs.trad$x + 1),
    t=(obs.trad$t + 1),
    #mr = obs.trad$mr_gaussian,
    exp_mr = obs.trad$eta,  # should be named log_mr
    nx = length(unique(obs.trad$x)),
    nt = length(unique(obs.trad$t))
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

run_stan_program_gaussian <- function(data, chains, warmup, iter, stan_program){
  obs <- data$obs
  
  library("rstan")
  
  input_stan.lc <- list(
    X=length(unique(obs$x)),
    T=length(unique(obs$t)),
    x=(obs$x + 1),
    t=(obs$t + 1),
    log_mr = obs$eta,
    nx = length(unique(obs$x)),
    nt = length(unique(obs$t))
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

store_stan_results_traditional <- function(fit, output.path, config, stan_program = "", chains = "", warmup = "", iter = "", cohort=F){
  figures.path <- output.path
  results.path <- output.path
  
  list_of_draws <- rstan::extract(fit)
  
  # extract draws for hyperparameters, and save histograms
  tau_alpha_draws <- list_of_draws$tau_alpha
  save.figure(qplot(data.frame(x = tau_alpha_draws)$x, geom="histogram"), "tau_alpha_draws", figures.path, pdf=F)
  
  tau_beta_draws <- list_of_draws$tau_beta
  save.figure(qplot(data.frame(x = tau_beta_draws)$x, geom="histogram"), "tau_beta_draws", figures.path, pdf=F)
  
  tau_kappa_draws <- list_of_draws$tau_kappa
  save.figure(qplot(data.frame(x = tau_kappa_draws)$x, geom="histogram"), "tau_kappa_draws", figures.path, pdf=F)
  
  if(cohort){
    tau_gamma_draws <- list_of_draws$tau_gamma
    save.figure(qplot(data.frame(x = tau_gamma_draws)$x, geom="histogram"), "tau_gamma_draws", figures.path, pdf=F)
  }
  
  tau_epsilon_draws <- list_of_draws$tau_epsilon
  save.figure(qplot(data.frame(x = tau_epsilon_draws)$x, geom="histogram"), "tau_epsilon_draws", figures.path, pdf=F)
  
  intercept_draws <- list_of_draws$intercept
  save.figure(qplot(data.frame(x = intercept_draws)$x, geom="histogram"), "intercept_draws", figures.path, pdf=F)
  
  # extract draws for random effects
  alpha_draws <- list_of_draws$alpha
  beta_draws <- list_of_draws$beta
  kappa_draws <- list_of_draws$kappa
  if(cohort){
    gamma_draws <- list_of_draws$gamma
  }
  
  # extract draws from predictor
  eta_draws <- list_of_draws$eta
  eta_draws_reduced <- eta_draws[seq(1,nrow(eta_draws),5), ]
  eta_draws_100 <- eta_draws[seq(1,nrow(eta_draws), 100), ]
  
  
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
  
  # extract summary of fit and save to dataframe
  fit_summary.stan.lc <- summary(fit)
  stan_lc_df <- as.data.frame(fit_summary.stan.lc$summary)
  save(stan_lc_df, file=file.path(results.path, paste('stan_', config, '.Rda', sep = "")))
  
  tryCatch({
    # save STAN draws of parameters
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
    save(eta_draws_100, file = file.path(output.path, "draws_eta_100.RData"))
    save(eta_draws_reduced, file = file.path(output.path, "draws_eta_reduced.RData"))
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
  },
  finally = {
    
    params <- row.names(stan_lc_df)
    alphas <- str_extract(params, "alpha\\[([0-9+])\\]"); alphas <- alphas[!is.na(alphas)]
    betas <- str_extract(params, "beta\\[([0-9+])\\]"); betas <- betas[!is.na(betas)]
    kappas <- str_extract(params, "kappa\\[([0-9+])\\]"); kappas <- kappas[!is.na(kappas)]
    
    tryCatch({
      
      save.figure(traceplot(fit, pars=kappas),
                  'trace_kappa', figures.path, png=F)
      
      save.figure(traceplot(fit, pars=betas),
                  'trace_beta', figures.path, png=F)
      
      
      save.figure(traceplot(fit, pars=alphas),
                  'trace_alpha', figures.path, png=F)
      
      if(cohort){
        save.figure(traceplot(fit,pars = c("tau_alpha", "tau_beta", "tau_kappa",
                                           "tau_gamma", "tau_epsilon",
                                           "intercept")), 'trace_hyperpars',
                    figures.path, png=F)
      } else {
        save.figure(traceplot(fit,pars = c("tau_alpha", "tau_beta", "tau_kappa",
                                           "tau_epsilon", "intercept")), 'trace_hyperpars',
                    figures.path, png=F)
      }
      
      trace.eta <- traceplot(fit, pars=c("eta[1]", "eta[2]", "eta[3]", "eta[4]", "eta[5]",
                                         "eta[6]", "eta[7]", "eta[8]", "eta[9]", "eta[10]" ))
      save.figure(trace.eta,'trace_eta', figures.path, png=F)
      
    },
    error = function(cond){
      message("Error in saving trace plot as pdf: ")
      message(cond)
      message(" ")
    },
    warning = function(cond){
      message("Warning when saving trace plot as pdf: ")
      message(cond)
      message(" ")
    },
    finally = {
      tryCatch(
        {
          
          save.figure(traceplot(fit, pars=kappas),
                      'trace_kappa', figures.path, pdf=F)
          
          save.figure(traceplot(fit, pars=betas),
                      'trace_beta', figures.path, pdf=F)
          
          
          save.figure(traceplot(fit, pars=alphas),
                      'trace_alpha', figures.path, pdf=F)
          
          if(cohort){
            save.figure(traceplot(fit,pars = c("tau_alpha", "tau_beta", "tau_kappa",
                                               "tau_gamma", "tau_epsilon",
                                               "intercept")), 'trace_hyperpars',
                        figures.path, pdf=F)
          } else {
            save.figure(traceplot(fit,pars = c("tau_alpha", "tau_beta", "tau_kappa",
                                               "tau_epsilon", "intercept")), 'trace_hyperpars',
                        figures.path, pdf=F)
          }
          
          trace.eta <- traceplot(fit, pars=c("eta[1]", "eta[2]", "eta[3]", "eta[4]", "eta[5]",
                                             "eta[6]", "eta[7]", "eta[8]", "eta[9]", "eta[10]" ))
          save.figure(trace.eta,'trace_eta', figures.path, pdf=F)
        },
        error = function(cond){
          message("Error in saving trace plot as png: ")
          message(cond)
          message(" ")
        },
        warning = function(cond){
          message("Warning in saving trace plot as png: ")
          message(cond)
          message(" ")
        }
      )
    })
  })
}

store_stan_results_gaus_linear <- function(fit, output.path, config, stan_program = "", chains = "", warmup = "", iter = ""){
  figures.path <- output.path
  results.path <- output.path
  
  list_of_draws <- rstan::extract(fit)
  
  # extract draws for hyperparameters, and save histograms
  tau_alpha_draws <- list_of_draws$tau_alpha
  save.figure(qplot(data.frame(x = tau_alpha_draws)$x, geom="histogram"), "tau_alpha_draws", figures.path, pdf=F)
  
  tau_beta_draws <- list_of_draws$tau_beta
  save.figure(qplot(data.frame(x = tau_beta_draws)$x, geom="histogram"), "tau_beta_draws", figures.path, pdf=F)
  
  tau_kappa_draws <- list_of_draws$tau_kappa
  save.figure(qplot(data.frame(x = tau_kappa_draws)$x, geom="histogram"), "tau_kappa_draws", figures.path, pdf=F)
  
  tau_epsilon_draws <- list_of_draws$tau_epsilon
  save.figure(qplot(data.frame(x = tau_epsilon_draws)$x, geom="histogram"), "tau_epsilon_draws", figures.path, pdf=F)
  
  intercept_draws <- list_of_draws$intercept
  save.figure(qplot(data.frame(x = intercept_draws)$x, geom="histogram"), "intercept_draws", figures.path, pdf=F)
  
  # extract draws for random effects
  alpha_draws <- list_of_draws$alpha
  beta_draws <- list_of_draws$beta
  kappa_draws <- list_of_draws$kappa
  
  # extract draws from predictor
  eta_draws <- list_of_draws$eta
  eta_draws_reduced <- eta_draws[seq(1,nrow(eta_draws),5), ]
  eta_draws_100 <- eta_draws[seq(1,nrow(eta_draws), 100), ]
  
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
  
  # extract summary of fit and save to dataframe
  fit_summary.stan.lc <- summary(fit)
  stan_lc_df <- as.data.frame(fit_summary.stan.lc$summary)
  save(stan_lc_df, file=file.path(results.path, paste('stan_', config, '.Rda', sep = "")))
  
  tryCatch({
    # save STAN draws of parameters
    save(tau_alpha_draws, file = file.path(output.path, 'draws_tau_alpha.RData'))
    save(alpha_draws, file = file.path(output.path, 'draws_alpha.RData'))
    save(tau_kappa_draws, file = file.path(output.path, 'draws_tau_kappa.RData'))
    save(kappa_draws, file = file.path(output.path, 'draws_kappa.RData'))
    save(tau_epsilon_draws, file = file.path(output.path, 'draws_tau_epsilon.RData'))
    save(intercept_draws, file = file.path(output.path, 'draws_intercept.RData'))
    save(eta_draws_100, file = file.path(output.path, "draws_eta_100.RData"))
    save(eta_draws_reduced, file = file.path(output.path, "draws_eta_reduced.RData"))
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
  },
  finally = {
    
    params <- row.names(stan_lc_df)
    alphas <- str_extract(params, "alpha\\[([0-9+])\\]"); alphas <- alphas[!is.na(alphas)]
    kappas <- str_extract(params, "kappa\\[([0-9+])\\]"); kappas <- kappas[!is.na(kappas)]
    
    tryCatch({
      
      save.figure(traceplot(fit, pars=kappas),
                  'trace_kappa', figures.path, png=F)
      
      
      save.figure(traceplot(fit, pars=alphas),
                  'trace_alpha', figures.path, png=F)
      
      
      save.figure(traceplot(
        fit,pars = c("tau_alpha", "tau_kappa","tau_epsilon", "intercept")),
        'trace_hyperpars', figures.path, png=F)
      
      trace.eta <- traceplot(fit, pars=c("eta[1]", "eta[2]", "eta[3]", "eta[4]", "eta[5]",
                                         "eta[6]", "eta[7]", "eta[8]", "eta[9]", "eta[10]" ))
      save.figure(trace.eta,'trace_eta', figures.path, png=F)
      
    },
    error = function(cond){
      message("Error in saving trace plot as pdf: ")
      message(cond)
      message(" ")
    },
    warning = function(cond){
      message("Warning when saving trace plot as pdf: ")
      message(cond)
      message(" ")
    },
    finally = {
      tryCatch(
        {
          save.figure(traceplot(fit, pars=kappas), 'trace_kappa', figures.path, pdf=F)
          
          save.figure(traceplot(fit, pars=alphas), 'trace_alpha', figures.path, pdf=F)
          
            save.figure(traceplot(
              fit, pars = c("tau_alpha", "tau_kappa", "tau_epsilon", "intercept")),
              'trace_hyperpars', figures.path, pdf=F)
          
          trace.eta <- traceplot(fit, pars=c("eta[1]", "eta[2]", "eta[3]", "eta[4]", "eta[5]",
                                             "eta[6]", "eta[7]", "eta[8]", "eta[9]", "eta[10]" ))
          save.figure(trace.eta,'trace_eta', figures.path, pdf=F)
        },
        error = function(cond){
          message("Error in saving trace plot as png: ")
          message(cond)
          message(" ")
        },
        warning = function(cond){
          message("Warning in saving trace plot as png: ")
          message(cond)
          message(" ")
        }
      )
    })
  })
}

store_stan_results_traditional.extra.error.term <- function(fit, output.path, config, stan_program = "", chains = "", warmup = "", iter = "", cohort=F){
  figures.path <- output.path
  results.path <- output.path
  
  list_of_draws <- rstan::extract(fit)
  
  # extract draws for hyperparameters, and save histograms
  tau_alpha_draws <- list_of_draws$tau_alpha
  save.figure(qplot(data.frame(x = tau_alpha_draws)$x, geom="histogram"), "tau_alpha_draws", figures.path, pdf=F)
  
  tau_beta_draws <- list_of_draws$tau_beta
  save.figure(qplot(data.frame(x = tau_beta_draws)$x, geom="histogram"), "tau_beta_draws", figures.path, pdf=F)
  
  tau_kappa_draws <- list_of_draws$tau_kappa
  save.figure(qplot(data.frame(x = tau_kappa_draws)$x, geom="histogram"), "tau_kappa_draws", figures.path, pdf=F)
  
  if(cohort){
    tau_gamma_draws <- list_of_draws$tau_gamma
    save.figure(qplot(data.frame(x = tau_gamma_draws)$x, geom="histogram"), "tau_gamma_draws", figures.path, pdf=F)
  }
  
  tau_epsilon_draws <- list_of_draws$tau_epsilon
  save.figure(qplot(data.frame(x = tau_epsilon_draws)$x, geom="histogram"), "tau_epsilon_draws", figures.path, pdf=F)
  
  intercept_draws <- list_of_draws$intercept
  save.figure(qplot(data.frame(x = intercept_draws)$x, geom="histogram"), "intercept_draws", figures.path, pdf=F)
  
  # extract draws for random effects
  alpha_draws <- list_of_draws$alpha
  beta_draws <- list_of_draws$beta
  kappa_draws <- list_of_draws$kappa
  if(cohort){
    gamma_draws <- list_of_draws$gamma
  }
  
  # extract draws from predictor
  eta_draws <- list_of_draws$eta
  eta_draws_reduced <- eta_draws[seq(1,nrow(eta_draws),5), ]
  eta_draws_100 <- eta_draws[seq(1,nrow(eta_draws), 100), ]
  
  # extract draws from error term:
  error_term_draws <- list_of_draws$error_term
  tau_error_term_draws <- list_of_draws$tau_error_term
  save.figure(qplot(data.frame(x = tau_error_term_draws)$x, geom="histogram"), "tau_error_term_draws", figures.path, pdf=F)
  
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
  
  # extract summary of fit and save to dataframe
  fit_summary.stan.lc <- summary(fit)
  stan_lc_df <- as.data.frame(fit_summary.stan.lc$summary)
  save(stan_lc_df, file=file.path(results.path, paste('stan_', config, '.Rda', sep = "")))
  
  tryCatch({
    # save STAN draws of parameters
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
    save(eta_draws_100, file = file.path(output.path, "draws_eta_100.RData"))
    save(eta_draws_reduced, file = file.path(output.path, "draws_eta_reduced.RData"))
    save(eta_draws, file = file.path(output.path, 'draws_eta.RData'))
    save(error_term_draws, file = file.path(output.path, 'draws_error_term.RData'))
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
  },
  finally = {
    
    params <- row.names(stan_lc_df)
    alphas <- str_extract(params, "alpha\\[([0-9+])\\]"); alphas <- alphas[!is.na(alphas)]
    betas <- str_extract(params, "beta\\[([0-9+])\\]"); betas <- betas[!is.na(betas)]
    kappas <- str_extract(params, "kappa\\[([0-9+])\\]"); kappas <- kappas[!is.na(kappas)]
    error_terms <- str_extract(params, "error_term\\[([0-9+])\\]"); error_terms <- error_terms[!is.na(error_terms)]
    
    tryCatch({
      
      save.figure(traceplot(fit, pars=kappas),
                  'trace_kappa', figures.path, png=F)
      
      save.figure(traceplot(fit, pars=betas),
                  'trace_beta', figures.path, png=F)
      
      
      save.figure(traceplot(fit, pars=alphas),
                  'trace_alpha', figures.path, png=F)
      
      if(cohort){
        save.figure(traceplot(fit,pars = c("tau_alpha", "tau_beta", "tau_kappa",
                                           "tau_gamma", "tau_epsilon",
                                           "intercept")), 'trace_hyperpars',
                    figures.path, png=F)
      } else {
        save.figure(traceplot(fit,pars = c("tau_alpha", "tau_beta", "tau_kappa",
                                           "tau_epsilon", "intercept")), 'trace_hyperpars',
                    figures.path, png=F)
      }
      
      trace.eta <- traceplot(fit, pars=c("eta[1]", "eta[2]", "eta[3]", "eta[4]", "eta[5]",
                                         "eta[6]", "eta[7]", "eta[8]", "eta[9]", "eta[10]" ))
      save.figure(trace.eta,'trace_eta', figures.path, png=F)
      
      save.figure(traceplot(fit, pars=error_terms),
                  'trace_error_term', figures.path, png=F)
      
    },
    error = function(cond){
      message("Error in saving trace plot as pdf: ")
      message(cond)
      message(" ")
    },
    warning = function(cond){
      message("Warning when saving trace plot as pdf: ")
      message(cond)
      message(" ")
    },
    finally = {
      tryCatch(
        {
          
          save.figure(traceplot(fit, pars=kappas),
                      'trace_kappa', figures.path, pdf=F)
          
          save.figure(traceplot(fit, pars=betas),
                      'trace_beta', figures.path, pdf=F)
          
          
          save.figure(traceplot(fit, pars=alphas),
                      'trace_alpha', figures.path, pdf=F)
          
          if(cohort){
            save.figure(traceplot(fit,pars = c("tau_alpha", "tau_beta", "tau_kappa",
                                               "tau_gamma", "tau_epsilon",
                                               "intercept")), 'trace_hyperpars',
                        figures.path, pdf=F)
          } else {
            save.figure(traceplot(fit,pars = c("tau_alpha", "tau_beta", "tau_kappa",
                                               "tau_epsilon", "intercept")), 'trace_hyperpars',
                        figures.path, pdf=F)
          }
          
          trace.eta <- traceplot(fit, pars=c("eta[1]", "eta[2]", "eta[3]", "eta[4]", "eta[5]",
                                             "eta[6]", "eta[7]", "eta[8]", "eta[9]", "eta[10]" ))
          save.figure(trace.eta,'trace_eta', figures.path, pdf=F)
          
          save.figure(traceplot(fit, pars=error_terms),
                      'trace_error_term', figures.path, pdf=F)
          
        },
        error = function(cond){
          message("Error in saving trace plot as png: ")
          message(cond)
          message(" ")
        },
        warning = function(cond){
          message("Warning in saving trace plot as png: ")
          message(cond)
          message(" ")
        }
      )
    })
  })
}

store_stan_results <- function(fit, output.path, config, stan_program = "", chains = "", warmup = "", iter = "", cohort=FALSE){
  
  figures.path <- output.path
  results.path <- output.path
  
  list_of_draws <- rstan::extract(fit)
  
  # extract draws for hyperparameters, and save histograms
  tau_alpha_draws <- list_of_draws$tau_alpha
  save.figure(qplot(data.frame(x = tau_alpha_draws)$x, geom="histogram"), "tau_alpha_draws", figures.path, pdf=F)
  
  tau_beta_draws <- list_of_draws$tau_beta
  save.figure(qplot(data.frame(x = tau_beta_draws)$x, geom="histogram"), "tau_beta_draws", figures.path, pdf=F)
  
  tau_kappa_draws <- list_of_draws$tau_kappa
  save.figure(qplot(data.frame(x = tau_kappa_draws)$x, geom="histogram"), "tau_kappa_draws", figures.path, pdf=F)
  
  if(cohort){
    tau_gamma_draws <- list_of_draws$tau_gamma
    save.figure(qplot(data.frame(x = tau_gamma_draws)$x, geom="histogram"), "tau_gamma_draws", figures.path, pdf=F)
  }
  
  tau_epsilon_draws <- list_of_draws$tau_epsilon
  save.figure(qplot(data.frame(x = tau_epsilon_draws)$x, geom="histogram"), "tau_epsilon_draws", figures.path, pdf=F)
  
  intercept_draws <- list_of_draws$intercept
  save.figure(qplot(data.frame(x = intercept_draws)$x, geom="histogram"), "intercept_draws", figures.path, pdf=F)
  
  # extract draws for random effects
  alpha_draws <- list_of_draws$alpha
  beta_draws <- list_of_draws$beta
  kappa_draws <- list_of_draws$kappa
  if(cohort){
    gamma_draws <- list_of_draws$gamma
  }
  epsilon_draws <- list_of_draws$epsilon
  epsilon_draws_100 <- epsilon_draws[seq(1, nrow(epsilon_draws), 100), ]
  
  # extract draws from predictor
  eta_draws <- list_of_draws$eta
  eta_draws_reduced <- eta_draws[seq(1,nrow(eta_draws),5), ]
  eta_draws_100 <- eta_draws[seq(1,nrow(eta_draws), 100), ]
  
  
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
  
  # extract summary of fit and save to dataframe
  fit_summary.stan.lc <- summary(fit)
  stan_lc_df <- as.data.frame(fit_summary.stan.lc$summary)
  save(stan_lc_df, file=file.path(results.path, paste('stan_', config, '.Rda', sep = "")))
  
  tryCatch({
    # save STAN draws of parameters
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
    save(eta_draws_100, file = file.path(output.path, "draws_eta_100.RData"))
    save(epsilon_draws_100, file = file.path(output.path, 'draws_epsilon_100.RData'))
    save(eta_draws_reduced, file = file.path(output.path, "draws_eta_reduced.RData"))
    save(eta_draws, file = file.path(output.path, 'draws_eta.RData'))
    save(epsilon_draws, file = file.path(output.path, 'draws_epsilon.RData'))
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
  },
  finally = {
    
    params <- row.names(stan_lc_df)
    alphas <- str_extract(params, "alpha\\[([0-9+])\\]"); alphas <- alphas[!is.na(alphas)]
    betas <- str_extract(params, "beta\\[([0-9+])\\]"); betas <- betas[!is.na(betas)]
    kappas <- str_extract(params, "kappa\\[([0-9+])\\]"); kappas <- kappas[!is.na(kappas)]
    
    tryCatch({
      
      save.figure(traceplot(fit, pars=kappas),
                  'trace_kappa', figures.path, png=F)
      
      save.figure(traceplot(fit, pars=betas),
                  'trace_beta', figures.path, png=F)
      
      
      save.figure(traceplot(fit, pars=alphas),
                  'trace_alpha', figures.path, png=F)

      if(cohort){
        save.figure(traceplot(fit,pars = c("tau_alpha", "tau_beta", "tau_kappa",
                                           "tau_gamma", "tau_epsilon",
                                           "intercept")), 'trace_hyperpars',
                    figures.path, png=F)
      } else {
        save.figure(traceplot(fit,pars = c("tau_alpha", "tau_beta", "tau_kappa",
                                           "tau_epsilon", "intercept")), 'trace_hyperpars',
                    figures.path, png=F)
      }
      
      trace.eta <- traceplot(fit, pars=c("eta[1]", "eta[2]", "eta[3]", "eta[4]", "eta[5]",
                                         "eta[6]", "eta[7]", "eta[8]", "eta[9]", "eta[10]" ))
      save.figure(trace.eta,'trace_eta', figures.path, png=F)
      
    },
    error = function(cond){
      message("Error in saving trace plot as pdf: ")
      message(cond)
      message(" ")
    },
    warning = function(cond){
      message("Warning when saving trace plot as pdf: ")
      message(cond)
      message(" ")
    },
    finally = {
      tryCatch(
        {
          
          save.figure(traceplot(fit, pars=kappas),
                     'trace_kappa', figures.path, pdf=F)
          
          save.figure(traceplot(fit, pars=betas),
                      'trace_beta', figures.path, pdf=F)
          
          
          save.figure(traceplot(fit, pars=alphas),
                      'trace_alpha', figures.path, pdf=F)
          
          if(cohort){
            save.figure(traceplot(fit,pars = c("tau_alpha", "tau_beta", "tau_kappa",
                                               "tau_gamma", "tau_epsilon",
                                               "intercept")), 'trace_hyperpars',
                        figures.path, pdf=F)
          } else {
            save.figure(traceplot(fit,pars = c("tau_alpha", "tau_beta", "tau_kappa",
                                               "tau_epsilon", "intercept")), 'trace_hyperpars',
                        figures.path, pdf=F)
          }
          
          trace.eta <- traceplot(fit, pars=c("eta[1]", "eta[2]", "eta[3]", "eta[4]", "eta[5]",
                                             "eta[6]", "eta[7]", "eta[8]", "eta[9]", "eta[10]" ))
          save.figure(trace.eta,'trace_eta', figures.path, pdf=F)
          },
        error = function(cond){
          message("Error in saving trace plot as png: ")
          message(cond)
          message(" ")
        },
        warning = function(cond){
          message("Warning in saving trace plot as png: ")
          message(cond)
          message(" ")
        }
      )
    })
  })
  
}

store_stan_results_pois_linear <- function(fit, output.path, config, stan_program = "", chains = "", warmup = "", iter = "", cohort=FALSE){
  
  figures.path <- output.path
  results.path <- output.path
  
  list_of_draws <- rstan::extract(fit)
  
  # extract draws for hyperparameters, and save histograms
  tau_alpha_draws <- list_of_draws$tau_alpha
  save.figure(qplot(data.frame(x = tau_alpha_draws)$x, geom="histogram"), "tau_alpha_draws", figures.path, pdf=F)
  
  tau_kappa_draws <- list_of_draws$tau_kappa
  save.figure(qplot(data.frame(x = tau_kappa_draws)$x, geom="histogram"), "tau_kappa_draws", figures.path, pdf=F)
  
  if(cohort){
    tau_gamma_draws <- list_of_draws$tau_gamma
    save.figure(qplot(data.frame(x = tau_gamma_draws)$x, geom="histogram"), "tau_gamma_draws", figures.path, pdf=F)
  }
  
  tau_epsilon_draws <- list_of_draws$tau_epsilon
  save.figure(qplot(data.frame(x = tau_epsilon_draws)$x, geom="histogram"), "tau_epsilon_draws", figures.path, pdf=F)
  
  intercept_draws <- list_of_draws$intercept
  save.figure(qplot(data.frame(x = intercept_draws)$x, geom="histogram"), "intercept_draws", figures.path, pdf=F)
  
  # extract draws for random effects
  alpha_draws <- list_of_draws$alpha
  kappa_draws <- list_of_draws$kappa
  if(cohort){
    gamma_draws <- list_of_draws$gamma
  }
  epsilon_draws <- list_of_draws$epsilon
  epsilon_draws_100 <- epsilon_draws[seq(1, nrow(epsilon_draws), 100), ]
  
  # extract draws from predictor
  eta_draws <- list_of_draws$eta
  eta_draws_reduced <- eta_draws[seq(1,nrow(eta_draws),5), ]
  eta_draws_100 <- eta_draws[seq(1,nrow(eta_draws), 100), ]
  
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
  
  # extract summary of fit and save to dataframe
  fit_summary.stan.lc <- summary(fit)
  stan_lc_df <- as.data.frame(fit_summary.stan.lc$summary)
  save(stan_lc_df, file=file.path(results.path, paste('stan_', config, '.Rda', sep = "")))
  
  tryCatch({
    # save STAN draws of parameters
    save(tau_alpha_draws, file = file.path(output.path, 'draws_tau_alpha.RData'))
    save(alpha_draws, file = file.path(output.path, 'draws_alpha.RData'))
    save(tau_kappa_draws, file = file.path(output.path, 'draws_tau_kappa.RData'))
    save(kappa_draws, file = file.path(output.path, 'draws_kappa.RData'))
    if(cohort){
      save(tau_gamma_draws, file = file.path(output.path, 'draws_tau_gamma.RData'))
      save(gamma_draws, file = file.path(output.path, 'draws_gamma.RData'))
    }
    save(tau_epsilon_draws, file = file.path(output.path, 'draws_tau_epsilon.RData'))
    save(intercept_draws, file = file.path(output.path, 'draws_intercept.RData'))
    save(eta_draws_100, file = file.path(output.path, "draws_eta_100.RData"))
    save(epsilon_draws_100, file = file.path(output.path, 'draws_epsilon_100.RData'))
    save(eta_draws_reduced, file = file.path(output.path, "draws_eta_reduced.RData"))
    save(eta_draws, file = file.path(output.path, 'draws_eta.RData'))
    save(epsilon_draws, file = file.path(output.path, 'draws_epsilon.RData'))
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
  },
  finally = {
    
    params <- row.names(stan_lc_df)
    alphas <- str_extract(params, "alpha\\[([0-9+])\\]"); alphas <- alphas[!is.na(alphas)]
    kappas <- str_extract(params, "kappa\\[([0-9+])\\]"); kappas <- kappas[!is.na(kappas)]
    
    tryCatch({
      
      save.figure(traceplot(fit, pars=kappas),
                  'trace_kappa', figures.path, png=F)
      
      save.figure(traceplot(fit, pars=alphas),
                  'trace_alpha', figures.path, png=F)
      
      if(cohort){
        save.figure(traceplot(fit,pars = c("tau_alpha", "tau_kappa",
                                           "tau_gamma", "tau_epsilon",
                                           "intercept")), 'trace_hyperpars',
                    figures.path, png=F)
      } else {
        save.figure(traceplot(fit,pars = c("tau_alpha", "tau_kappa",
                                           "tau_epsilon", "intercept")), 'trace_hyperpars',
                    figures.path, png=F)
      }
      
      trace.eta <- traceplot(fit, pars=c("eta[1]", "eta[2]", "eta[3]", "eta[4]", "eta[5]",
                                         "eta[6]", "eta[7]", "eta[8]", "eta[9]", "eta[10]" ))
      save.figure(trace.eta,'trace_eta', figures.path, png=F)
      
    },
    error = function(cond){
      message("Error in saving trace plot as pdf: ")
      message(cond)
      message(" ")
    },
    warning = function(cond){
      message("Warning when saving trace plot as pdf: ")
      message(cond)
      message(" ")
    },
    finally = {
      tryCatch(
        {
          
          save.figure(traceplot(fit, pars=kappas),
                      'trace_kappa', figures.path, pdf=F)
          
          
          save.figure(traceplot(fit, pars=alphas),
                      'trace_alpha', figures.path, pdf=F)
          
          if(cohort){
            save.figure(traceplot(fit,pars = c("tau_alpha", "tau_kappa",
                                               "tau_gamma", "tau_epsilon",
                                               "intercept")), 'trace_hyperpars',
                        figures.path, pdf=F)
          } else {
            save.figure(traceplot(fit,pars = c("tau_alpha", "tau_kappa",
                                               "tau_epsilon", "intercept")), 'trace_hyperpars',
                        figures.path, pdf=F)
          }
          
          trace.eta <- traceplot(fit, pars=c("eta[1]", "eta[2]", "eta[3]", "eta[4]", "eta[5]",
                                             "eta[6]", "eta[7]", "eta[8]", "eta[9]", "eta[10]" ))
          save.figure(trace.eta,'trace_eta', figures.path, pdf=F)
        },
        error = function(cond){
          message("Error in saving trace plot as png: ")
          message(cond)
          message(" ")
        },
        warning = function(cond){
          message("Warning in saving trace plot as png: ")
          message(cond)
          message(" ")
        }
      )
    })
  })
  
}