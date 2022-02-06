#' stand-alone script for 
#' Poisson model, with gamma priors for the precisions, with kappa and alpha as rw1s. 
#' 
#' The name of this configuration (the investigation name) mean:
#' poiss: Poisson model
#' fh: Fixed hyperparameters
#' rw1: kappa (and alpha) modelled as Random walk of order 1

#   ----   Load libraries and set workspace   ----
library("tidyverse")
library("inlabru")
library("ggplot2")
library("INLA")
library("patchwork")
library("rstan")

# TODO: The investigation.name should be the name of the folder where this file is stored, so you should change it to your folder name
# I have it stored in a folder named "stand_alone_poiss_fh_rw1_4".
# Note: the stan results will be stored in a file named <investigation.name>.Rda

investigation.name <- "stand_alone_poiss_gp_rw1_4"  # TODO: change!
investigation.path <- investigation.name

# TODO: Change the first part of this path it to where you have stored your folder containing these files
stan.output  <- file.path("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Step_by_step_results", investigation.path)  # TODO: Change!
output.path <- stan.output

setwd(stan.output)
# Note: at this step, the working directory should be at the location of the folder that contains this file


#   ----    Retrieve the data   ----

synthetic.male.lung.v4 <- function(){
  obs <- read.csv("synthetic_male_lung_4.csv")
  
  obs.trad <- obs %>% 
    select(c(x, t, xt, age.int, year, x.c, alpha, beta, kappa, intercept, epsilon,
             eta, tau.alpha, tau.beta, tau.kappa, tau.epsilon, E)) %>%
    mutate(eta.no.error = intercept + alpha + beta*kappa) %>%
    mutate(mr_gaussian = exp(eta)) %>%
    mutate(Y_gaussian  = mr_gaussian * E)
  
  underlying.effects <- list(obs = obs.trad, nx = 18, nt = 18,
                             alpha.true = {obs %>% filter(t == 0)}$alpha,
                             beta.true = {obs %>% filter(t == 0)}$beta,
                             kappa.true = {obs %>% filter(x == 0)}$kappa,
                             intercept = unique(obs$intercept),
                             age.intercept.true = unique(obs$intercept),
                             tau.alpha.true = unique(obs$tau.alpha),
                             tau.beta.true = unique(obs$tau.beta),
                             tau.kappa.true = unique(obs$tau.kappa),
                             tau.epsilon.true = unique(obs$tau.epsilon))
  
  return(list(obs = obs, underlying.effects = underlying.effects))
}

# We use this data for both inlabru and stan
config.data <- synthetic.male.lung.v4()
obs <- config.data$obs
underlying.effects <- config.data$underlying.effects

#   ----   Run STAN analysis   ----

#source("Scripts/Synthetic\ data/run_stan_functions.R")

#   ----   Functions for running stan   ----
# this actually runs the stan program stan_program with the data in data
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

# this functions just stores the stan results
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
      
      trace.eta <- traceplot(fit, pars=c("eta[1]", "eta[64]", "eta[128]", "eta[192]", "eta[256]",
                                         "eta[54]", "eta[324]", "eta[3]", "eta[9]", "eta[32]" ))
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

# this is a function for saving figures in the same format
save.figure <- function(plot, name, path, pdf=TRUE, png=FALSE){
  #'@param plot gg object
  #'@param name string:  on the format '<name>.png'
  #'@param path string: path to location where figure should be stored
  #'@param pdf (boolean): whether to save as pdf
  #'@param png (boolean) whether to save as png
  
  if (png){
    ggsave(paste(name, '.png', sep=""),
           plot = plot,
           device = "png",
           path = path,
           height = 5, width = 8,
           dpi = "retina"
    )
  }
  if (pdf) {
    ggsave(paste(name, '.pdf', sep=""),
           plot = plot,
           device = "pdf",
           path = path,
           height = 5, width = 8,
           dpi = "retina"
    )
  }
}

# This is a wrapper around run_stan_program_lc, which sorts the data, runs stan and saves the stan results to the correct location
run_stan <- function(stan_program, obs, chains, warmup, iter, output.path, config.name){
  
  stan_fit <- run_stan_program_lc(
    list(obs = obs), chains=chains,warmup=warmup,
    iter=iter, stan_program=stan_program)
  
  store_stan_results(
    fit=stan_fit, output.path=output.path, config=config.name,
    chains=chains, warmup=warmup, iter=iter, stan_program=stan_program)
  
  return(stan_fit)
}

#   ----   Run Stan   ----

# TODO: Say how many chains, warm-up runs and total iterations you want!
stan.fit <- run_stan(
  stan_program="stan_pois_gp_rw1_sc.stan",
  obs = obs, chains=1, warmup = 30, iter = 300, output.path = stan.output,
  config.name = investigation.name)  ## TODO: Change!

#   ----   Functions for running Inlabru   ----

# this function configures the model in inlabru and runs it
inlabru.pois.gp.rw1 <- function(obs, max_iter=30){
  #'Implements inlabru analysis for lc model, fixing the precisions and modelling all random effects as iid
  #'
  #'@param obs: Contains the observed data and the real underlying random effects
  #'@param max_iter (int): maximum number of iterations in inlabru
  
  nx = length(unique(obs$x))
  nt = length(unique(obs$t))
  
  # constraints for the age effect beta
  A.beta = matrix(1, nrow = 1, ncol = nx)
  e.beta = 1
  
  # fixed.theta.alpha <- list(prec = list(initial = log(1.96), fixed = T))
  # fixed.theta.beta <- list(prec = list(initial = log(202), fixed = T))
  # fixed.theta.kappa <- list(prec = list(initial = log(336), fixed = T))
  # fixed.theta.epsilon <- list(prec = list(initial = log(420), fixed = T))
  
  loggamma.prior <- list(prec = list(prior = 'loggamma', param = c(1,0.00005), initial = log(1)))
  loggamma.prior.high.variance <- list(prec = list(prior = 'loggamma', param = c(1,0.005), initial = log(1)))
  
  comp = ~ -1 +
    Int(1, prec.linear = 0.001, mean.linear = 0) +
    alpha(x, model = "rw1", hyper = loggamma.prior, constr = TRUE) +
    beta(x.c, model = "iid", hyper = loggamma.prior, extraconstr = list(A = A.beta, e = e.beta)) +
    kappa(t, model = "rw1", hyper = loggamma.prior.high.variance, constr = TRUE) + 
    epsilon(xt, model = "iid", hyper = loggamma.prior, constr = FALSE)
  
  formula = Y ~ Int + alpha + beta*kappa + epsilon
  
  likelihood = like(formula = formula, family = "poisson", data = obs, E = obs$E)
  
  c.compute <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE, return.marginals.predictor = TRUE)  # control.compute
  
  res.inlabru = bru(components = comp,
                    likelihood, 
                    options = list(verbose = F,
                                   bru_verbose = 4, 
                                   num.threads = "1:1",
                                   control.compute = c.compute,
                                   bru_max_iter=max_iter,
                                   control.predictor = list(link = 1)
                    ))
  
  return(res.inlabru)
}

#   ----   Run inlabru ----

# here, we run inlabru
res.inlabru <- inlabru.pois.gp.rw1(obs, max_iter = 100)

#  save inlabru object
save(res.inlabru, file = file.path(output.path, "res_inlabru.RData"))

#   ----   Functions for plotting results   ----

# nice color palette
palette <- c('#70A4D4', '#ECC64B', '#93AD80', '#da9124', '#696B8D',
             '#3290c1', '#5d8060', '#D7B36A', '#826133', '#A85150')

# makes plot of inlabru results, as well as summarizing the inlabru results for further plotting
plot.inlabru.vs.underlying.synthetic.cancer<- function(
  res.inlabru, underlying.effects, path.to.storage="",
  cutoff_alpha = 1000, cutoff_beta = 1000, cutoff_kappa = 1000, cutoff_epsilon=1000,
  save=FALSE, pdf = TRUE, png = TRUE, cohort = FALSE){
  
  obs <- underlying.effects$obs
  nx <- underlying.effects$nx
  nt <- underlying.effects$nt
  
  data.alpha = cbind(res.inlabru$summary.random$alpha,
                     alpha.true = underlying.effects$alpha.true[res.inlabru$summary.random$alpha$ID + 1])
  p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = alpha.true, color = "True value", fill = "True value"), size = 1) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 1) + 
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(title="Alpha - inlabru", x = "x", y='')
  
  data.beta = cbind(res.inlabru$summary.random$beta,
                    beta.true = underlying.effects$beta.true[res.inlabru$summary.random$beta$ID + 1])
  print(data.beta)
  p.beta <- ggplot(data = data.beta, aes(x = ID)) + 
    geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated", color = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = beta.true, color = "True value", fill = "True value"), size = 1) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 1) + 
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(x = "x", y = "beta", title = "Beta - inlabru")
  
  data.kappa = cbind(res.inlabru$summary.random$kappa,
                     kappa.true = underlying.effects$kappa.true[res.inlabru$summary.random$kappa$ID + 1])
  p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = kappa.true, color = "True value", fill = "True value"), size = 1) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 1) + 
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(x = "t", y = "kappa", title = "Kappa - inlabru")
  
  data.fixed  <- data.frame(res.inlabru$marginals.fixed)
  
  p.intercept <- ggplot(data.fixed) + 
    geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Estimated", fill="Estimated")) + 
    geom_vline(aes(xintercept = underlying.effects$intercept, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    theme_classic() + 
    labs(x = "Value of phi", y = " ", title = "Intercept - inlabru")
  
  p.random.effects <- (p.intercept | p.alpha ) / (p.beta | p.kappa ) + 
    plot_layout(guides = "collect")
  
  if(save){
    save.figure(p.random.effects, name = "random_effects_inlabru", path = path.to.storage,
                pdf=pdf, png=png)
  }
  
  data.eta <- data.frame(eta.sim = res.inlabru$summary.linear.predictor$mean[1:length(obs$eta)],
                         `0.025quant` = res.inlabru$summary.linear.predictor$`0.025quant`[1:length(obs$eta)],
                         `0.975quant` = res.inlabru$summary.linear.predictor$`0.975quant`[1:length(obs$eta)]) %>%
    mutate(true.eta = obs$eta) %>%
    mutate(xt = obs$xt, x = obs$x, t = obs$t)
  
  p.eta <- ggplot(data = data.eta) +
    geom_point(aes(x = eta.sim, y = true.eta), color = palette[1]) + 
    labs(x="Estimated eta", y="True value for eta", title = "Eta")
  
  p.eta.2 <- ggplot(data = data.eta) +
    geom_line(aes(x=xt, y = eta.sim, color="Estimated")) +
    geom_line(aes(x=xt, y = true.eta, color="True")) +
    labs(x=" ", y="Eta", title="Eta- inlabru")
  
  p.eta.t <- ggplot(data = data.eta) + 
    geom_line(aes(x = x, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = x, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.eta.x <- ggplot(data = data.eta) + 
    geom_line(aes(x = t, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = t, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each age") + 
    facet_wrap(~x)
  
  p.eta.xt <- (p.eta | p.eta.2) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.xt, name="eta_xt_inlabru", path=path.to.storage, pdf=pdf, png=png)
  }
  
  p.eta.facet <- (p.eta.t | p.eta.x) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.facet, name="eta_facet_inlabru", path=path.to.storage, pdf=pdf, png = png)
  }
  
  data.mr <- data.frame(mr.sim = res.inlabru$summary.fitted.values$mean[1:length(obs$eta)]) %>%
    mutate(true.mr = obs$mr) %>%
    mutate(xt = obs$xt, x = obs$x, t = obs$t)
  
  p.mr <- ggplot(data = data.mr) +
    geom_point(aes(x = mr.sim, y = true.mr), color = palette[1]) + 
    labs(x="Estimated mortality rate", y="Observed mortality rate", title = "Mortality rate")
  
  p.mr.2 <- ggplot(data = data.mr) +
    geom_line(aes(x=xt, y = mr.sim, color="Estimated")) +
    geom_line(aes(x=xt, y = true.mr, color="True")) +
    labs(x=" ", y="Mortality rate", title="Mortality rate - inlabru")
  
  p.mr.t <- ggplot(data = data.mr) + 
    geom_line(aes(x = x, y = mr.sim, color = "Estimated")) +
    geom_line(aes(x = x, y = true.mr, color = "True")) +
    labs(x = " ", y = " ", title = "Mortality rate - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.mr.x <- ggplot(data = data.mr) + 
    geom_line(aes(x = t, y = mr.sim, color = "Estimated")) +
    geom_line(aes(x = t, y = true.mr, color = "True")) +
    labs(x = " ", y = " ", title = "Mortality rate - inlabru, for each age") + 
    facet_wrap(~x)
  
  p.mr.xt <- (p.mr | p.mr.2) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.mr.xt, name="mr_xt_inlabru", path=path.to.storage, pdf=pdf, png=png)
  }
  
  p.mr.facet <- (p.mr.t | p.mr.x) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.mr.facet, name="mr_facet_inlabru", path=path.to.storage, pdf = pdf, png=png)
  }
  
  data.hyperpar <- data.frame(res.inlabru$marginals.hyperpar)
  
  p.alpha.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.alpha.x < cutoff_alpha)) + 
    geom_area(aes(x = Precision.for.alpha.x, y = Precision.for.alpha.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_vline(aes(xintercept = res.inlabru$summary.hyperpar$mean[1], color = "Inlabru", fill = "Inlabru")) +
    geom_vline(aes(xintercept = underlying.effects$tau.alpha.true, color = "Observed", fill = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = " ", y = " ", title = "Precision of alpha")
  
  if(save){
    save.figure(p.alpha.prec, name="alpha_prec_inlabru", path=path.to.storage, pdf = pdf, png=png)
  }
  
  p.beta.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.beta.x < cutoff_beta)) + 
    geom_area(aes(x = Precision.for.beta.x, y = Precision.for.beta.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_vline(aes(xintercept = res.inlabru$summary.hyperpar$mean[2], color = "Inlabru", fill = "Inlabru")) + 
    geom_vline(aes(xintercept = underlying.effects$tau.beta.true, color = "Observed", fill = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = " ", y = " ", title = "Precision of beta")
  
  if(save){
    save.figure(p.beta.prec, name="beta_prec_inlabru", path=path.to.storage, pdf = pdf, png=png)
  }
  
  p.kappa.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.kappa.x < cutoff_kappa)) + 
    geom_area(aes(x = Precision.for.kappa.x, y = Precision.for.kappa.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_vline(aes(xintercept = res.inlabru$summary.hyperpar$mean[3], color = "Inlabru", fill = "Inlabru")) + 
    geom_vline(aes(xintercept = underlying.effects$tau.kappa.true, color = "Observed", fill = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = " ", y = " ", title = "Precision of kappa")
  
  if(cohort){
    p.gamma.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.gamma.x < cutoff_gamma)) + 
      geom_area(aes(x = Precision.for.gamma.x, y = Precision.for.gamma.y), color = palette[1], fill = palette[1], alpha = 0.5) + 
      geom_vline(aes(xintercept = res.inlabru$summary.hyperpar$mean[4]), color = palette[1]) + 
      labs(x = " ", y = " ", title = "Precision of gamma")
  }
  
  p.epsilon.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.epsilon.x < cutoff_epsilon)) + 
    geom_area(aes(x = Precision.for.epsilon.x, y = Precision.for.epsilon.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$tau.epsilon.true, color = "Observed", fill = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = " ", y = " ", title = "Precision of epsilon")
  
  if(save){
    save.figure(p.epsilon.prec, name="beta_epsilon_inlabru", path=path.to.storage, pdf = pdf, png=png)
  }
  
  if(cohort){
    p.hyperpars <- (p.alpha.prec | p.beta.prec | p.kappa.prec)/(p.gamma.prec | p.epsilon.prec) + plot_layout(guides = "collect")
  } else {
    p.hyperpars <- (p.alpha.prec | p.beta.prec )/(p.kappa.prec | p.epsilon.prec) + plot_layout(guides = "collect")
  }
  
  if(save){
    save.figure(p.hyperpars, name="hypers_inlabru", path=path.to.storage, pdf = pdf, png=png)
  }
  
  
  plots <- list(p.alpha = p.alpha, p.beta = p.beta, p.kappa = p.kappa,
                p.eta = p.eta,
                p.eta.2 = p.eta.2, p.eta.x = p.eta.x, 
                p.eta.t = p.eta.t, p.intercept = p.intercept,
                p.random.effects = p.random.effects,
                p.eta.all = p.eta.xt,
                p.eta.facet = p.eta.facet)
  
  summaries <- list(data.alpha = data.alpha,
                    data.beta = data.beta,
                    data.kappa = data.kappa,
                    data.eta = data.eta,
                    intercept = res.inlabru$summary.fixed$mean[1]
  )
  return(list(plots=plots, summaries=summaries))
}

# the following three functions makes summaries from the stan results and plots them
produce.stan.plots <- function(stan_df, underlying.effects, plot.func, save.func, 
                               path.to.storage, summaries.func = produce.summaries.stan){
  obs <- underlying.effects$obs
  
  summaries <- summaries.func(stan_df, obs, underlying.effects)
  
  plots <- plot.func(stan_df, obs, underlying.effects, summaries)
  
  save.func(plots, path.to.storage)
  
  res = list(plots = plots, summaries = summaries)
  return(res)
}

plot.stan.vs.underlying.synthetic.cancer <- function(stan_df, obs, underlying.effects, summaries){
  
  plot_intercept <- ggplot(data=stan_lc_df) +
    geom_vline(aes(xintercept = mean[6], color = "estimated")) + 
    geom_vline(aes(xintercept = `2.5%`[6], color = "estimated"), alpha = 0.5) +
    geom_vline(aes(xintercept = `97.5%`[6], color = "estimated"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$intercept, color = "true")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) + 
    labs(title="Stan - intercept", x = "Value of intercept", y = " ")
  
  summary_alpha <- summaries$summary_alpha
  
  plot_alpha <- ggplot(data=summary_alpha) +
    geom_point(aes(x=index, y=mean, color="estimated")) + 
    geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
    geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha=0.5) +
    geom_point(aes(x=index, y=true_alpha, color="true")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) + 
    labs(title="Stan - alpha", x = "x", y = " ")
  
  summary_beta <- summaries$summary_beta
  
  plot_beta <- ggplot(data=summary_beta) +
    geom_point(aes(x=index, y=mean, color="estimated")) + 
    geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
    geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha = 0.5) +
    geom_point(aes(x=index, y=true_beta, color="true")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) + 
    labs(title="Stan - beta", x = "x", y = " ")
  
  summary_kappa <- summaries$summary_kappa
  print(summary_kappa)
  
  plot_kappa <- ggplot(data=summary_kappa) +
    geom_point(aes(x=index, y=mean, color="estimated")) + 
    geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
    geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha=0.5) +
    geom_point(aes(x=index, y=true_kappa, color="true")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) + 
    labs(title="Stan - kappa", x = "t", y = " ")
  
  summary_eta <- summaries$summary_eta
  
  plot_eta <- ggplot(data=summary_eta) +
    geom_line(aes(x=index, y=mean, color="estimated")) + 
    geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
    geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha=0.5) +
    geom_line(aes(x=index, y=true_eta, color="true"), alpha = 0.5) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) + 
    labs(title="Stan - intercept", x = "x,t", y = " ")
  
  plots <- list(p.alpha = plot_alpha,
                p.beta = plot_beta, 
                p.kappa = plot_kappa,
                p.intercept = plot_intercept,
                p.eta = plot_eta)
  return(plots)
  
}

save.stan.plots.lc.rw2 <- function(plots, path.to.storage=""){
  
  p.alpha <- plots$p.alpha
  p.beta <- plots$p.beta
  p.eta <- plots$p.eta
  p.kappa <- plots$p.kappa
  p.intercept <- plots$p.intercept
  
  p.random.effects <- (p.intercept | p.alpha) / (p.beta | p.kappa) +
    plot_layout(guides="collect")
  
  save.figure(p.random.effects, name = "random_effects_stan", path = path.to.storage)
  save.figure(p.eta, name = "eta_stan", path = path.to.storage)
}

produce.summaries.stan.traditional <- function(stan_df, obs, underlying.effects){
  
  summary_fixed <- stan_df[5:6,]
  
  summary_alpha <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('alpha', parameter)) %>%
    filter(!grepl('tau_alpha', parameter)) %>%
    filter(!grepl('theta_alpha', parameter)) %>%
    filter(!grepl('alpha_raw', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_alpha = obs$alpha[index])
  
  summary_beta <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('beta', parameter)) %>%
    filter(!grepl('tau_beta', parameter)) %>%
    filter(!grepl('theta_beta', parameter)) %>%
    filter(!grepl('beta_raw', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_beta = underlying.effects$beta.true[index])
  
  summary_beta_raw <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('beta_raw', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_beta = underlying.effects$beta.true[index])
  
  summary_kappa <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('kappa', parameter)) %>%
    filter(!grepl('tau_kappa', parameter)) %>%
    filter(!grepl('theta_kappa', parameter)) %>%
    filter(!grepl('kappa_raw', parameter)) %>%
    filter(!grepl('kappa_0', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_kappa = underlying.effects$kappa.true[index]) %>%
    mutate(kappa_drifted = underlying.effects$kappa.drifted[index])
  
  summary_eta <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('eta', parameter)) %>%
    filter(!grepl('beta', parameter)) %>%
    filter(!grepl('theta', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_eta = obs$eta.no.error) %>%
    mutate(true_eta_error  = obs$eta) %>%
    mutate(xt = obs$xt, x = obs$x, t = obs$t)
  
  summaries <- list(
    summary_alpha=summary_alpha,
    summary_beta=summary_beta,
    summary_beta_raw=summary_beta,
    summary_kappa=summary_kappa,
    summary_eta=summary_eta,
    summary_fixed=summary_fixed
  )
  
  return(summaries)
}

# produces plots for comparison of stan and inlabru results, and saves them 
produce.compared.plots <- function(
  stan.summaries, stan.marginals, inlabru.summaries, res.inlabru, underlying.effects, plot.func, save.func,
  path.to.storage){
  
  plots <- plot.func(stan.summaries, stan.marginals, inlabru.summaries, res.inlabru, underlying.effects)
  save.func(plots, path.to.storage)
  return(plots)
}

# this function takes output that the previous functions have produced, and produces plots that compare inlabru and stan - mainly by mean and confidence bounds
plot.inlabru.stan.compared <- function(stan.summaries,
                                           stan.marginals,
                                           inlabru.summaries,
                                           res.inlabru,
                                           underlying.effects,
                                           tau.alpha.cutoff = 20,
                                           tau.beta.cutoff = 50000,
                                           tau.kappa.cutoff = 50000,
                                           tau.epsilon.cutoff = 1500
){
  #' Produces plots with comparison of estimation results from inlabru and STAN
  #' 
  #'@param stan.summaries (list<data.frame>) summaries of STAN results
  #'@param stan.marginals (list<array>) Hamiltonian MC samples from STAN
  #'@param inlabru.summaries (list<data.frame>) summaries of inlabru results
  #'@param res.inlabru (bru object) raw inlabru results
  #'@param underlying.effects (list) underlying data for which analysis is run
  #'@param cohort (boolean) whether analysis includes cohort effect 
  #'
  obs <- underlying.effects$obs
  
  intercept.marginal <- data.frame(int = stan.marginals$intercept_draws)
  
  inlabru.data.fixed = data.frame(res.inlabru$marginals.fixed)
  
  #  ----   intercept   ----
  if(length(stan.marginals$intercept_draws) > 0){
    p.intercept <- ggplot() + 
      geom_histogram(data = intercept.marginal, aes(x = int, y = after_stat(density), color = "Stan", fill = "Stan"), bins=200, alpha = 0.5) + 
      geom_area(data=inlabru.data.fixed, aes(x = Int.x, y = Int.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.4, size = 0.5) + 
      geom_vline(aes(xintercept = underlying.effects$age.intercept.true, color="True", fill="True")) +
      scale_color_manual(name = " ", values = palette) + 
      scale_fill_manual(name = " ", values = palette) +
      theme_classic() + 
      labs(x = "Value of intercept", y = " ", title = "Intercept")
  } else {
    p.intercept <- ggplot(data = data.frame(a = 1, b = 2)) + geom_point(aes(x = a, y = b)) + labs(title = "no intercept data for stan")
  }
  
  
  # ---   alpha   ----
  
  p.alpha <- ggplot() + 
    geom_ribbon(data=inlabru.summaries$data.alpha, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
    geom_point(data=inlabru.summaries$data.alpha, aes(x = ID, y = mean, color = "Inlabru", fill = "Inlabru"), size = 0.5) + 
    geom_point(data=stan.summaries$summary_alpha, aes(x=index - 1, y = mean, fill="Stan", color="Stan"), size=0.5) + 
    geom_line(data=stan.summaries$summary_alpha, aes(x=index - 1, y = `2.5%`, fill = "Stan", color="Stan"), alpha=0.5) + 
    geom_line(data=stan.summaries$summary_alpha, aes(x=index - 1, y = `97.5%`, fill="Stan", color="Stan"), alpha=0.5) +
    
    geom_point(data=inlabru.summaries$data.alpha, 
               aes(x = ID, y = underlying.effects$alpha.true, color = "True", fill = "True"), size = 0.5) + 
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(title="Alpha", x = "x", y='')
  
  
  # ---   beta   ----
  
  p.beta <- ggplot() + 
    geom_errorbar(data = inlabru.summaries$data.beta, aes(x = ID + 0.1, ymin = `0.025quant`, ymax = `0.975quant`, color = "Inlabru", fill = "Inlabru")) +
    #geom_point(data=inlabru.summaries$data.beta, aes(x = ID, y = mean, color = "Inlabru", fill = "Inlabru"), size = 0.5) +
    #geom_point(data=stan.summaries$summary_beta, aes(x=index - 1, y=mean, fill="Stan", color="Stan"), size=0.5) + 
    geom_errorbar(data = stan.summaries$summary_beta, aes(x = index - 1 - 0.1, ymin = `2.5%`, ymax = `97.5%`, color = "Stan", fill = "Stan")) +
    geom_point(data = inlabru.summaries$data.beta, aes(x = ID, y = underlying.effects$beta.true, color = "True", fill = "True"), shape = 4) +
    
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    scale_shape_manual(name = "") + 
    theme_classic() + 
    labs(title="Beta", x = "x", y='')
  
  #   ----   kappa   ----  
  p.kappa <- ggplot() + 
    geom_ribbon(data=inlabru.summaries$data.kappa, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
    geom_point(data=inlabru.summaries$data.kappa, aes(x = ID, y = mean, color = "Inlabru", fill = "Inlabru"), size = 0.5) + 
    
    geom_point(data=stan.summaries$summary_kappa, aes(x=index - 1, y=mean, fill="Stan",color="Stan"), size=0.5) + 
    geom_line(data=stan.summaries$summary_kappa, aes(x=index - 1, y=`2.5%`, fill="Stan", color="Stan"), alpha=0.5) + 
    geom_line(data=stan.summaries$summary_kappa, aes(x=index - 1, y=`97.5%`, fill="Stan", color="Stan"), alpha=0.5) +
    
    geom_point(data=inlabru.summaries$data.kappa,
               aes(x = ID, y = underlying.effects$kappa.true, color = "True", fill = "True"), size = 0.5) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(title="Kappa", x = "t", y='')
  
  #   ----   eta   ----
  eta.stan.inlabru<- inlabru.summaries$data.eta %>% left_join(stan.summaries$summary_eta, by = c("xt" = "xt")) %>%
    mutate(diff.eta = eta.sim - mean)
  
  p.eta.stan.inlabru <- ggplot(data = eta.stan.inlabru) + 
    geom_point(aes(x = xt, y = diff.eta)) + 
    labs(title = "Difference between Inlabru and Stan estimation of eta", x = "xt", y = "Difference")
  
  p.eta.stan.v.inlabru <- ggplot(data = eta.stan.inlabru) + 
    geom_point(aes(x = eta.sim, y = mean), color = palette[1]) +
    theme_classic() + 
    labs(x= "Eta, by inlabru", y = "Eta, by Stan", title = "Estimated eta by Inlabru (x) and Stan (y)")
  
  
  p.eta <- ggplot() +
    geom_point(data=inlabru.summaries$data.eta, aes(x = eta.sim, y = true.eta, color = "Inlabru"), alpha = 0.5) + 
    geom_point(data=stan.summaries$summary_eta, aes(x = mean, y = true_eta, color = "Stan"), alpha = 0.5) + 
    scale_color_manual(name = " ", values = palette) + 
    theme_classic() + 
    labs(x="Estimated eta", y="True value for eta", title = "Eta")
  
  p.eta.2 <- ggplot() +
    geom_point(data = inlabru.summaries$data.eta, aes(x=xt, y = eta.sim, color="Inlabru", fill="Inlabru"), size=0.5, alpha = 0.5) +
    #geom_ribbon(data = inlabru.summaries$data.eta, aes(x = xt, ymin = `0.025quant`, ymax=`0.975quant`, fill = "Inlabru"), alpha = 0.5)  +
    geom_ribbon(data = inlabru.summaries$data.eta, aes(x = xt, ymin = X0.025quant, ymax=X0.975quant, fill = "Inlabru"), alpha = 0.5)  +
    geom_point(data = stan.summaries$summary_eta, aes(x=xt, y = mean, color="Stan", fill="Stan"), size=0.5, alpha = 0.5) +
    geom_ribbon(data = stan.summaries$summary_eta, aes(x = xt, ymin = `2.5%`, ymax=`97.5%`, fill="Stan"), alpha=0.5)  +
    geom_line(data = inlabru.summaries$data.eta, aes(x=xt, y = true.eta, color="True", fill="True")) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(x=" ", y="Eta", title="Eta")
  
  p.eta.t <- ggplot() + 
    geom_point(data = inlabru.summaries$data.eta, aes(x=x, y = eta.sim, color = "Inlabru", fill="Inlabru"), size=0.5) +
    #geom_ribbon(data = inlabru.summaries$data.eta, aes(x = x, ymin = `0.025quant`, ymax=`0.975quant`, fill = "Inlabru"), alpha = 0.5)  +
    geom_ribbon(data = inlabru.summaries$data.eta, aes(x = x, ymin = X0.025quant, ymax=X0.975quant, fill = "Inlabru"), alpha = 0.5)  +
    geom_point(data = stan.summaries$summary_eta, aes(x = x, y = mean, color = "Stan", fill="Stan"), size=0.5) +
    geom_ribbon(data = stan.summaries$summary_eta, aes(x = x, ymin = `2.5%`, ymax=`97.5%`, fill="Stan"), alpha=0.5)  +
    geom_line(data = inlabru.summaries$data.eta, aes(x = x, y = true.eta, color = "True", fill="True")) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(x = " ", y = " ", title = "Eta - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.eta.x <- ggplot() + 
    geom_point(data = inlabru.summaries$data.eta, aes(x = t, y = eta.sim, color = "Inlabru", fill="Inlabru"), size=0.5) +
    #geom_ribbon(data = inlabru.summaries$data.eta, aes(x = t, ymin = `0.025quant`, ymax=`0.975quant`, fill = "Inlabru"), alpha = 0.5)  +
    geom_ribbon(data = inlabru.summaries$data.eta, aes(x = t, ymin = X0.025quant, ymax=X0.975quant, fill = "Inlabru"), alpha = 0.5)  +
    geom_point(data=stan.summaries$summary_eta, aes(x = t, y = mean, color = "Stan", fill="Stan"), size=0.5) +
    geom_ribbon(data = stan.summaries$summary_eta, aes(x = t, ymin = `2.5%`, ymax=`97.5%`, fill="Stan"), alpha=0.5)  +
    geom_line(data = inlabru.summaries$data.eta, aes(x = t, y = true.eta, color = "True", fill="True")) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(x = " ", y = " ", title = "Eta - inlabru, for each age") + 
    facet_wrap(~x)
  
  #   ----   hyperparameters: precisions   ----   
  #  tau alpha 
  
  tau.alpha.stan <- data.frame(tau = stan.marginals$tau_alpha_draws) %>%
    filter(tau < tau.alpha.cutoff)
  tau.alpha.inlabru <- data.frame(res.inlabru$marginals.hyperpar$`Precision for alpha`) %>%
    filter(x < tau.alpha.cutoff)
  
  p.tau.alpha <- ggplot() +
    geom_area(data = tau.alpha.inlabru, aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) +
    geom_histogram(data = tau.alpha.stan, aes(x = tau, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins=100) +
    geom_vline(data = tau.alpha.inlabru, aes(xintercept = underlying.effects$tau.alpha.true, color = "Observed", fill  = "Observed")) +
    scale_color_manual(name = " ", values = palette) +
    scale_fill_manual(name = " ", values = palette) +
    theme_classic() +
    labs(x = "Value of precision of alpha", y = " ", title = "Precision of Alpha")
  
  #  tau beta
  
  tau.beta.stan <- data.frame(tau = stan.marginals$tau_beta_draws) %>%
    filter(tau < tau.beta.cutoff)
  tau.beta.inlabru <- data.frame(res.inlabru$marginals.hyperpar$`Precision for beta`) %>%
    filter(x < tau.beta.cutoff)
  
  p.tau.beta <- ggplot() +
    geom_area(data = tau.beta.inlabru, aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) +
    geom_histogram(data = tau.beta.stan, aes(x = tau, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins=100) +
    geom_vline(data = tau.beta.inlabru, aes(xintercept = underlying.effects$tau.beta.true, color = "Observed", fill = "Observed")) +
    scale_color_manual(name = " ", values = palette) +
    scale_fill_manual(name = " ", values = palette) +
    theme_classic() +
    labs(x = "Value of precision of beta", y = " ", title = "Precision of Beta")
  
  # tau kappa
  tau.kappa.stan <- data.frame(tau = stan.marginals$tau_kappa_draws) %>%
    filter(tau < tau.kappa.cutoff)
  tau.kappa.inlabru <- data.frame(res.inlabru$marginals.hyperpar$`Precision for kappa`) %>%
    filter(x < tau.kappa.cutoff)
  
  p.tau.kappa <- ggplot() +
    geom_area(data = tau.kappa.inlabru, aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) +
    geom_histogram(data = tau.kappa.stan, aes(x = tau, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins = 100) +
    geom_vline(data = tau.kappa.inlabru, aes(xintercept = underlying.effects$tau.kappa.true, color = "Observed", fill = "Observed")) +
    scale_color_manual(name = " ", values = palette) +
    scale_fill_manual(name = " ", values = palette) +
    theme_classic() +
    labs(x = "Value of precision of kappa", y = " ", title = "Precision of Kappa")
  
  # # tau epsilon
  if(length(stan.marginals$tau_epsilon_draws) > 0){
    tau.epsilon.stan <- data.frame(tau = stan.marginals$tau_epsilon_draws) %>%
      filter(tau < tau.epsilon.cutoff)
    tau.epsilon.inlabru <- data.frame(res.inlabru$marginals.hyperpar$`Precision for epsilon`) %>%
      filter(x < tau.epsilon.cutoff)
    
    p.tau.epsilon <- ggplot() +
      geom_area(data = tau.epsilon.inlabru, aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) +
      geom_histogram(data = tau.epsilon.stan, aes(x = tau, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins = 100) +
      geom_vline(data = tau.epsilon.inlabru, aes(xintercept = underlying.effects$tau.epsilon.true, color = "Observed", fill = "Observed")) +
      scale_color_manual(name = " ", values = palette) +
      scale_fill_manual(name = " ", values = palette) +
      theme_classic() +
      labs(x = "Value of precision of epsilon", y = " ", title = "Precision of Epsilon")
  } else {
    print("No tau epsilon!")
    dummy.data  <- data.frame(a = c(1.0, 2.0), b = c(1.0, 2.0))
    print(dummy.data)
    p.tau.epsilon <- ggplot(data = dummy.data) + geom_point(aes(x = a, y = b)) + labs(title = "No available tau epsilon")
  }
  #   ----   Returns   ----
  
  plots <- list(p.intercept = p.intercept, 
                p.alpha = p.alpha, 
                p.beta = p.beta,
                p.kappa = p.kappa,
                p.eta = p.eta,
                p.eta.2 = p.eta.2,
                p.eta.t = p.eta.t,
                p.eta.x = p.eta.x,
                p.tau.alpha = p.tau.alpha,
                p.tau.beta = p.tau.beta,
                p.tau.kappa = p.tau.kappa,
                p.tau.epsilon = p.tau.epsilon,
                p.eta.stan.inlabru = p.eta.stan.inlabru,
                p.eta.stan.v.inlabru = p.eta.stan.v.inlabru)
  
  return(plots)
}

# saves plots for comparioson of stan and inlabru results
save.compared <- function(plots, path.to.storage, cohort=TRUE, pdf=T, png=F){
  p.random.effects <- (plots$p.intercept | plots$p.alpha )/(plots$p.beta | plots$p.kappa) + 
    plot_layout(guides="collect")
  
  save.figure(p.random.effects, name = "random_effects_comparison", path = path.to.storage, pdf=pdf, png=png)
  
  p.hypers <- (plots$p.tau.alpha | plots$p.tau.beta) / (plots$p.tau.kappa | plots$p.tau.gamma | plots$p.tau.epsilon) + plot_layout(guides = "collect")
  
  save.figure(p.hypers, name = "hypers_comparison", path = path.to.storage, pdf=pdf, png=png)
  
  p.eta.xt <- (plots$p.eta | plots$p.eta.2) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  save.figure(p.eta.xt, name = "eta_xt_comparison", path = path.to.storage)
  
  p.eta.facet <- (plots$p.eta.x | plots$p.eta.t) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  save.figure(p.eta.facet, name = "eta_facet_comparison", path = path.to.storage)
  
  save.figure(plots$p.eta.x, "eta_x_comparison", path = path.to.storage, pdf=pdf, png=png)
  save.figure(plots$p.eta.t, "eta_t_comparison", path = path.to.storage, pdf=pdf, png=png)
  
  if("p.eta.stan.inlabru" %in% names(plots)){
    save.figure(plots$p.eta.stan.inlabru, "eta_difference", path = path.to.storage, pdf=pdf, png=png)
  }
  
  if("p.eta.stan.v.inlabru" %in% names(plots)){
    save.figure(plots$p.eta.stan.v.inlabru, "eta_stan_v_inlabru", path = path.to.storage, pdf=pdf, png=png)
  }
}

# functions plotting compared marginals of some values of the predictor, beta and kappa
plot.predictor.inlabru.stan.compared <- function(
  res.inlabru, stan.predictor.df,
  path.to.storage, columns, pdf=T, png=F) {
  
  pred.1.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.001)
  
  p.predictor.1 <- ggplot() + 
    geom_area(data = pred.1.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.predictor.df, aes(x = X1, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Predictor[1], Inlabru and stan", x = " ", y = " ")
  
  pred.2.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.064)
  
  p.predictor.2 <- ggplot() + 
    geom_area(data = pred.2.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.predictor.df, aes(x = X64, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Predictor[64], Inlabru and stan", x = " ", y = " ")
  
  pred.3.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.128)
  
  p.predictor.3 <- ggplot() + 
    geom_area(data = pred.3.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.predictor.df, aes(x = X128, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Predictor[128], Inlabru and stan", x = " ", y = " ")
  
  pred.4.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.192)
  
  p.predictor.4 <- ggplot() + 
    geom_area(data = pred.4.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.predictor.df, aes(x = X192, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Predictor[192], Inlabru and stan", x = " ", y = " ")
  
  pred.5.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.256)
  
  p.predictor.5 <- ggplot() + 
    geom_area(data = pred.5.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.predictor.df, aes(x = X256, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Predictor[256], Inlabru and stan", x = " ", y = " ")
  
  pred.6.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.324)
  
  p.predictor.6 <- ggplot() + 
    geom_area(data = pred.6.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.predictor.df, aes(x = X324, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Predictor[324], Inlabru and stan", x = " ", y = " ")
  
  p.predictor <- (p.predictor.1 | p.predictor.2 | p.predictor.3)/(p.predictor.4 | p.predictor.5 | p.predictor.6) + 
    plot_layout(guides = "collect")
  
  save.figure(p.predictor, "predictor_marginals_comparison", path = path.to.storage, png = png, pdf = pdf)
}

plot.beta.inlabru.stan.compared <- function(
  res.inlabru, stan.predictor.df,
  path.to.storage, columns, pdf=T, png=F, a45=F) {
  
  pred.1.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.1)
  
  p.beta.1 <- ggplot() + 
    geom_area(data = pred.1.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.predictor.df, aes(x = X1, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Beta[1], Inlabru and stan", x = " ", y = " ")
  
  pred.3.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.3)
  
  p.beta.3 <- ggplot() + 
    geom_area(data = pred.3.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.predictor.df, aes(x = X3, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Beta[3], Inlabru and stan", x = " ", y = " ")
  
  pred.5.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.5)
  
  p.beta.5 <- ggplot() + 
    geom_area(data = pred.5.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.predictor.df, aes(x = X5, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Beta[5], Inlabru and stan", x = " ", y = " ")
  
  pred.7.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.7)
  
  p.beta.7 <- ggplot() + 
    geom_area(data = pred.7.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.predictor.df, aes(x = X7, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Beta[7], Inlabru and stan", x = " ", y = " ")
  
  pred.9.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.9)
  
  p.beta.9 <- ggplot() + 
    geom_area(data = pred.9.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.predictor.df, aes(x = X9, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Beta[9], Inlabru and stan", x = " ", y = " ")
  
  pred.11.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.11)
  
  p.beta.11 <- ggplot() + 
    geom_area(data = pred.11.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.predictor.df, aes(x = X11, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Beta[11], Inlabru and stan", x = " ", y = " ")
  
  pred.13.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.13)
  
  p.beta.13 <- ggplot() + 
    geom_area(data = pred.13.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.predictor.df, aes(x = X13, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Beta[13], Inlabru and stan", x = " ", y = " ")
  
  pred.15.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.15)
  
  p.beta.15 <- ggplot() + 
    geom_area(data = pred.15.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.predictor.df, aes(x = X15, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Beta[15], Inlabru and stan", x = " ", y = " ")
  
  pred.17.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.17)
  
  p.beta.17 <- ggplot() + 
    geom_area(data = pred.17.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.predictor.df, aes(x = X17, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Beta[17], Inlabru and stan", x = " ", y = " ")
  
  p.beta <- (p.beta.1 | p.beta.3 | p.beta.5) / (p.beta.7 | p.beta.9 | p.beta.11) / (p.beta.13 | p.beta.15 | p.beta.17) + 
    plot_layout(guides = "collect")
  
  save.figure(p.beta, "beta_marginals_comparison", path = path.to.storage, png = png, pdf = pdf)
}

plot.kappa.inlabru.stan.compared <- function(
  res.inlabru, stan.predictor.df,
  path.to.storage, columns, pdf=T, png=F) {
  
  pred.1.inlabru <- data.frame(res.inlabru$marginals.random$kappa$index.1)
  
  p.kappa.1 <- ggplot() + 
    geom_area(data = pred.1.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.predictor.df, aes(x = X1, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Kappa[1], Inlabru and stan", x = " ", y = " ")
  
  pred.3.inlabru <- data.frame(res.inlabru$marginals.random$kappa$index.3)
  
  p.kappa.3 <- ggplot() + 
    geom_area(data = pred.3.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.predictor.df, aes(x = X3, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Kappa[3], Inlabru and stan", x = " ", y = " ")
  
  pred.5.inlabru <- data.frame(res.inlabru$marginals.random$kappa$index.5)
  
  p.kappa.5 <- ggplot() + 
    geom_area(data = pred.5.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.predictor.df, aes(x = X5, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Kappa[5], Inlabru and stan", x = " ", y = " ")
  
  pred.7.inlabru <- data.frame(res.inlabru$marginals.random$kappa$index.7)
  
  p.kappa.7 <- ggplot() + 
    geom_area(data = pred.7.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.predictor.df, aes(x = X7, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Kappa[7], Inlabru and stan", x = " ", y = " ")
  
  
  pred.9.inlabru <- data.frame(res.inlabru$marginals.random$kappa$index.9)
  
  p.kappa.9 <- ggplot() + 
    geom_area(data = pred.9.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.predictor.df, aes(x = X9, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Kappa[9], Inlabru and stan", x = " ", y = " ")
  
  pred.11.inlabru <- data.frame(res.inlabru$marginals.random$kappa$index.11)
  
  p.kappa.11 <- ggplot() + 
    geom_area(data = pred.11.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.predictor.df, aes(x = X11, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Kappa[11], Inlabru and stan", x = " ", y = " ")
  
  pred.13.inlabru <- data.frame(res.inlabru$marginals.random$kappa$index.13)
  
  p.kappa.13 <- ggplot() + 
    geom_area(data = pred.13.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.predictor.df, aes(x = X13, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Kappa[13], Inlabru and stan", x = " ", y = " ")
  
  pred.15.inlabru <- data.frame(res.inlabru$marginals.random$kappa$index.15)
  
  p.kappa.15 <- ggplot() + 
    geom_area(data = pred.15.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.predictor.df, aes(x = X15, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Kappa[15], Inlabru and stan", x = " ", y = " ")
  
  pred.17.inlabru <- data.frame(res.inlabru$marginals.random$kappa$index.17)
  
  p.kappa.17 <- ggplot() + 
    geom_area(data = pred.17.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.predictor.df, aes(x = X17, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Kappa[17], Inlabru and stan", x = " ", y = " ")
  
  p.kappa <- (p.kappa.1 | p.kappa.3 | p.kappa.5) / (p.kappa.7 | p.kappa.9 | p.kappa.11) / (p.kappa.13 | p.kappa.15 | p.kappa.17) + 
    plot_layout(guides = "collect")
  
  save.figure(p.kappa, "kappa_marginals_comparison", path = path.to.storage, png = png, pdf = pdf)
}

plot.epsilon.inlabru.stan.compared <- function(
  res.inlabru, stan.epsilon.df,
  path.to.storage, pdf=T, png=F) {
  
  pred.1.inlabru <- data.frame(res.inlabru$marginals.random$epsilon$index.1)
  
  p.epsilon.1 <- ggplot() + 
    geom_area(data = pred.1.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.epsilon.df, aes(x = X1, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Epsilon[1], Inlabru and stan", x = " ", y = " ")
  
  pred.2.inlabru <- data.frame(res.inlabru$marginals.random$epsilon$index.64)
  
  p.epsilon.2 <- ggplot() + 
    geom_area(data = pred.2.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.epsilon.df, aes(x = X64, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Epsilon[64], Inlabru and stan", x = " ", y = " ")
  
  pred.3.inlabru <- data.frame(res.inlabru$marginals.random$epsilon$index.128)
  
  p.epsilon.3 <- ggplot() + 
    geom_area(data = pred.3.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.epsilon.df, aes(x = X128, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Epsilon[128], Inlabru and stan", x = " ", y = " ")
  
  pred.4.inlabru <- data.frame(res.inlabru$marginals.random$epsilon$index.192)
  
  p.epsilon.4 <- ggplot() + 
    geom_area(data = pred.4.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.epsilon.df, aes(x = X192, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Epsilon[192], Inlabru and stan", x = " ", y = " ")
  
  pred.5.inlabru <- data.frame(res.inlabru$marginals.random$epsilon$index.256)
  
  p.epsilon.5 <- ggplot() + 
    geom_area(data = pred.5.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.epsilon.df, aes(x = X256, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Epsilon[256], Inlabru and stan", x = " ", y = " ")
  
  pred.6.inlabru <- data.frame(res.inlabru$marginals.random$epsilon$index.324)
  
  p.epsilon.6 <- ggplot() + 
    geom_area(data = pred.6.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = stan.epsilon.df, aes(x = X324, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Epsilon[324], Inlabru and stan", x = " ", y = " ")
  
  p.epsilon <- (p.epsilon.1 | p.epsilon.2 | p.epsilon.3)/(p.epsilon.4 | p.epsilon.5 | p.epsilon.6) + 
    plot_layout(guides = "collect")
  
  save.figure(p.epsilon, "epsilon_marginals_comparison", path = path.to.storage, png = png, pdf = pdf)
}

#   ----   Load and plot results ----

plots.summaries.inlabru <- plot.inlabru.vs.underlying.synthetic.cancer(
  res.inlabru,
  underlying.effects,
  path.to.storage = output.path,
  save=F)

load(file.path(stan.output, paste("stan_", investigation.name, ".Rda", sep = "")))

load(file=file.path(stan.output, "draws_intercept.RData"))
load(file=file.path(stan.output, "draws_tau_epsilon.RData"))
load(file.path(stan.output, "draws_tau_alpha.RData"))
load(file.path(stan.output, "draws_tau_beta.RData"))
load(file.path(stan.output, "draws_tau_kappa.RData"))
load(file.path(stan.output, "draws_alpha.RData"))
load(file.path(stan.output, "draws_beta.RData"))
load(file.path(stan.output, "draws_kappa.RData"))
load(file.path(stan.output, "draws_eta_100.RData"))
load(file.path(stan.output, "draws_eta.RData"))
load(file.path(stan.output, "draws_eta_reduced.RData"))
load(file.path(stan.output, "draws_epsilon.RData"))

stan.marginals <- list(intercept_draws = intercept_draws,
                       tau_epsilon_draws = tau_epsilon_draws,
                       tau_alpha_draws = tau_alpha_draws,
                       tau_beta_draws = tau_beta_draws,
                       tau_kappa_draws = tau_kappa_draws,
                       alpha_draws = alpha_draws,
                       beta_draws = beta_draws,
                       kappa_draws = kappa_draws,
                       eta_draws = eta_draws)

stan.res <- produce.stan.plots(stan_df=stan_lc_df,
                               underlying.effects=underlying.effects,
                               plot.func=plot.stan.vs.underlying.synthetic.cancer,
                               save.func=save.stan.plots.lc.rw2,
                               path.to.storage=output.path,
                               summaries.func=produce.summaries.stan.traditional)

plots_compared <- produce.compared.plots(
  stan.summaries = stan.res$summaries,
  stan.marginals = stan.marginals,
  inlabru.summaries = plots.summaries.inlabru$summaries,
  res.inlabru = res.inlabru,
  underlying.effects = underlying.effects,
  plot.func = function(...) {plot.inlabru.stan.compared(...,tau.beta.cutoff = 700, tau.kappa.cutoff = 500, tau.alpha.cutoff = 10)},
  save.func = save.compared,
  path.to.storage=output.path)

#   ----   Sample predictor   ----

stan.predictor.df <- data.frame(eta_draws)

plot.predictor.inlabru.stan.compared(res.inlabru, stan.predictor.df, path.to.storage = output.path)

#   ----   Plot marginals of random effects   ----

stan.beta.df <- data.frame(beta_draws)
plot.beta.inlabru.stan.compared(res.inlabru, stan.beta.df, path.to.storage = output.path)

stan.kappa.df <- data.frame(kappa_draws)
plot.kappa.inlabru.stan.compared(res.inlabru, stan.kappa.df, path.to.storage = output.path)

stan.epsilon.df <- data.frame(epsilon_draws)
plot.epsilon.inlabru.stan.compared(res.inlabru, stan.epsilon.df, path.to.storage = output.path)



