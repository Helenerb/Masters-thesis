# file containing functions for plotting stan results compared to true
# underlying effects

library(ggplot2)
library(patchwork)

#   ----   Source relevant functions

# source("../Functions/plotters.R")
# source("../Misc/palette.R")

# assume wokring directory at ---Master Thesis Code
source("Scripts/Functions/plotters.R")
source("Scripts/Misc/palette.R")

produce.summaries.stan <- function(stan_df, obs, underlying.effects){
  
  summary_fixed <- stan_df[6:7,]
  
  summary_alpha <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('alpha', parameter)) %>%
    filter(!grepl('tau_alpha', parameter)) %>%
    filter(!grepl('alpha_raw', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_alpha = obs$alpha[index])
  
  summary_beta <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('beta', parameter)) %>%
    filter(!grepl('tau_beta', parameter)) %>%
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
    filter(!grepl('kappa_raw', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_kappa = underlying.effects$kappa.true[index]) %>%
    mutate(kappa_drifted = true_kappa + underlying.effects$phi.true*(index-1))
  
  summary_gamma <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('gamma', parameter)) %>%
    filter(!grepl('tau_gamma', parameter)) %>%
    filter(!grepl('gamma_raw', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_gamma = underlying.effects$gamma.true[index]) 
  
  summary_eta <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('eta', parameter)) %>%
    filter(!grepl('beta', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_eta = obs$eta) %>%
    mutate(xt = obs$xt) %>%
    mutate(x = obs$x, t = obs$t)
  
  summaries <- list(
    summary_alpha=summary_alpha,
    summary_beta=summary_beta,
    summary_beta_raw=summary_beta,
    summary_kappa=summary_kappa,
    summary_gamma=summary_gamma,
    summary_eta=summary_eta,
    summary_fixed=summary_fixed
  )
  
  return(summaries)
}

produce.summaries.stan.lc <- function(stan_df, obs, underlying.effects){
  
  summary_fixed <- stan_df[5:6,]
  
  summary_alpha <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('alpha', parameter)) %>%
    filter(!grepl('tau_alpha', parameter)) %>%
    filter(!grepl('alpha_raw', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_alpha = obs$alpha[index])
  
  summary_beta <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('beta', parameter)) %>%
    filter(!grepl('tau_beta', parameter)) %>%
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
    filter(!grepl('kappa_raw', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_kappa = underlying.effects$kappa.true[index]) %>%
    mutate(kappa_drifted = true_kappa + underlying.effects$phi.true*(index-1))
  
  summary_eta <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('eta', parameter)) %>%
    filter(!grepl('beta', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_eta = obs$eta) %>%
    mutate(xt = obs$xt) %>%
    mutete(x = obs$x, t = obs$t)
  
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

produce.summaries.stan.hard <- function(stan_df, obs, underlying.effects){
  
  summary_fixed <- stan_df[6:7,]
  
  summary_alpha <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('alpha', parameter)) %>%
    filter(!grepl('tau_alpha', parameter)) %>%
    filter(!grepl('alpha_raw', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_alpha = obs$alpha[index])
  
  summary_beta <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('beta', parameter)) %>%
    filter(!grepl('tau_beta', parameter)) %>%
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
    filter(!grepl('kappa_raw', parameter)) %>%
    filter(!grepl('kappa_0', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_kappa = underlying.effects$kappa.true[index]) %>%
    mutate(kappa_drifted = true_kappa + underlying.effects$phi.true*(index-1))
  
  summary_gamma <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('gamma', parameter)) %>%
    filter(!grepl('tau_gamma', parameter)) %>%
    filter(!grepl('gamma_raw', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_gamma = underlying.effects$gamma.true[index]) 
  
  summary_eta <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('eta', parameter)) %>%
    filter(!grepl('beta', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_eta = obs$eta) %>%
    mutate(xt = obs$xt) %>%
    mutate(x = obs$x, t = obs$t)
  
  summaries <- list(
    summary_alpha=summary_alpha,
    summary_beta=summary_beta,
    summary_beta_raw=summary_beta,
    summary_kappa=summary_kappa,
    summary_gamma=summary_gamma,
    summary_eta=summary_eta,
    summary_fixed=summary_fixed
  )
  
  return(summaries)
}

produce.summaries.stan.hard.lc <- function(stan_df, obs, underlying.effects){
  
  summary_fixed <- stan_df[5:6,]
  
  summary_alpha <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('alpha', parameter)) %>%
    filter(!grepl('tau_alpha', parameter)) %>%
    filter(!grepl('alpha_raw', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_alpha = obs$alpha[index])
  
  summary_beta <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('beta', parameter)) %>%
    filter(!grepl('tau_beta', parameter)) %>%
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
    filter(!grepl('kappa_raw', parameter)) %>%
    filter(!grepl('kappa_0', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_kappa = underlying.effects$kappa.true[index]) %>%
    mutate(kappa_drifted = true_kappa + underlying.effects$phi.true*(index-1))
  
  summary_eta <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('eta', parameter)) %>%
    filter(!grepl('beta', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_eta = obs$eta) %>%
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

produce.summaries.stan.cohort.rw2 <- function(stan_df, obs, underlying.effects){
  
  summary_fixed <- stan_df[6:7,]
  
  summary_alpha <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('alpha', parameter)) %>%
    filter(!grepl('tau_alpha', parameter)) %>%
    filter(!grepl('alpha_raw', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_alpha = obs$alpha[index])
  
  summary_beta <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('beta', parameter)) %>%
    filter(!grepl('tau_beta', parameter)) %>%
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
    filter(!grepl('kappa_raw', parameter)) %>%
    filter(!grepl('kappa_0', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_kappa = underlying.effects$kappa.true[index]) %>%
    mutate(kappa_drifted = underlying.effects$kappa.drifted[index])
  
  summary_gamma <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('gamma', parameter)) %>%
    filter(!grepl('tau_gamma', parameter)) %>%
    filter(!grepl('gamma_raw', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_gamma = underlying.effects$gamma.true[index]) 
  
  summary_eta <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('eta', parameter)) %>%
    filter(!grepl('beta', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_eta = obs$eta) %>%
    mutate(xt = obs$xt) %>%
    mutate(x = obs$x, t = obs$t)
  
  summaries <- list(
    summary_alpha=summary_alpha,
    summary_beta=summary_beta,
    summary_beta_raw=summary_beta,
    summary_kappa=summary_kappa,
    summary_gamma=summary_gamma,
    summary_eta=summary_eta,
    summary_fixed=summary_fixed
  )
  
  return(summaries)
}

produce.summaries.stan.lc.rw2 <- function(stan_df, obs, underlying.effects){
  
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
    mutate(true_eta = obs$eta) %>%
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

plot.stan.vs.underlying.cohort.undrifted <- function(stan_df, obs, underlying.effects, summaries){
  
  plot_intercept <- ggplot(data=stan_lc_df) +
    geom_vline(aes(xintercept = mean[7], color = "estimated")) + 
    geom_vline(aes(xintercept = `2.5%`[7], color = "estimated"), alpha = 0.5) +
    geom_vline(aes(xintercept = `97.5%`[7], color = "estimated"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$age.intercept.true, color = "true")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) + 
    labs(title="Stan - intercept", x = "Value of intercept", y = " ")
  
  plot_phi <- ggplot(data=stan_lc_df) +
    geom_vline(aes(xintercept = mean[6], color = "estimated")) + 
    geom_vline(aes(xintercept = `2.5%`[6], color = "estimated"), alpha = 0.5) +
    geom_vline(aes(xintercept = `97.5%`[6], color = "estimated"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$phi.true, color = "true")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) + 
    labs(title="Stan - phi", x = "Value of phi", y = " ")
  
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
  
  plot_kappa <- ggplot(data=summary_kappa) +
    geom_point(aes(x=index, y=mean, color="estimated")) + 
    geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
    geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha=0.5) +
    geom_point(aes(x=index, y=true_kappa, color="true")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) + 
    labs(title="Stan - kappa", x = "t", y = " ")
  
  summary_gamma <- summaries$summary_gamma
  
  plot_gamma <- ggplot(data=summary_gamma) +
    geom_point(aes(x=index, y=mean, color="estimated")) + 
    geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
    geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha=0.5) +
    geom_point(aes(x=index, y=true_gamma, color="true")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) + 
    labs(title="Stan - gamma", x = "c", y = " ")
  
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
              p.gamma = plot_gamma,
              p.intercept = plot_intercept,
              p.eta = plot_eta,
              p.phi = plot_phi)
  return(plots)
  
}

plot.stan.vs.underlying.lc.undrifted <- function(stan_df, obs, underlying.effects, summaries){
  
  plot_intercept <- ggplot(data=stan_lc_df) +
    geom_vline(aes(xintercept = mean[6], color = "estimated")) + 
    geom_vline(aes(xintercept = `2.5%`[6], color = "estimated"), alpha = 0.5) +
    geom_vline(aes(xintercept = `97.5%`[6], color = "estimated"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$age.intercept.true, color = "true")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) + 
    labs(title="Stan - intercept", x = "Value of intercept", y = " ")
  
  plot_phi <- ggplot(data=stan_lc_df) +
    geom_vline(aes(xintercept = mean[5], color = "estimated")) + 
    geom_vline(aes(xintercept = `2.5%`[5], color = "estimated"), alpha = 0.5) +
    geom_vline(aes(xintercept = `97.5%`[5], color = "estimated"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$phi.true, color = "true")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) + 
    labs(title="Stan - phi", x = "Value of phi", y = " ")
  
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
                p.eta = plot_eta,
                p.phi = plot_phi)
  return(plots)
  
}

plot.stan.vs.underlying.cohort.drifted <- function(stan_df, obs, underlying.effects, summaries){
  
  plot_intercept <- ggplot(data=stan_lc_df) +
    geom_vline(aes(xintercept = mean[7], color = "estimated")) + 
    geom_vline(aes(xintercept = `2.5%`[7], color = "estimated"), alpha = 0.5) +
    geom_vline(aes(xintercept = `97.5%`[7], color = "estimated"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$age.intercept.true, color = "true")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) + 
    labs(title="Stan - intercept", x = "Value of intercept", y = " ")
  
  plot_phi <- ggplot(data=stan_lc_df) +
    geom_vline(aes(xintercept = mean[6], color = "estimated")) + 
    geom_vline(aes(xintercept = `2.5%`[6], color = "estimated"), alpha = 0.5) +
    geom_vline(aes(xintercept = `97.5%`[6], color = "estimated"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$phi.true, color = "true")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) + 
    labs(title="Stan - phi", x = "Value of phi", y = " ")
  
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
  
  plot_kappa <- ggplot(data=summary_kappa) +
    geom_point(aes(x=index, y=mean, color="estimated")) + 
    geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
    geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha=0.5) +
    geom_point(aes(x=index, y=kappa_drifted, color="true")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) + 
    labs(title="Stan - kappa", x = "t", y = " ")
  
  summary_gamma <- summaries$summary_gamma
  
  plot_gamma <- ggplot(data=summary_gamma) +
    geom_point(aes(x=index, y=mean, color="estimated")) + 
    geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
    geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha=0.5) +
    geom_point(aes(x=index, y=true_gamma, color="true")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) + 
    labs(title="Stan - gamma", x = "c", y = " ")
  
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
                p.gamma = plot_gamma,
                p.intercept = plot_intercept,
                p.eta = plot_eta,
                p.phi = plot_phi)
  return(plots)
  
}

plot.stan.vs.underlying.lc.drifted <- function(stan_df, obs, underlying.effects, summaries){
  
  plot_intercept <- ggplot(data=stan_lc_df) +
    geom_vline(aes(xintercept = mean[6], color = "estimated")) + 
    geom_vline(aes(xintercept = `2.5%`[6], color = "estimated"), alpha = 0.5) +
    geom_vline(aes(xintercept = `97.5%`[6], color = "estimated"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$age.intercept.true, color = "true")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) + 
    labs(title="Stan - intercept", x = "Value of intercept", y = " ")
  
  plot_phi <- ggplot(data=stan_lc_df) +
    geom_vline(aes(xintercept = mean[5], color = "estimated")) + 
    geom_vline(aes(xintercept = `2.5%`[5], color = "estimated"), alpha = 0.5) +
    geom_vline(aes(xintercept = `97.5%`[5], color = "estimated"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$phi.true, color = "true")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) + 
    labs(title="Stan - phi", x = "Value of phi", y = " ")
  
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
  
  plot_kappa <- ggplot(data=summary_kappa) +
    geom_point(aes(x=index, y=mean, color="estimated")) + 
    geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
    geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha=0.5) +
    geom_point(aes(x=index, y=kappa_drifted, color="true")) +
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
                p.eta = plot_eta,
                p.phi = plot_phi)
  return(plots)
  
}

plot.stan.vs.underlying.lc.rw2 <- function(stan_df, obs, underlying.effects, summaries){
  
  plot_intercept <- ggplot(data=stan_lc_df) +
    geom_vline(aes(xintercept = mean[6], color = "estimated")) + 
    geom_vline(aes(xintercept = `2.5%`[6], color = "estimated"), alpha = 0.5) +
    geom_vline(aes(xintercept = `97.5%`[6], color = "estimated"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$age.intercept.true, color = "true")) +
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
  
  plot_kappa <- ggplot(data=summary_kappa) +
    geom_point(aes(x=index, y=mean, color="estimated")) + 
    geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
    geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha=0.5) +
    geom_point(aes(x=index, y=kappa_drifted, color="true")) +
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

plot.stan.vs.underlying.synthetic.cancer.no.beta <- function(stan_df, obs, underlying.effects, summaries){
  
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
  
  # summary_beta <- summaries$summary_beta
  # 
  # plot_beta <- ggplot(data=summary_beta) +
  #   geom_point(aes(x=index, y=mean, color="estimated")) + 
  #   geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
  #   geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha = 0.5) +
  #   geom_point(aes(x=index, y=true_beta, color="true")) +
  #   scale_color_manual(name = " ", values = palette) + 
  #   scale_fill_manual(name = " ", values = palette) + 
  #   labs(title="Stan - beta", x = "x", y = " ")
  
  plot_beta <- ggplot(data.frame(x = 1, y = 2)) + geom_point(aes(x = x, y = y)) + labs(title = "No beta available")
  
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

plot.stan.vs.underlying.cohort.rw2 <- function(stan_df, obs, underlying.effects, summaries){
  
  plot_intercept <- ggplot(data=stan_lc_df) +
    geom_vline(aes(xintercept = mean[5], color = "estimated")) + 
    geom_vline(aes(xintercept = `2.5%`[5], color = "estimated"), alpha = 0.5) +
    geom_vline(aes(xintercept = `97.5%`[5], color = "estimated"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$age.intercept.true, color = "true")) +
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
  
  plot_kappa <- ggplot(data=summary_kappa) +
    geom_point(aes(x=index, y=mean, color="estimated")) + 
    geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
    geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha=0.5) +
    geom_point(aes(x=index, y=kappa_drifted, color="true")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) + 
    labs(title="Stan - kappa", x = "t", y = " ")
  
  summary_gamma <- summaries$summary_gamma
  
  plot_gamma <- ggplot(data=summary_gamma) +
    geom_point(aes(x=index, y=mean, color="estimated")) + 
    geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
    geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha=0.5) +
    geom_point(aes(x=index, y=true_gamma, color="true")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) + 
    labs(title="Stan - gamma", x = "c", y = " ")
  
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
                p.eta = plot_eta,
                p.gamma = plot_gamma)
  return(plots)
  
}

# save.stan.plots.cohort.undrifted <- function(plots, path.to.storage=""){
#   
#   p.alpha <- plots$p.alpha
#   p.beta <- plots$p.beta
#   p.eta <- plots$p.eta
#   p.gamma <- plots$p.gamma
#   p.kappa <- plots$p.kappa
#   p.intercept <- plots$p.intercept
#   p.phi <- plots$p.phi
#   
#   p.random.effects <- (p.intercept | p.alpha | p.beta) / (p.phi | p.kappa | p.gamma) +
#     plot_layout(guides="collect")
#   
#   save.figure(p.random.effects, name = "random_effects_stan", path = path.to.storage)
#   save.figure(p.eta, name = "eta_stan", path = path.to.storage)
# }

# save.stan.plots.lc <- function(plots, path.to.storage=""){
#   
#   p.alpha <- plots$p.alpha
#   p.beta <- plots$p.beta
#   p.eta <- plots$p.eta
#   p.kappa <- plots$p.kappa
#   p.intercept <- plots$p.intercept
#   p.phi <- plots$p.phi
#   
#   p.random.effects <- (p.intercept | p.alpha | p.beta) / (p.phi | p.kappa) +
#     plot_layout(guides="collect")
#   
#   save.figure(p.random.effects, name = "random_effects_stan", path = path.to.storage)
#   save.figure(p.eta, name = "eta_stan", path = path.to.storage)
# }

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

save.stan.plots.cohort.rw2 <- function(plots, path.to.storage=""){
  
  p.alpha <- plots$p.alpha
  p.beta <- plots$p.beta
  p.eta <- plots$p.eta
  p.kappa <- plots$p.kappa
  p.intercept <- plots$p.intercept
  p.gamma <- plots.p.gamma
  
  p.random.effects <- (p.intercept | p.alpha | p.beta) / (p.kappa | p.gamma) +
    plot_layout(guides="collect")
  
  save.figure(p.random.effects, name = "random_effects_stan", path = path.to.storage)
  save.figure(p.eta, name = "eta_stan", path = path.to.storage)
}

produce.stan.plots <- function(stan_df, underlying.effects, plot.func, save.func, 
                               path.to.storage, summaries.func = produce.summaries.stan){
  obs <- underlying.effects$obs
  
  summaries <- summaries.func(stan_df, obs, underlying.effects)
  
  plots <- plot.func(stan_df, obs, underlying.effects, summaries)
  
  save.func(plots, path.to.storage)
  
  res = list(plots = plots, summaries = summaries)
  return(res)
}