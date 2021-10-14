# TODO: Script running and comparing inference with the lc-cohort model
# by inlabru and STAN. 

setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data")
source("configurations_synthetic_data.R")

# cohfiguration with cohort effect
#underlying.effects.lc.cohort <- configuration.v7()  ##  Fine grid
#underlying.effects.lc.cohort <- configuration.v17()  ## First attempt at coarser grid
underlying.effects.lc.cohort <- configuration.v18()  ## more erratic beta
#underlying.effects.lc.cohort <- configuration.v19()
#underlying.effects.lc.cohort <- configuration.v20()
#underlying.effects.lc.cohort <- configuration.v21()
#underlying.effects.lc.cohort <- configuration.v22()
#underlying.effects.lc.cohort <- configuration.v23()

figures.folder = "/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Inlabru analyses/Output/Figures"

#storage_path = file.path(figures.folder, "v7")
#storage_path = file.path(figures.folder, "v17")
storage_path = file.path(figures.folder, "v18")
#storage_path = file.path(figures.folder, "v19")
#storage_path = file.path(figures.folder, "v20")
#storage_path = file.path(figures.folder, "v21")
#storage_path = file.path(figures.folder, "v22")
#storage_path = file.path(figures.folder, "v23")

obs.cohort <- underlying.effects.lc.cohort$obs

source("inlabru_analyses.R")
runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.lc.cohort.kappa_high_variance_prior(obs.cohort)})


###   ----   Plot results from inlabru inference   ----  

source("plot_inlabru_vs_underlying.R")

#plot.period.posterior <- plot.posterior.period.effects(res.inlabru.lc.1, underlying.effects.lc.cohort)

# plotting results from run with cohort effects:
plots.lc.cohort <- plot.inlabru.vs.underlying.cohort(
  res.inlabru.lc.1, 
  underlying.effects.lc.cohort,path.to.storage = storage_path,
  save=TRUE,
  phi.plus.kappa.func = phi.plus.kappa.v17)
plots.lc.cohort$p.alpha
plots.lc.cohort$p.beta
plots.lc.cohort$p.phi
plots.lc.cohort$p.intercept
plots.lc.cohort$p.kappa
plots.lc.cohort$p.eta
plots.lc.cohort$p.eta.2
plots.lc.cohort$p.eta.t
plots.lc.cohort$p.eta.x
plots.lc.cohort$p.gamma

###    ----   Configure and run inference with STAN    ----

#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v17/stan_v17.Rda")
load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v18/stan_v18.Rda")

plot.stan.vs.underlying.cohort <- function(stan_df, obs, underlying.effects){
  
  plot_intercept <- ggplot(data=stan_lc_df) +
    geom_vline(aes(xintercept = mean[7], color = "estimated")) + 
    geom_vline(aes(xintercept = `2.5%`[7], color = "estimated"), alpha = 0.5) +
    geom_vline(aes(xintercept = `97.5%`[7], color = "estimated"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$age.intercept.true, color = "true"))
  
  plot_phi <- ggplot(data=stan_lc_df) +
    geom_vline(aes(xintercept = mean[6], color = "estimated")) + 
    geom_vline(aes(xintercept = `2.5%`[6], color = "estimated"), alpha = 0.5) +
    geom_vline(aes(xintercept = `97.5%`[6], color = "estimated"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$phi.true, color = "true"))
  
  summary_alpha <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('alpha', parameter)) %>%
    filter(!grepl('tau_alpha', parameter)) %>%
    filter(!grepl('alpha_raw', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_alpha = obs$alpha[index])
  
  plot_alpha <- ggplot(data=summary_alpha) +
    geom_point(aes(x=index, y=mean, color="estimated")) + 
    geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
    geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha=0.5) +
    geom_point(aes(x=index, y=true_alpha, color="true value")) +
    ggtitle("Alpha")
  
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
  
  plot_beta <- ggplot(data=summary_beta) +
    geom_point(aes(x=index, y=mean, color="estimated")) + 
    geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
    geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha = 0.5) +
    geom_point(aes(x=index, y=true_beta, color="true value")) +
    ggtitle("Beta -STAN")
  
  summary_kappa <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('kappa', parameter)) %>%
    filter(!grepl('tau_kappa', parameter)) %>%
    filter(!grepl('kappa_raw', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_kappa = underlying.effects$kappa.true[index]) %>%
    mutate(kappa_drifted = true_kappa + underlying.effects$phi.true*(index-1))
  
  plot_kappa <- ggplot(data=summary_kappa) +
    geom_point(aes(x=index, y=mean, color="estimated")) + 
    geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
    geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha=0.5) +
    geom_point(aes(x=index, y=true_kappa, color="true undrifted")) +
    geom_point(aes(x=index, y=kappa_drifted, color="true drifted")) +
    ggtitle("Kappa - STAN")
  
  summary_gamma <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('gamma', parameter)) %>%
    filter(!grepl('tau_gamma', parameter)) %>%
    filter(!grepl('gamma_raw', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_gamma = underlying.effects$gamma.true[index])  # kanskje du må skru den opp...
  
  plot_gamma <- ggplot(data=summary_gamma) +
    geom_point(aes(x=index, y=mean, color="estimated")) + 
    geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
    geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha=0.5) +
    geom_point(aes(x=index, y=true_gamma, color="true value")) +
    ggtitle("Gamma - STAN")
  
  summary_eta <- stan_df %>%
    rownames_to_column("parameter") %>%
    filter(grepl('eta', parameter)) %>%
    filter(!grepl('beta', parameter)) %>%
    mutate(index = parse_number(parameter)) %>%
    mutate(true_eta = obs$eta)
  
  plot_eta <- ggplot(data=summary_eta) +
    geom_line(aes(x=index, y=mean, color="estimated")) + 
    geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
    geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha=0.5) +
    geom_line(aes(x=index, y=true_eta, color="true value"), alpha = 0.5) +
    ggtitle("Eta - STAN")
  
  res <- list(p.alpha = plot_alpha,
              p.beta = plot_beta, 
              p.kappa = plot_kappa,
              p.gamma = plot_gamma,
              p.intercept = plot_intercept,
              p.eta = plot_eta,
              p.phi = plot_phi)
}

save.stan.plots <- function(plots, path.to.storage=""){
  
  p.alpha <- plots$p.alpha
  p.beta <- plots$p.beta
  p.eta <- plots$p.eta
  p.gamma <- plots$p.gamma
  p.kappa <- plots$p.kappa
  p.intercept <- plots$p.intercept
  p.phi <- plots$p.phi
  
  p.random.effects <- (p.intercept | p.alpha | p.beta) / (p.phi | p.kappa | p.gamma) +
    plot_layout(guides="collect")
  
  save.figure(p.random.effects, name = "random_effects_stan", path = path.to.storage)
  save.figure(p.eta, name = "eta_stan", path = path.to.storage)
}

stan.plots <- plot.stan.vs.underlying.cohort(stan_df = stan_lc_df, obs=obs.cohort, underlying.effects = underlying.effects.lc.cohort)
save.stan.plots(plots=stan.plots, path.to.storage=storage_path)
stan.plots$p.alpha
stan.plots$p.beta
stan.plots$p.eta
stan.plots$p.gamma
stan.plots$p.kappa
