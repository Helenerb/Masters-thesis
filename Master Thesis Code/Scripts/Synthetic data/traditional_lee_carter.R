# script with traditional lc-model, data based on the 
# synthetic male lung configuration, v4

library("tidyverse")

setwd("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master\ Thesis\ Code")

synthetic.male.lung.v4 <- function(){
  obs <- read.csv("Data/synthetic_male_lung_4.csv")
  
  underlying.effects <- list(obs = obs, nx = 9, nt = 18,
                             alpha.true = {obs %>% filter(t == 0)}$alpha,
                             beta.true = {obs %>% filter(t == 0)}$beta,
                             kappa.true = {obs %>% filter(x == 9)}$kappa,
                             intercept = unique(obs$intercept),
                             tau.alpha.true = unique(obs$tau.alpha),
                             tau.beta.true = unique(obs$tau.beta),
                             tau.kappa.true = unique(obs$tau.kappa),
                             tau.epsilon.true = unique(obs$tau.epsilon))
  return(list(obs = obs, underlying.effects = underlying.effects))
}

# TODO: Move this to inlabru_analyses.R when this runs ok
inlabru.traditional.lc <- function(obs, max_iter=30){
  #' TODO: Implement traditional lee-carter
  #'Implements inlabru analysis for lc model using an ar1c to model the period effect
  #'
  #'@param obs: Contains the observed data and the real underlying random effects
  #'@param max_iter (int): maximum number of iterations in inlabru
  
  nx = length(unique(obs$x))
  nt = length(unique(obs$t))
  
  # constraints for the age effect beta
  A.beta = matrix(1, nrow = 1, ncol = nx)  
  e.beta = 1  
  
  #pc.prior <- list(prec = list(prior = "pc.prec", param = c(1,0.05)))
  loggamma.prior <- list(prec = list(prior = 'loggamma', param = c(1,0.00005), initial = log(1)))
  loggamma.prior.high.variance <- list(prec = list(prior = 'loggamma', param = c(1,0.005), initial = log(1)))
  
  comp = ~ -1 +
    Int(1, prec.linear = 0.001, mean.linear = 0) +
    alpha(x, model = "rw1", values=unique(obs$x), hyper = loggamma.prior, constr = TRUE) +
    beta(x.c, model = "iid", extraconstr = list(A = A.beta, e = e.beta), hyper = loggamma.prior) +
    kappa(t, model = "rw2", values = unique(obs$t), constr = TRUE, hyper = loggamma.prior.high.variance)
  
  formula = eta ~ Int + alpha + beta*kappa 
  
  likelihood = like(formula = formula, family = "gaussian", data = obs)
  
  c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)  # control.compute
  c.family <- list(hyper = list(
    theta = list(
      initial = log(1), fixed = F, prior = "loggamma", param = c(1,0.00005)
    )
  ))
  #c.fixed <- list(mean = list(Int = 0), prec = list(Int = 1))
  
  res.inlabru = bru(components = comp,
                    likelihood, 
                    options = list(verbose = F,
                                   bru_verbose = 1, 
                                   num.threads = "1:1",
                                   control.compute = c.c,
                                   bru_max_iter=max_iter,
                                   # control.predictor = list(link = 1), # Dont know if this still applies
                                   control.family  = c.family
                    ))
  return(res.inlabru)
}

config_data <- synthetic.male.lung.v4()

underlying.effects.lc.poiss <- config_data$underlying.effects

obs.lc.poiss <- underlying.effects.lc.poiss$obs

obs.trad <- obs.lc.poiss %>% 
  select(c(x, t, xt, age.int, year, x.c, alpha, beta, kappa, intercept, epsilon,
           eta, tau.alpha, tau.beta, tau.kappa, tau.epsilon)) %>%
  mutate(mu = exp(eta)) %>%
  mutate(eta.no.error = intercept + alpha + beta*kappa)

underlying.effects.trad <- list(obs = obs.trad, nx = 18, nt = 18,
                                alpha.true = underlying.effects.lc.poiss$alpha.true,
                                beta.true = underlying.effects.lc.poiss$beta.true,
                                kappa.true = underlying.effects.lc.poiss$kappa.true,
                                intercept = underlying.effects.lc.poiss$intercept,
                                age.intercept.true = underlying.effects.lc.poiss$intercept,
                                tau.alpha.true = underlying.effects.lc.poiss$tau.alpha.true,
                                tau.beta.true = underlying.effects.lc.poiss$tau.beta.true,
                                tau.kappa.true = underlying.effects.lc.poiss$tau.kappa.true,
                                tau.epsilon.true = underlying.effects.lc.poiss$tau.epsilon.true)

runtime.inlabru <- system.time({res.inlabru.trad <- inlabru.traditional.lc(obs.trad)})

#    ----   plot inlabru results   ----
ggplot(data = data.frame(res.inlabru.trad$marginals.fixed)) + 
  geom_area(aes(x = Int.x, y = Int.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_vline(aes(xintercept = unique(obs.trad$intercept), color = "True", fill = "True"))

ggplot(data = data.frame(res.inlabru.trad$summary.random$alpha)) + 
  geom_ribbon(aes(x = ID, ymin = X0.025quant, ymax = X0.975quant, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_point(aes(x = ID, y = mean, color = "Inlabru", fill = "Inlabru")) + 
  geom_point(aes(x = ID, y = underlying.effects.lc.poiss$alpha.true[1:18], color = "True", fill = "True"))

ggplot(data = data.frame(res.inlabru.trad$summary.random$beta)) + 
  geom_errorbar(aes(x = ID, ymin = X0.025quant, ymax = X0.975quant, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_point(aes(x = ID, y = mean, color = "Inlabru", fill = "Inlabru")) + 
  geom_point(aes(x = ID, y = underlying.effects.lc.poiss$beta.true[1:18], color = "True", fill = "True"))

ggplot(data = data.frame(res.inlabru.trad$summary.random$kappa)) + 
  geom_errorbar(aes(x = ID, ymin = X0.025quant, ymax = X0.975quant, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_point(aes(x = ID, y = mean, color = "Inlabru", fill = "Inlabru")) + 
  geom_point(aes(x = ID, y = underlying.effects.lc.poiss$kappa.true[1:18], color = "True", fill = "True"))

ggplot(data = data.frame(mean = res.inlabru.trad$summary.fitted.values$mean[1:324],
                         X0.025quant = res.inlabru.trad$summary.fitted.values$`0.025quant`[1:324],
                         X0.975quant = res.inlabru.trad$summary.fitted.values$`0.975quant`[1:324]) %>%
         mutate(ID = seq_along(mean))) + 
  geom_point(data = obs.trad, aes(x = xt + 1, y = eta.no.error, color = "True", fill  ="True")) + 
  geom_point(aes(x = ID, y = mean, color = "Inlabru", fill = "Inlabru")) + 
  geom_errorbar(aes(x = ID, ymin = X0.025quant, ymax = X0.975quant, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) 

#   ----    run stan analysis   ----

# library("rstan")
# 
# input_stan <- list(
#   X=length(unique(obs.trad$x)),
#   T=length(unique(obs.trad$t)),
#   x=(obs.trad$x + 1),
#   t=(obs.trad$t + 1),
#   exp_mortality_rate = obs.trad$eta,
#   nx = length(unique(obs.trad$x)),
#   nt = length(unique(obs.trad$t))
# )
# 
# 
# fit <- stan(
#   file="Scripts/Synthetic\ data/Stan\ analyses/stan_programs/stan_lc_traditional.stan",
#   data = input_stan,
#   chains=2,
#   warmup = 100,
#   iter = 400,
#   refresh = 400%/%10,
#   seed=123
# )
# 
# return(fit)

#  ----   Compare stan and inlabru   ----
source("Scripts/Functions/plotters.R")
source("Scripts/Synthetic data/plot_inlabru_vs_underlying.R")
source("Scripts/Synthetic data/plot_inlabru_stan_compared.R")
source("Scripts/Synthetic data/plot_stan_vs_underlying.R")

figures.folder = "Scripts/Synthetic data/Output/Figures"
storage_path = file.path(figures.folder, "traditional_lc")

# only inlabru
#plots.summaries.inlabru <- plot.inlabru.vs.underlying.synthetic.cancer(
plots.summaries.inlabru <- plot.inlabru.vs.underlying.traditional.lc(
  res.inlabru.trad,
  underlying.effects.trad,
  path.to.storage = storage_path,
  save=TRUE, cutoff_alpha = 500)
  
# when results are produced locally
load("Scripts/Synthetic data/Stan analyses/traditional_lc/stan_traditional_lc.Rda")
path.to.stan.results = "Scripts/Synthetic\ data/Stan analyses/traditional_lc"


load(file=file.path(path.to.stan.results, "draws_intercept.RData"))
load(file=file.path(path.to.stan.results, "draws_tau_epsilon.RData"))
load(file.path(path.to.stan.results, "draws_tau_alpha.RData"))
load(file.path(path.to.stan.results, "draws_tau_beta.RData"))
load(file.path(path.to.stan.results, "draws_tau_kappa.RData"))
load(file.path(path.to.stan.results, "draws_alpha.RData"))
load(file.path(path.to.stan.results, "draws_beta.RData"))
load(file.path(path.to.stan.results, "draws_kappa.RData"))
load(file.path(path.to.stan.results, "draws_eta_100.RData"))
load(file.path(path.to.stan.results, "draws_eta.RData"))
load(file.path(path.to.stan.results, "draws_eta_reduced.RData"))

stan.marginals <- list(intercept_draws = intercept_draws,
                       tau_epsilon_draws = tau_epsilon_draws,
                       tau_alpha_draws = tau_alpha_draws,
                       tau_beta_draws = tau_beta_draws,
                       tau_kappa_draws = tau_kappa_draws,
                       alpha_draws = alpha_draws,
                       beta_draws = beta_draws,
                       kappa_draws = kappa_draws,
                       eta_draws = eta_draws_reduced)

stan.res <- produce.stan.plots(stan_df=stan_lc_df,
                               underlying.effects=underlying.effects.trad,
                               plot.func=plot.stan.vs.underlying.synthetic.cancer,
                               save.func=save.stan.plots.lc.rw2,
                               path.to.storage=storage_path,
                               summaries.func=produce.summaries.stan.lc.rw2)


plots_compared <- produce.compared.plots(
  stan.summaries = stan.res$summaries,
  stan.marginals = stan.marginals,
  inlabru.summaries = plots.summaries.inlabru$summaries,
  res.inlabru = res.inlabru.trad,
  underlying.effects = underlying.effects.trad,
  plot.func = function(...) {plot.inlabru.stan.traditional.lc(..., cohort=FALSE, tau.beta.cutoff = 5000, tau.kappa.cutoff = 5000, tau.alpha.cutoff = 100)},
  save.func = function(...) {save.compared.rw2(..., cohort=FALSE)},
  path.to.storage=storage_path)
  
  
  