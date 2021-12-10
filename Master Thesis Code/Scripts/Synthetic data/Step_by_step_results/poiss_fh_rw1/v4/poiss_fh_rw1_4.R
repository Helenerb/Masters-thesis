#' This program investigates inference with:
#' - The traditional_lc_log_prec version - Gaussian lc-model. 
#' All random effects are modelled as iid, in inlabru and in stan
#' Constraints are imposed as usual in inlabru, and with soft constraints in stan
#' They hyperparameters are fixed. 
#' 

#   ----   Load libraries and set workspace   ----
library("tidyverse")
library("inlabru")
library("ggplot2")
library("INLA")
library("patchwork")
library("rstan")

#TODO: Set this to your local working directory
setwd("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master\ Thesis\ Code")

investigation.name <- "poiss_fh_rw1"

# Path to where files and results are stored
# If you run locally - change this to where you have stored the code
investigation.path <- file.path(investigation.name, "v4")

#   ----    Retrieve the data   ----

synthetic.male.lung.v4 <- function(){
  #TODO: If you run this locally - change to where you have stored the data
  obs <- read.csv("Data/synthetic_male_lung_4.csv")
  
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
# Running traditional lc version of stan, 
# implemented with log-precisions, random effects as iid, and no constraints

# TODO: If you run this locally - change it to where you have your folder containing these files
stan.output  <- file.path("Scripts/Synthetic data/Step_by_step_results", investigation.path)
source("Scripts/Synthetic\ data/run_stan_functions.R")

run_stan <- function(stan_program, obs, chains, warmup, iter, output.path, config.name, markov=TRUE){
  
  stan_fit <- run_stan_program_lc(
    list(obs = obs), chains=chains,warmup=warmup,
    iter=iter, stan_program=stan_program)
  
  store_stan_results(
    fit=stan_fit, output.path=output.path, config=config.name,
    chains=chains, warmup=warmup, iter=iter, stan_program=stan_program)
}

#TODO: if you run locally, change path to where you have stored stan program
run_stan(
  stan_program="Scripts/Synthetic data/Stan analyses/stan_programs/step_by_step_results/stan_pois_fh_rw1_sc_4.stan",
  obs = obs, chains=4, warmup = 8000, iter = 80000, output.path = stan.output,
  config.name = investigation.name, markov=F)

inlabru.pois.fh.rw1 <- function(obs, max_iter=30){
  #'Implements inlabru analysis for lc model, fixing the precisions and modelling all random effects as iid
  #'
  #'@param obs: Contains the observed data and the real underlying random effects
  #'@param max_iter (int): maximum number of iterations in inlabru
  
  nx = length(unique(obs$x))
  nt = length(unique(obs$t))
  
  # constraints for the age effect beta
  A.beta = matrix(1, nrow = 1, ncol = nx)
  e.beta = 1
  
  fixed.theta.alpha <- list(prec = list(initial = log(1.96), fixed = T))
  fixed.theta.beta <- list(prec = list(initial = log(202), fixed = T))
  fixed.theta.kappa <- list(prec = list(initial = log(336), fixed = T))
  fixed.theta.epsilon <- list(prec = list(initial = log(420), fixed = T))
  
  comp = ~ -1 +
    Int(1, prec.linear = 0.001, mean.linear = 0) +
    alpha(x, model = "rw1", hyper = fixed.theta.alpha, constr = TRUE) +
    beta(x.c, model = "iid", hyper = fixed.theta.beta, extraconstr = list(A = A.beta, e = e.beta)) +
    kappa(t, model = "rw1", hyper = fixed.theta.kappa, constr = TRUE) + 
    epsilon(xt, model = "iid", hyper = fixed.theta.epsilon, constr = FALSE)
  
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

res.inlabru <- inlabru.pois.fh.rw1(obs, max_iter = 100)

#  save inlabru object
save(res.inlabru, file = file.path(output.path, "res_inlabru.RData"))

#TODO: If you run these files locally - change sources to the local location
source("Scripts/Functions/plotters.R")
source("Scripts/Synthetic data/plot_inlabru_vs_underlying.R")
source("Scripts/Synthetic data/plot_inlabru_stan_compared.R")
source("Scripts/Synthetic data/plot_stan_vs_underlying.R")

output.path <- stan.output

plots.summaries.inlabru <- plot.inlabru.vs.underlying.synthetic.cancer.fixed.effects(
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
                               save.func=function(...) {save.stan.plots.lc.rw2(..., save=F)},
                               path.to.storage=output.path,
                               summaries.func=produce.summaries.stan.poiss.lc)

plots_compared <- produce.compared.plots(
  stan.summaries = stan.res$summaries,
  stan.marginals = stan.marginals,
  inlabru.summaries = plots.summaries.inlabru$summaries,
  res.inlabru = res.inlabru,
  underlying.effects = underlying.effects,
  plot.func = function(...) {plot.inlabru.stan.compared.poisson.lc.fixed.hypers(..., tau.beta.cutoff = 700, tau.kappa.cutoff = 500, tau.alpha.cutoff = 100)},
  #plot.func = function(...) {plot.inlabru.stan.traditional.lc(..., cohort=FALSE, tau.beta.cutoff = 700, tau.kappa.cutoff = 500, tau.alpha.cutoff = 10, a45=F)},
  #plot.func = function(...) {plot.inlabru.stan.traditional.lc.fixed.hypers(..., cohort=FALSE, tau.beta.cutoff = 700, tau.kappa.cutoff = 500, tau.alpha.cutoff = 10, a45=F)},
  #plot.func = function(...) {plot.inlabru.stan.traditional.lc.no.beta(..., cohort=FALSE, tau.beta.cutoff = 5000, tau.kappa.cutoff = 5000, tau.alpha.cutoff = 100, a45=F)},
  #plot.func = function(...) {plot.inlabru.stan.traditional.lc.fixed.hypers.no.beta(..., cohort=FALSE, tau.beta.cutoff = 5000, tau.kappa.cutoff = 5000, tau.alpha.cutoff = 100, a45=F)},
  save.func = function(...) {save.compared.rw2(..., cohort=FALSE, png = F)},
  path.to.storage=output.path)

#   ----   Sample predictor   ----

stan.predictor.df <- data.frame(eta_draws)

plot.predictor.inlabru.stan.compared(res.inlabru, stan.predictor.df, path.to.storage = output.path, a45=F)

#   ----   Plot marginals of random effects   ----

stan.beta.df <- data.frame(beta_draws)
plot.beta.inlabru.stan.compared(res.inlabru, stan.beta.df, path.to.storage = output.path, a45=F)

stan.kappa.df <- data.frame(kappa_draws)
plot.kappa.inlabru.stan.compared(res.inlabru, stan.kappa.df, path.to.storage = output.path)

#   ----   Specifically check the predictors at xt = 54:   ----

pred.54.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.054)

p.pred.54 <- ggplot(pred.54.inlabru) + 
  geom_area(aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
  geom_histogram(data = stan.predictor.df, aes(x = X54, y = after_stat(density), fill = "Stan", color = "Stan"), alpha = 0.5, bins = 100) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Predictor at xt=54", x = "", y = "")
p.pred.54

save.figure(p.pred.54, name = "predictor_54", path = output.path, png= F)

#   ----   Specifically save the first values of the predictor   ----

pred.1.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.001)
p.pred.1 <- ggplot(pred.1.inlabru) + 
  geom_area(aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
  geom_histogram(data = stan.predictor.df, aes(x = X1, y = after_stat(density), fill = "Stan", color = "Stan"), alpha = 0.5, bins = 50) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Pred[1]", x = "", y = "")

pred.3.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.003)
p.pred.3 <- ggplot(pred.3.inlabru) + 
  geom_area(aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
  geom_histogram(data = stan.predictor.df, aes(x = X3, y = after_stat(density), fill = "Stan", color = "Stan"), alpha = 0.5, bins = 50) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Pred[3]", x = "", y = "")

pred.5.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.005)
p.pred.5 <- ggplot(pred.5.inlabru) + 
  geom_area(aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
  geom_histogram(data = stan.predictor.df, aes(x = X5, y = after_stat(density), fill = "Stan", color = "Stan"), alpha = 0.5, bins = 50) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Pred[5]", x = "", y = "")

pred.8.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.008)
p.pred.8 <- ggplot(pred.8.inlabru) + 
  geom_area(aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
  geom_histogram(data = stan.predictor.df, aes(x = X8, y = after_stat(density), fill = "Stan", color = "Stan"), alpha = 0.5, bins = 50) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Pred[8]", x = "", y = "")

pred.11.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.011)
p.pred.11 <- ggplot(pred.11.inlabru) + 
  geom_area(aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
  geom_histogram(data = stan.predictor.df, aes(x = X11, y = after_stat(density), fill = "Stan", color = "Stan"), alpha = 0.5, bins = 50) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Pred[11]", x = "", y = "")

pred.14.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.014)
p.pred.14 <- ggplot(pred.14.inlabru) + 
  geom_area(aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
  geom_histogram(data = stan.predictor.df, aes(x = X14, y = after_stat(density), fill = "Stan", color = "Stan"), alpha = 0.5, bins = 50) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Pred[14]", x = "", y = "")

pred.18.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.018)
p.pred.18 <- ggplot(pred.18.inlabru) + 
  geom_area(aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
  geom_histogram(data = stan.predictor.df, aes(x = X18, y = after_stat(density), fill = "Stan", color = "Stan"), alpha = 0.5, bins = 50) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Pred[18]", x = "", y = "")

pred.22.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.022)
p.pred.22 <- ggplot(pred.22.inlabru) + 
  geom_area(aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
  geom_histogram(data = stan.predictor.df, aes(x = X22, y = after_stat(density), fill = "Stan", color = "Stan"), alpha = 0.5, bins = 50) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Pred[22]", x = "", y = "")

pred.26.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.026)
p.pred.26 <- ggplot(pred.26.inlabru) + 
  geom_area(aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
  geom_histogram(data = stan.predictor.df, aes(x = X26, y = after_stat(density), fill = "Stan", color = "Stan"), alpha = 0.5, bins = 50) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Pred[26]", x = "", y = "")

pred.30.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.030)
p.pred.30 <- ggplot(pred.30.inlabru) + 
  geom_area(aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
  geom_histogram(data = stan.predictor.df, aes(x = X30, y = after_stat(density), fill = "Stan", color = "Stan"), alpha = 0.5, bins = 50) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Pred[30]", x = "", y = "")

pred.34.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.034)
p.pred.34 <- ggplot(pred.34.inlabru) + 
  geom_area(aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
  geom_histogram(data = stan.predictor.df, aes(x = X34, y = after_stat(density), fill = "Stan", color = "Stan"), alpha = 0.5, bins = 50) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Pred[34]", x = "", y = "")

pred.38.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.038)
p.pred.38 <- ggplot(pred.38.inlabru) + 
  geom_area(aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
  geom_histogram(data = stan.predictor.df, aes(x = X38, y = after_stat(density), fill = "Stan", color = "Stan"), alpha = 0.5, bins = 50) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Pred[38]", x = "", y = "")

p.predictor.start <- (p.pred.1 | p.pred.3 | p.pred.5 | p.pred.8)/(p.pred.11 | p.pred.14 | p.pred.18 | p.pred.22) / (p.pred.26 | p.pred.30 | p.pred.34 | p.pred.38) + 
  plot_layout(guides = "collect")
p.predictor.start

save.figure(p.predictor.start, name = "predictor_start", path = output.path, png = F)




