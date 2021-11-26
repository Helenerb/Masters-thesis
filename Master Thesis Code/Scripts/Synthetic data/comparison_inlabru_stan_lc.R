library(patchwork)

setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code")
source("Scripts/Synthetic data/configurations_synthetic_data.R")

# configuration without cohort effect
#underlying.effects.lc <- configuration.v5()
#underlying.effects.lc <- configuration.v9()  #  config with coarser grid
#underlying.effects.lc <- configuration.v10()
#underlying.effects.lc <- configuration.v10.1()# config with coarser grid and larger variance in beta
#underlying.effects.lc <- configuration.v10.2()
#underlying.effects.lc <- configuration.v10.3()
#underlying.effects.lc <- configuration.v11()  # config with coarser grid, larger variance in beta and steeper phi (compared to alpha)
#underlying.effects.lc <- configuration.v11.1()
#underlying.effects.lc <- configuration.v11.3()
#underlying.effects.lc <- configuration.v12()  #  config with coarser grid, smaller variance in beta, steeper phi (compared to alpha)
#underlying.effects.lc <- configuration.v12.3()
#underlying.effects.lc <- configuration.v13()  # config with coarser grid, even smaller variance in beta, a bit less steep phi, higher variance in kappa
#underlying.effects.lc <- configuration.v14()
#underlying.effects.lc <- configuration.v15()  # half the grid of the original v5

#source("Scripts/Real\ data/synthetic_male_stomach_lc.R")
#underlying.effects.lc <- synthetic.male.stomach.lc()

# source("Scripts/Synthetic\ data/config_synthetic_male_lung_v6.R")
# config_data <- synthetic.male.lung.v6()

source("Scripts/Synthetic\ data/config_synthetic_male_lung_v4.R")
config_data <- synthetic.male.lung.v4()

#source("Scripts/Synthetic\ data/config_synthetic_male_lung_v7.R")
#config_data <- synthetic.male.lung.a45.v7()

underlying.effects.lc <- config_data$underlying.effects


#underlying.effects.lc <- configuration.test.1()

figures.folder = "/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Output/Figures"

#storage_path = file.path(figures.folder, "v10_2")
#storage_path = file.path(figures.folder, "v10_2_only_kappa")
#storage_path = file.path(figures.folder, "v10_3_only_kappa")
#storage_path = file.path(figures.folder, "v10d")
#storage_path = file.path(figures.folder, "v10dh")
#storage_path = file.path(figures.folder, "v10_2_ar1c")
#storage_path = file.path(figures.folder, "v10_3_ar1c")
#storage_path = file.path(figures.folder, "v10_3_rw2")
#storage_path = file.path(figures.folder, "v11_3_rw2")
#storage_path = file.path(figures.folder, "v12_3_rw2")

#storage_path = file.path(figures.folder, "synthetic_male_stomach_lc")

#storage_path = file.path(figures.folder, "synthetic_male_lung_lc/v6")
#storage_path = file.path(figures.folder, "synthetic_male_lung_lc/v4")
#storage_path = file.path(figures.folder, "synthetic_male_lung_lc/v4/rw1")
storage_path = file.path(figures.folder, "synthetic_male_lung_lc/v4/rw1/soft_constraints_local")
#storage_path = file.path(figures.folder, "synthetic_male_lung_lc/v7/rw1")

obs.lc <- underlying.effects.lc$obs


source("Scripts/Functions/inlabru_analyses.R")
#runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.lc.kappa_high_variance_prior(obs.lc)})
#runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.undrifted.period.2(obs.lc)})
#runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.ar1c.lc(obs.lc)})
#runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.ar1c.lc.2(obs.lc)})
#runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.rw2.lc.2(obs.lc)})

#runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.rw2.lc.2(obs.lc)})
runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.rw1.lc(obs.lc)})

#res.inlabru.no.int <- inlabru.rw2.lc.no.intercept(obs.lc, max_iter = 100)

source("Scripts/Synthetic data/plot_inlabru_vs_underlying.R")


# plots.summaries.inlabru <- plot.inlabru.vs.underlying.lc.only.kappa.2(
#   res.inlabru.lc.1,
#   underlying.effects.lc,
#   path.to.storage = storage_path,
#   save=TRUE,
#   phi.plus.kappa.func = phi.plus.kappa.v17)

# plots.inlabru.no.int <- plot.inlabru.vs.underlying.lc.only.kappa.no.intercept(
#   res.inlabru.no.int,
#   underlying.effects.lc,
#   path.to.storage = storage_path,
#   save = TRUE
# )

plots.summaries.inlabru <- plot.inlabru.vs.underlying.synthetic.cancer(
  res.inlabru.lc.1,
  underlying.effects.lc,
  path.to.storage = storage_path,
  save=TRUE)


print("Runtime for inlabru: ")
print(runtime.inlabru)


#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v10/stan_results/stan_v10.Rda")
#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v10d/stan_results/stan_v10d.Rda")
#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v10dh/stan_results/stan_v10dh.Rda")
#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v10_3/stan_results/stan_v10_3.Rda")
#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v11_3/stan_results/stan_v11_3.Rda")
#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v12_3/stan_results/stan_v12_3.Rda")

#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/synthetic_male_stomach_lc/stan_results/stan_synthetic_male_stomach_lc.Rda")

#load("Scripts/Synthetic data/Stan analyses/synthetic_male_lung_6/stan_results/stan_synthetic_male_lung_6.Rda")
#load("Scripts/Synthetic data/Stan analyses/synthetic_male_lung_4/stan_results/stan_synthetic_male_lung_4.Rda")
#load("Scripts/Synthetic data/Stan analyses/synthetic_male_lung_7/stan_results/stan_synthetic_male_lung_7.Rda")
#load("Scripts/Synthetic data/Stan analyses/synthetic_male_lung_7/stan_synthetic_male_lung_7.Rda")
load("Scripts/Synthetic data/Stan analyses/synthetic_male_lung_4/rw1/stan_synthetic_male_lung_4.Rda")


#   ----   load STAN marginals   ---- 

#path.to.stan.results = "/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master\ Thesis\ Code/Scripts/Synthetic\ data/Stan analyses/v10_3/stan_results"
#path.to.stan.results = "/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master\ Thesis\ Code/Scripts/Synthetic\ data/Stan analyses/v12_3/stan_results"

#path.to.stan.results = "Scripts/Synthetic\ data/Stan analyses/synthetic_male_stomach_lc/stan_results"
#path.to.stan.results = "Scripts/Synthetic\ data/Stan analyses/synthetic_male_lung_6/stan_results"
#path.to.stan.results = "Scripts/Synthetic\ data/Stan analyses/synthetic_male_lung_4/stan_results"
#path.to.stan.results = "Scripts/Synthetic\ data/Stan analyses/synthetic_male_lung_7/stan_results"
#path.to.stan.results = "Scripts/Synthetic\ data/Stan analyses/synthetic_male_lung_7"
path.to.stan.results = "Scripts/Synthetic\ data/Stan analyses/synthetic_male_lung_4/rw1"


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

source("Scripts/Functions/plotters.R")

#   ----   trace plots for stan results   ----
trace.intercept <- trace_plot(intercept_draws, chains = 1, iterations = 3400000, warmup = 340000, title= "Trace plot intercept - lsynthetic male lung cancer")
save.figure(trace.intercept, name="trace_intercept", path=storage_path, pdf=F)

trace.tau.beta <- trace_plot(tau_beta_draws, chains = 1, iterations = 3400000, warmup = 340000, title= "Trace plot tau beta - synthetic male lung cancer")
save.figure(trace.tau.beta, name="trace_tau_beta", path=storage_path, pdf=F)

trace.tau.beta.first <- trace_plot(tau_beta_draws[1:5000], chains = 1, iterations = 5000, warmup = 0, title = "Trace plot tau beta, first 5000 - synthetic male lung cancer")
save.figure(trace.tau.beta.first, name = "trace_tau_beta_first", path = storage_path, pdf = F)

trace.tau.beta.last <- trace_plot(tau_beta_draws[3055001:3060000], chains = 1, iterations = 5000, warmup = 0, title = "Trace plot tau beta, last 5000 - synthetic male lung cancer")
save.figure(trace.tau.beta.last, name = "trace_tau_beta_last", path = storage_path, pdf = F)


trace.tau.kappa <- trace_plot(tau_kappa_draws, chains = 1, iterations = 3400000, warmup = 340000, title= "Trace plot tau kappa - synthetic male lung cancer")
save.figure(trace.tau.kappa, name="trace_tau_kappa", path=storage_path, pdf=F)

trace.tau.kappa.first <- trace_plot(tau_kappa_draws[1:5000], chains = 1, iterations = 5000, warmup = 0, title = "Trace plot tau kappa, first 5000 - synthetic male lung cancer")
save.figure(trace.tau.kappa.first, name = "trace_tau_kappa_first", path = storage_path, pdf = F)

trace.tau.kappa.last <- trace_plot(tau_kappa_draws[3055001:3060000], chains = 1, iterations = 5000, warmup = 0, title = "Trace plot tau kappa, last 5000 - synthetic male lung cancer")
save.figure(trace.tau.kappa.last, name = "trace_tau_kappa_last", path = storage_path, pdf = F)


trace.tau.alpha <- trace_plot(tau_alpha_draws, chains = 1, iterations = 3400000, warmup = 340000, title= "Trace plot tau alpha - synthetic male lung cancer")
save.figure(trace.tau.alpha, name="trace_tau_alpha", path=storage_path, pdf=F)

trace.tau.alpha.first <- trace_plot(tau_alpha_draws[1:5000], chains = 1, iterations = 5000, warmup = 0, title = "Trace plot tau alpha, first 5000 - synthetic male lung cancer")
save.figure(trace.tau.alpha.first, name = "trace_tau_alpha_first", path = storage_path, pdf = F)

trace.tau.alpha.last <- trace_plot(tau_alpha_draws[3055001:3060000], chains = 1, iterations = 5000, warmup = 0, title = "Trace plot tau alpha, last 5000 - synthetic male lung cancer")
save.figure(trace.tau.alpha.last, name = "trace_tau_alpha_last", path = storage_path, pdf = F)


trace.tau.epsilon <- trace_plot(tau_epsilon_draws, chains = 1, iterations = 3400000, warmup = 340000, title= "Trace plot tau epsilon - synthetic male lung cancer")
save.figure(trace.tau.epsilon, name="trace_tau_epsilon", path=storage_path, pdf=F)

trace.tau.epsilon.first <- trace_plot(tau_epsilon_draws[1:5000], chains = 1, iterations = 5000, warmup = 0, title = "Trace plot tau epsilon, first 5000 - synthetic male lung cancer")
save.figure(trace.tau.epsilon.first, name = "trace_tau_epsilon_first", path = storage_path, pdf = F)

trace.tau.epsilon.last <- trace_plot(tau_epsilon_draws[3055001:3060000], chains = 1, iterations = 5000, warmup = 0, title = "Trace plot tau epsilon, last 5000 - synthetic male lung cancer")
save.figure(trace.tau.epsilon.last, name = "trace_tau_epsilon_last", path = storage_path, pdf = F)



trace.beta.1 <- trace_plot(beta_draws[,1], chains = 1, iterations = 3400000, warmup = 340000, title= "Trace plot beta[1] - synthetic male lung cancer")
save.figure(trace.beta.1, name="trace_beta_1", path=storage_path, pdf=F)

trace.beta.1.first <- trace_plot(beta_draws[1:5000,1], chains = 1, iterations = 5000, warmup = 0, title = "Trace plot beta[1], first 5000 - synthetic male lung cancer")
save.figure(trace.beta.1.first, name = "trace_beta_1_first", path = storage_path, pdf = F)

trace.beta.1.last <- trace_plot(beta_draws[3055001:3060000,1], chains = 1, iterations = 5000, warmup = 0, title = "Trace plot beta[1], last 5000 - synthetic male lung cancer")
save.figure(trace.beta.1.last, name = "trace_beta_1_last", path = storage_path, pdf = F)

# beta[12] - low neff
trace.beta.12 <- trace_plot(beta_draws[,12], chains = 1, iterations = 3400000, warmup = 340000, title= "Trace plot beta[12] - synthetic male lung cancer")
save.figure(trace.beta.12, name="trace_beta_12", path=storage_path, pdf=F)

trace.beta.12.first <- trace_plot(beta_draws[1:5000,12], chains = 1, iterations = 5000, warmup = 0, title = "Trace plot beta[12], first 5000 - synthetic male lung cancer")
save.figure(trace.beta.12.first, name = "trace_beta_12_first", path = storage_path, pdf = F)

trace.beta.12.last <- trace_plot(beta_draws[3055001:3060000,12], chains = 1, iterations = 5000, warmup = 0, title = "Trace plot beta[12], last 5000 - synthetic male lung cancer")
save.figure(trace.beta.12.last, name = "trace_beta_12_last", path = storage_path, pdf = F)

trace.beta.9 <- trace_plot(beta_draws[,9], chains = 1, iterations = 3400000, warmup = 340000, title= "Trace plot beta[9] - synthetic male lung cancer")
save.figure(trace.beta.9, name="trace_beta_9", path=storage_path, pdf=F)

trace.beta.18 <- trace_plot(beta_draws[,18], chains = 1, iterations = 3400000, warmup = 340000, title= "Trace plot beta[18] - synthetic male lung cancer")
save.figure(trace.beta.18, name="trace_beta_18", path=storage_path, pdf=F)

trace.kappa.1.first <- trace_plot(kappa_draws[1:5000,1], chains = 1, iterations = 5000, warmup = 0, title = "Trace plot kappa[1], first 5000 - synthetic male lung cancer")
save.figure(trace.kappa.1.first, name = "trace_kappa_1_first", path = storage_path, pdf = F)

trace.kappa.1.last <- trace_plot(kappa_draws[3055001:3060000,1], chains = 1, iterations = 5000, warmup = 0, title = "Trace plot kappa[1], last 5000 - synthetic male lung cancer")
save.figure(trace.kappa.1.last, name = "trace_kappa_1_last", path = storage_path, pdf = F)

trace.kappa.9.first <- trace_plot(kappa_draws[1:5000,9], chains = 1, iterations = 5000, warmup = 0, title = "Trace plot kappa[9], first 5000 - synthetic male lung cancer")
save.figure(trace.kappa.9.first, name = "trace_kappa_9_first", path = storage_path, pdf = F)

trace.kappa.9.last <- trace_plot(kappa_draws[3055001:3060000,9], chains = 1, iterations = 5000, warmup = 0, title = "Trace plot kappa[9], last 5000 - synthetic male lung cancer")
save.figure(trace.kappa.9.last, name = "trace_kappa_9_last", path = storage_path, pdf = F)


#   ----   Plot stan results   ----

stan.marginals <- list(intercept_draws = intercept_draws,
                  tau_epsilon_draws = tau_epsilon_draws,
                  tau_alpha_draws = tau_alpha_draws,
                  tau_beta_draws = tau_beta_draws,
                  tau_kappa_draws = tau_kappa_draws,
                  alpha_draws = alpha_draws,
                  beta_draws = beta_draws,
                  kappa_draws = kappa_draws,
                  eta_draws = eta_draws_reduced)

# intercept.marginal = data.frame(int = marginals$intercept_draws)
# library("bayesplot")
# test_int_dens <- mcmc_dens(intercept.marginal, pars=c("int"))
# test_int_dens
# 
# gg_test <- ggplot() + geom_histogram(data=intercept.marginal, aes(x = int)); gg_test

source("Scripts/Synthetic data/plot_stan_vs_underlying.R")


# when period effect is modelled as rw2, summed to zero


# stan.res <- produce.stan.plots(stan_df=stan_lc_df,
#                                underlying.effects=underlying.effects.lc,
#                                plot.func=plot.stan.vs.underlying.lc.rw2,
#                                save.func=save.stan.plots.lc.rw2,
#                                path.to.storage=storage_path,
#                                summaries.func=produce.summaries.stan.lc.rw2)

# synthetic cancer configs:
stan.res <- produce.stan.plots(stan_df=stan_lc_df,
                               underlying.effects=underlying.effects.lc,
                               plot.func=plot.stan.vs.underlying.synthetic.cancer,
                               save.func=save.stan.plots.lc.rw2,
                               path.to.storage=storage_path,
                               summaries.func=produce.summaries.stan.lc.rw2)


#   ----    Plot stan and inlabru-results together   ----
source("Scripts/Synthetic data/plot_inlabru_stan_compared.R")

# undrifted stan
# plots_compared <- produce.compared.plots(
#   stan.summaries = stan.res$summaries,
#   inlabru.summaries = plots.summaries.inlabru$summaries,
#   underlying.effects = underlying.effects.lc,
#   plot.func = plot.inlabru.stan.compared.lc,
#   save.func = save.compared.undrifted.lc,
#   path.to.storage=storage_path)

# drifted stan
underlying.effects.lc <- c(underlying.effects.lc, age.intercept.true = underlying.effects.lc$intercept)

plots_compared <- produce.compared.plots(
  stan.summaries = stan.res$summaries,
  stan.marginals = stan.marginals,
  inlabru.summaries = plots.summaries.inlabru$summaries,
  res.inlabru = res.inlabru.lc.1,
  underlying.effects = underlying.effects.lc,
  plot.func = function(...) {plot.inlabru.stan.compared.rw2(..., cohort=FALSE, tau.beta.cutoff = 300, tau.kappa.cutoff = 500, tau.alpha.cutoff = 100)},
  save.func = function(...) {save.compared.rw2(..., cohort=FALSE)},
  path.to.storage=storage_path)


#   ----    generate and plot posterior counts   ----

# inlabru

inlabru.samps <- generate(
  res.inlabru.lc.1,
  data = data.frame(x = obs.lc$x, t = obs.lc$t, x.c = obs.lc$x.c, xt = obs.lc$xt),
  formula = ~ Int + alpha + beta*kappa + epsilon,
  n.sample = 1000)

inlabru.lambda <- obs.lc$E * exp(inlabru.samps)

inlabru.Y <- matrix(rpois(324*1000, lambda = inlabru.lambda), nrow = 324, ncol = 1000)

inlabru.Y.df <- data.frame(x = obs.lc$x, t = obs.lc$t, xt = obs.lc$xt,
                           mean = apply(inlabru.Y, 1, mean),
                           X0.975 = apply(inlabru.Y, 1, quantile, 0.975),
                           X0.025 = apply(inlabru.Y, 1, quantile, 0.025))

# stan

stan.samps <- eta_draws_reduced[sample(nrow(eta_draws_reduced), size = 1000, replace = F),]
stan.samps <- t(stan.samps)  # transpose to get on same format as inlabru samples

stan.lambda <- obs.lc$E * exp(stan.samps)

stan.Y <- matrix(rpois(324*1000, lambda = stan.lambda), nrow = 324, ncol = 1000)
stan.Y.df <- data.frame(x = obs.lc$x, t = obs.lc$t, xt = obs.lc$xt,
                        mean = apply(stan.Y, 1, mean),
                        X0.975 = apply(stan.Y, 1, quantile, 0.975),
                        X0.025 = apply(stan.Y, 1, quantile, 0.025))

comparison.Y <- inlabru.Y.df %>%
  left_join(stan.Y.df, by = c("x" = "x", "t" = "t", "xt" = "xt"), suffix = c(".inlabru", ".stan")) %>%
  mutate(Y.observed = obs.lc$Y)

source("Scripts/Functions/plotters.R")

plot.counts.inlabru.stan.compared(comparison.Y, path.to.storage = storage_path, png = F)

#   ----   Compare marginals of predictor from stan and inlabru   ----

# inlabru

inlabru.samps.predictor <- generate(
  res.inlabru.lc.1,
  data = data.frame(x = obs.lc$x, t = obs.lc$t, x.c = obs.lc$x.c, xt = obs.lc$xt),
  formula = ~ Int + alpha + beta*kappa + epsilon,
  n.sample = 10000)

inlabru.predictor.df <- data.frame(t(inlabru.samps.predictor))

# stan

stan.samps.predictor <- eta_draws[sample(nrow(eta_draws), size = 10000, replace = F),]
#stan.samps <- eta_draws_reduced[sample(nrow(eta_draws_reduced), size = 1000, replace = F),]

stan.predictor.df <- data.frame(stan.samps.predictor)


comparison.Y <- inlabru.Y.df %>%
  left_join(stan.Y.df, by = c("x" = "x", "t" = "t", "xt" = "xt"), suffix = c(".inlabru", ".stan")) %>%
  mutate(Y.observed = obs.lc$Y)

source("Scripts/Functions/plotters.R")

plot.predictor.inlabru.stan.compared(inlabru.predictor.df, stan.predictor.df, path.to.storage = storage_path, a45=T)


  
