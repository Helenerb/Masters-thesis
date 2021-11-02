library(patchwork)

setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data")
source("configurations_synthetic_data.R")

# configuration without cohort effect
#underlying.effects.lc <- configuration.v5()
#underlying.effects.lc <- configuration.v9()  #  config with coarser grid
#underlying.effects.lc <- configuration.v10()
#underlying.effects.lc <- configuration.v10.1()# config with coarser grid and larger variance in beta
#underlying.effects.lc <- configuration.v10.2()
underlying.effects.lc <- configuration.v10.3()
#underlying.effects.lc <- configuration.v11()  # config with coarser grid, larger variance in beta and steeper phi (compared to alpha)
#underlying.effects.lc <- configuration.v11.1()
#underlying.effects.lc <- configuration.v11.3()
#underlying.effects.lc <- configuration.v12()  #  config with coarser grid, smaller variance in beta, steeper phi (compared to alpha)
#underlying.effects.lc <- configuration.v12.3()
#underlying.effects.lc <- configuration.v13()  # config with coarser grid, even smaller variance in beta, a bit less steep phi, higher variance in kappa
#underlying.effects.lc <- configuration.v14()
#underlying.effects.lc <- configuration.v15()  # half the grid of the original v5

#underlying.effects.lc <- configuration.test.1()

figures.folder = "/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Output/Figures"

#storage_path = file.path(figures.folder, "v10_2")
#storage_path = file.path(figures.folder, "v10_2_only_kappa")
#storage_path = file.path(figures.folder, "v10_3_only_kappa")
#storage_path = file.path(figures.folder, "v10d")
#storage_path = file.path(figures.folder, "v10dh")
#storage_path = file.path(figures.folder, "v10_2_ar1c")
#storage_path = file.path(figures.folder, "v10_3_ar1c")
storage_path = file.path(figures.folder, "v10_3_rw2")
#storage_path = file.path(figures.folder, "v11_3_rw2")
#storage_path = file.path(figures.folder, "v12_3_rw2")

obs.lc <- underlying.effects.lc$obs

source("Inlabru\ analyses/inlabru_analyses.R")
#runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.lc.kappa_high_variance_prior(obs.lc)})
#runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.undrifted.period.2(obs.lc)})
#runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.ar1c.lc(obs.lc)})
#runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.ar1c.lc.2(obs.lc)})
runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.rw2.lc.2(obs.lc)})

source("plot_inlabru_vs_underlying.R")

# plots.summaries.inlabru <- plot.inlabru.vs.underlying.lc.ar1c.2(
#   res.inlabru.lc.1,
#   underlying.effects.lc,
#   path.to.storage = storage_path,
#   save=TRUE,
#   phi.plus.kappa.func = phi.plus.kappa.v17)

# plots.summaries.inlabru <- plot.inlabru.vs.underlying.lc.ar1c(
#   res.inlabru.lc.1, 
#   underlying.effects.lc,
#   path.to.storage = storage_path,
#   save=TRUE,
#   phi.plus.kappa.func = phi.plus.kappa.v17)

plots.summaries.inlabru <- plot.inlabru.vs.underlying.lc.only.kappa.2(
  res.inlabru.lc.1,
  underlying.effects.lc,
  path.to.storage = storage_path,
  save=TRUE,
  phi.plus.kappa.func = phi.plus.kappa.v17)


print("Runtime for inlabru: ")
print(runtime.inlabru)

#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v10/stan_results/stan_v10.Rda")
#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v10d/stan_results/stan_v10d.Rda")
#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v10dh/stan_results/stan_v10dh.Rda")
load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v10_3/stan_results/stan_v10_3.Rda")
#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v11_3/stan_results/stan_v11_3.Rda")
#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v12_3/stan_results/stan_v12_3.Rda")

#   ----   load STAN marginals   ---- 

path.to.stan.results = "/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master\ Thesis\ Code/Scripts/Synthetic\ data/Stan analyses/v10_3/stan_results"

load(file=file.path(path.to.stan.results, "draws_intercept.RData"))
load(file=file.path(path.to.stan.results, "draws_tau_epsilon.RData"))
load(file.path(path.to.stan.results, "draws_tau_alpha.RData"))
load(file.path(path.to.stan.results, "draws_tau_beta.RData"))
load(file.path(path.to.stan.results, "draws_tau_kappa.RData"))
load(file.path(path.to.stan.results, "draws_alpha.RData"))
load(file.path(path.to.stan.results, "draws_beta.RData"))
load(file.path(path.to.stan.results, "draws_kappa.RData"))
load(file.path(path.to.stan.results, "draws_eta.RData"))

stan.marginals <- list(intercept_draws = intercept_draws,
                  tau_epsilon_draws = tau_epsilon_draws,
                  tau_alpha_draws = tau_alpha_draws,
                  tau_beta_draws = tau_beta_draws,
                  tau_kappa_draws = tau_kappa_draws,
                  alpha_draws = alpha_draws,
                  beta_draws = beta_draws,
                  kappa_draws = kappa_draws,
                  eta_draws = eta_draws)

# intercept.marginal = data.frame(int = marginals$intercept_draws)
# library("bayesplot")
# test_int_dens <- mcmc_dens(intercept.marginal, pars=c("int"))
# test_int_dens
# 
# gg_test <- ggplot() + geom_histogram(data=intercept.marginal, aes(x = int)); gg_test

source("plot_stan_vs_underlying.R")

#undrifted version of stan program:
# stan.results <- produce.stan.plots(stan_df=stan_lc_df,
#                    underlying.effects=underlying.effects.lc,
#                    plot.func=plot.stan.vs.underlying.lc.undrifted,
#                    save.func=save.stan.plots.lc,
#                    path.to.storage=storage_path,
#                    summaries.func=produce.summaries.stan.lc)

# drifted version, soft constraints
# stan.res <- produce.stan.plots(stan_df=stan_lc_df,
#                                underlying.effects=underlying.effects.lc,
#                                plot.func=plot.stan.vs.underlying.lc.drifted,
#                                save.func=save.stan.plots.lc,
#                                path.to.storage=storage_path,
#                                summaries.func=produce.summaries.stan.lc)

# drifted version of stan program, hard constraints
# stan.res <- produce.stan.plots(stan_df=stan_lc_df,
#                                underlying.effects=underlying.effects.lc,
#                                plot.func=plot.stan.vs.underlying.lc.drifted,
#                                save.func=save.stan.plots.lc,
#                                path.to.storage=storage_path,
#                                summaries.func=produce.summaries.stan.hard.lc)

# when period effect is modelled as rw2, summed to zero
stan.res <- produce.stan.plots(stan_df=stan_lc_df,
                               underlying.effects=underlying.effects.lc,
                               plot.func=plot.stan.vs.underlying.lc.rw2,
                               save.func=save.stan.plots.lc.rw2,
                               path.to.storage=storage_path,
                               summaries.func=produce.summaries.stan.lc.rw2)


#   ----    Plot stan and inlabru-results together   ----
source("plot_inlabru_stan_compared.R")

# undrifted stan
# plots_compared <- produce.compared.plots(
#   stan.summaries = stan.res$summaries,
#   inlabru.summaries = plots.summaries.inlabru$summaries,
#   underlying.effects = underlying.effects.lc,
#   plot.func = plot.inlabru.stan.compared.lc,
#   save.func = save.compared.undrifted.lc,
#   path.to.storage=storage_path)

# drifted stan
plots_compared <- produce.compared.plots(
  stan.summaries = stan.res$summaries,
  stan.marginals = stan.marginals,
  inlabru.summaries = plots.summaries.inlabru$summaries,
  res.inlabru = res.inalbru.lc.1,
  underlying.effects = underlying.effects.lc,
  plot.func = function(...) {plot.inlabru.stan.compared.rw2(..., cohort=FALSE)},
  save.func = function(...) {save.compared.rw2(..., cohort=FALSE)},
  path.to.storage=storage_path)
