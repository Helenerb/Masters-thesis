#Script running and comparing inference with the lc-cohort model
# by inlabru and STAN. 

setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data")
source("configurations_synthetic_data.R")

# cohfiguration with cohort effect
#underlying.effects.lc.cohort <- configuration.v7()  ##  Fine grid
#underlying.effects.lc.cohort <- configuration.v17()  ## First attempt at coarser grid
#underlying.effects.lc.cohort <- configuration.v17.1()
#underlying.effects.lc.cohort <- configuration.v17.3()
#underlying.effects.lc.cohort <- configuration.v17.4()
#underlying.effects.lc.cohort <- configuration.v18()  ## more erratic beta
#underlying.effects.lc.cohort <- configuration.v18.1()
underlying.effects.lc.cohort <- configuration.v18.3()
#underlying.effects.lc.cohort <- configuration.v19()
#underlying.effects.lc.cohort <- configuration.v20()
#underlying.effects.lc.cohort <- configuration.v21()
#underlying.effects.lc.cohort <- configuration.v22()
#underlying.effects.lc.cohort <- configuration.v22.1()
#underlying.effects.lc.cohort <- configuration.v22.3()
#underlying.effects.lc.cohort <- configuration.v23()
#underlying.effects.lc.cohort <- configuration.v23.1()

figures.folder = "/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Output/Figures"

#storage_path = file.path(figures.folder, "v7")
#storage_path = file.path(figures.folder, "v17")
#storage_path = file.path(figures.folder, "v17d")
#storage_path = file.path(figures.folder, "v17dh")
#storage_path = file.path(figures.folder, "v17_only_kappa")
#storage_path = file.path(figures.folder, "v17_3_only_kappa")
#storage_path = file.path(figures.folder, "v17_3_ar1c")
#storage_path = file.path(figures.folder, "v17_3_ar1c_extraconstr_gamma")
#storage_path = file.path(figures.folder, "v17_4")
#storage_path = file.path(figures.folder, "v17_3_extraconstr_gamma")
#storage_path = file.path(figures.folder, "v18")
#storage_path = file.path(figures.folder, "v18d")
#storage_path = file.path(figures.folder, "v18dh")
storage_path = file.path(figures.folder, "v18_3")
#storage_path = file.path(figures.folder, "v18_3_extraconstr")
#storage_path = file.path(figures.folder, "v19")
#storage_path = file.path(figures.folder, "v20")
#storage_path = file.path(figures.folder, "v21")
#storage_path = file.path(figures.folder, "v22")
#storage_path = file.path(figures.folder, "v22_1")
#storage_path = file.path(figures.folder, "v22_3")
#storage_path = file.path(figures.folder, "v22_3_extraconstr_gamma")
#storage_path = file.path(figures.folder, "v23")
#storage_path = file.path(figures.folder, "v23_1")


obs.cohort <- underlying.effects.lc.cohort$obs

source("Inlabru\ analyses/inlabru_analyses.R")
#runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.undrifted.period.cohort.2(obs.cohort)})
#runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.undrifted.period.cohort(obs.cohort)})
#runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.undrifted.period.cohort.2.gamma.extraconstr(obs.cohort)})
#runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.lc.cohort.kappa_high_variance_prior(obs.cohort)})
#runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.undrifted.period.cohort.2.beta.rw(obs.cohort)})
#runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.ar1c.cohort.2(obs.cohort)})
#runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.ar1c.cohort.2.gamma.extraconstr(obs.cohort)})
#runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.rw2.cohort.2(obs.cohort)})
runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.rw2.cohort.2.gamma.extraconstr(obs.cohort)})


###   ----   Plot results from inlabru inference   ----  

source("plot_inlabru_vs_underlying.R")

# plot with ar1c configuration
# plots.summaries.inlabru <- plot.inlabru.vs.underlying.cohort.ar1c.2(
#   res.inlabru.lc.1,
#   underlying.effects.lc.cohort,
#   path.to.storage = storage_path,
#   save=TRUE,
#   phi.plus.kappa.func = phi.plus.kappa.v17)

#plot.period.posterior <- plot.posterior.period.effects(res.inlabru.lc.1, underlying.effects.lc.cohort)

# plotting results from run with cohort effects:

# only kappa rw as period effect, with kappa_0 = 0 constraints
# plots.summaries.inlabru <- plot.inlabru.vs.underlying.cohort.only.kappa(
#   res.inlabru.lc.1,
#   underlying.effects.lc.cohort,path.to.storage = storage_path,
#   save=TRUE,
#   phi.plus.kappa.func = phi.plus.kappa.v17)

# only kappa rw as period effect, with kappa summed to zero
plots.summaries.inlabru <- plot.inlabru.vs.underlying.cohort.only.kappa.2(
  res.inlabru.lc.1,
  underlying.effects.lc.cohort,
  path.to.storage = storage_path,
  save=TRUE,
  phi.plus.kappa.func = phi.plus.kappa.v17)


plots.summaries.inlabru$plots$p.alpha
plots.summaries.inlabru$plots$p.beta
plots.summaries.inlabru$plots$p.phi
plots.summaries.inlabru$plots$p.intercept
plots.summaries.inlabru$plots$p.kappa
plots.summaries.inlabru$plots$p.eta
plots.summaries.inlabru$plots$p.eta.2
plots.summaries.inlabru$plots$p.eta.t
plots.summaries.inlabru$plots$p.eta.x
plots.summaries.inlabru$plots$p.gamma

###    ----   Configure and run inference with STAN    ----

#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v17/stan_results/stan_v17.Rda")
#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v18/stan_results/stan_v18.Rda")
#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v18d/stan_results/stan_v18d.Rda")
#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v17d/stan_results/stan_v17d.Rda")
#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v17dh/stan_results/stan_v17dh.Rda")
#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v18dh/stan_results/stan_v18dh.Rda")

load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v18_3/stan_results/stan_v18_3.Rda")
load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v18_3/stan_results/draws_tau_kappa.RData")

stan_fit_full <- readRDS("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v18_3/stan_results/stan_fit.rds")


# stan results df loaded under name "stan_lc_df"

source("plot_stan_vs_underlying.R")

# undrifted version of stan program:
# produce.stan.plots(stan_df=stan_lc_df,
#                    underlying.effects=underlying.effects.lc.cohort,
#                    plot.func=plot.stan.vs.underlying.cohort.undrifted,
#                    save.func=save.stan.plots.cohort.undrifted,
#                    path.to.storage=storage_path)

# drifted version, soft constraints
# stan.res <- produce.stan.plots(stan_df=stan_lc_df,
#                                underlying.effects=underlying.effects.lc.cohort,
#                                plot.func=plot.stan.vs.underlying.cohort.drifted,
#                                save.func=save.stan.plots.cohort.undrifted,
#                                path.to.storage=storage_path)

# drifted version of stan program, hard constraints
stan.res <- produce.stan.plots(stan_df=stan_lc_df,
                   underlying.effects=underlying.effects.lc.cohort,
                   plot.func=plot.stan.vs.underlying.cohort.drifted,
                   save.func=save.stan.plots.cohort.undrifted,
                   path.to.storage=storage_path,
                   summaries.func=produce.summaries.stan.hard)


stan.plots <- plot.stan.vs.underlying.cohort(stan_df = stan_lc_df, obs=obs.cohort, underlying.effects = underlying.effects.lc.cohort)
save.stan.plots(plots=stan.plots, path.to.storage=storage_path)
stan.plots$p.alpha
stan.plots$p.beta
stan.plots$p.eta
stan.plots$p.gamma
stan.plots$p.kappa


#   ----   Plot comparison   ----   
source("plot_inlabru_stan_compared.R")

# undrifted stan
# plots_compared <- produce.compared.plots(
#   stan.summaries = stan.res$summaries,
#   inlabru.summaries = plots.summaries.inlabru$summaries,
#   underlying.effects = underlying.effects.lc.cohort,
#   plot.func = plot.inlabru.stan.compared.cohort,
#   save.func = save.compared.undrifted.cohort,
#   path.to.storage=storage_path)

# drifted stan
plots_compared <- produce.compared.plots(
  stan.summaries = stan.res$summaries,
  inlabru.summaries = plots.summaries.inlabru$summaries,
  underlying.effects = underlying.effects.lc.cohort,
  plot.func = plot.inlabru.stan.compared.cohort,
  save.func = save.compared.drifted.cohort,
  path.to.storage=storage_path)
