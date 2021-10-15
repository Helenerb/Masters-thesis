# TODO: Script running and comparing inference with the lc-cohort model
# by inlabru and STAN. 

setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data")
source("configurations_synthetic_data.R")

# cohfiguration with cohort effect
#underlying.effects.lc.cohort <- configuration.v7()  ##  Fine grid
underlying.effects.lc.cohort <- configuration.v17()  ## First attempt at coarser grid
#underlying.effects.lc.cohort <- configuration.v18()  ## more erratic beta
#underlying.effects.lc.cohort <- configuration.v19()
#underlying.effects.lc.cohort <- configuration.v20()
#underlying.effects.lc.cohort <- configuration.v21()
#underlying.effects.lc.cohort <- configuration.v22()
#underlying.effects.lc.cohort <- configuration.v23()

figures.folder = "/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Output/Figures"

#storage_path = file.path(figures.folder, "v7")
#storage_path = file.path(figures.folder, "v17")
storage_path = file.path(figures.folder, "v17d")
#storage_path = file.path(figures.folder, "v18")
#storage_path = file.path(figures.folder, "v18d")
#storage_path = file.path(figures.folder, "v19")
#storage_path = file.path(figures.folder, "v20")
#storage_path = file.path(figures.folder, "v21")
#storage_path = file.path(figures.folder, "v22")
#storage_path = file.path(figures.folder, "v23")

obs.cohort <- underlying.effects.lc.cohort$obs

source("Inlabru\ analyses/inlabru_analyses.R")
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

#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v17/stan_results/stan_v17.Rda")
#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v18/stan_results/stan_v18.Rda")
#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v18d/stan_results/stan_v18_drift.Rda")
load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v17d/stan_results/stan_v17d.Rda")

# stan results df loaded under name "stan_lc_df"

source("plot_stan_vs_underlying.R")

produce.stan.plots(stan_df=stan_lc_df,
                   underlying.effects=underlying.effects.lc.cohort,
                   plot.func=plot.stan.vs.underlying.cohort.undrifted,
                   save.func=save.stan.plots.cohort.undrifted,
                   path.to.storage=storage_path)



stan.plots <- plot.stan.vs.underlying.cohort(stan_df = stan_lc_df, obs=obs.cohort, underlying.effects = underlying.effects.lc.cohort)
save.stan.plots(plots=stan.plots, path.to.storage=storage_path)
stan.plots$p.alpha
stan.plots$p.beta
stan.plots$p.eta
stan.plots$p.gamma
stan.plots$p.kappa
