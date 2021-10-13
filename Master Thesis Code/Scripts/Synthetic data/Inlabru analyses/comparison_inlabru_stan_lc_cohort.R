# TODO: Script running and comparing inference with the lc-cohort model
# by inlabru and STAN. 

setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data")
source("configurations_synthetic_data.R")

# cohfiguration with cohort effect
#underlying.effects.lc.cohort <- configuration.v7()  ##  Fine grid
#underlying.effects.lc.cohort <- configuration.v17()  ## First attempt at coarser grid
#underlying.effects.lc.cohort <- configuration.v18()  ## more erratic beta
#underlying.effects.lc.cohort <- configuration.v19()
#underlying.effects.lc.cohort <- configuration.v20()
#underlying.effects.lc.cohort <- configuration.v21()
underlying.effects.lc.cohort <- configuration.v22()
#underlying.effects.lc.cohort <- configuration.v23()

figures.folder = "/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Output/Figures"

#storage_path = file.path(figures.folder, "v7")
#storage_path = file.path(figures.folder, "v17")
#storage_path = file.path(figures.folder, "v18")
#storage_path = file.path(figures.folder, "v19")
#storage_path = file.path(figures.folder, "v20")
#storage_path = file.path(figures.folder, "v21")
storage_path = file.path(figures.folder, "v22")
#storage_path = file.path(figures.folder, "v23")

obs.lc.cohort <- underlying.effects.lc.cohort$obs

source("inlabru_analyses.R")
runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.lc.cohort.kappa_high_variance_prior(obs.lc.cohort)})


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
