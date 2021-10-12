# TODO: Script running and comparing inference with the lc-cohort model
# by inlabru and STAN. 

setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data")
source("configurations_synthetic_data.R")

# cohfiguration with cohort effect
underlying.effects.lc.cohort <- configuration.v7()  ##  TODO: change to config with coarser grid?

storage_path = file.path("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Output/Figures", "v7")

obs.lc.cohort <- underlying.effects.lc.cohort$obs

source("inlabru_analyses.R")
runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.lc.cohort.kappa_high_variance_prior(obs.lc.cohort)})


###   ----   Plot results from inlabru inference   ----  

source("plot_inlabru_vs_underlying.R")

#plot.period.posterior <- plot.posterior.period.effects(res.inlabru.lc.1, underlying.effects.lc.cohort)

# plotting results from run with cohort effects:
plots.lc.cohort <- plot.inlabru.vs.underlying.cohort(
  res.inlabru.lc.1, 
  underlying.effects.lc.cohort,path.to.storage = "/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Output/Figures/v7",
  save=TRUE)
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
