# TODO: Script running and comparing inference with the lc-cohort model
# by inlabru and STAN. 

setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data")
source("configurations_synthetic_data.R")

# cohfiguration with cohort effect
underlying.effects.lc.cohort <- configuration.v7()  ##  TODO: change to config with coarser grid?

obs.lc.cohort <- underlying.effects.lc.cohort$obs

source("inlabru_analyses.R")
res.inlabru.lc.cohort.1 <- inlabru.lc.cohort.1(obs.lc.cohort)

###   ----   Plot results from inlabru inference   ----  

source("plot_inlabru_vs_underlying.R")

# plotting results from run with cohort effects:
plots.lc.cohort <- plot.inlabru.vs.underlying.v1(res.inlabru.lc.cohort.1,
                                                 underlying.effects.lc.cohort)
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
