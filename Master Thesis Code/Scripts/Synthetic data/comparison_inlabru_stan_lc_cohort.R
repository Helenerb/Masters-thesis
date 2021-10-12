# TODO: Script running and comparing inference with the lc-cohort model
# by inlabru and STAN. 

plot.posterior.period.effects <- function(res.inlabru, underlying.effects){
  
  nt <- underlying.effects$nt
  samps = inla.posterior.sample(res.inlabru.lc.1, n = 1000)
  
  phi.plus.kappa <- function(){
    t = 0:(nt-1)
    res = kappa + phi*t
    return(res)
  }
  
  posterior.phi.kappa <- inla.posterior.sample.eval(fun = phi.plus.kappa, samples=samps)
  
  mean.post = apply(posterior.phi.kappa, 1, mean)
  q1.post <- apply(posterior.phi.kappa, 1, quantile, 0.025)
  q2.post <- apply(posterior.phi.kappa, 1, quantile, 0.975)
  
  posterior.phi.kappa.df <- data.frame(t = 1:nt,
                                       mean = apply(posterior.phi.kappa, 1, mean),
                                       q1 = apply(posterior.phi.kappa, 1, quantile, 0.025),
                                       q2 = apply(posterior.phi.kappa, 1, quantile, 0.975)) %>%
    mutate(kappa = underlying.effects$kappa.true[t]) %>%
    mutate(phi.t = underlying.effects$phi.true*(t-1)) %>%
    mutate(kappa.phi = kappa + phi.t)
  
  gg.posterior <- ggplot(data = posterior.phi.kappa.df) +
    geom_ribbon(aes(x = t, ymin = q1, ymax = q2, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_point(aes(x = t, y = mean, color = "Inlabru", fill = "Inlabru")) +
    geom_point(aes(x = t, y = kappa.phi, color = "True", fill = "True")) +
    scale_color_manual(name = " ", values = palette.basis) + 
    scale_fill_manual(name = " ", values = palette.basis) + 
    labs(title = "Phi*t + kappa", x = "t", y = "")
  
  return(gg.posterior)
}

setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data")
source("configurations_synthetic_data.R")

# cohfiguration with cohort effect
underlying.effects.lc.cohort <- configuration.v7()  ##  TODO: change to config with coarser grid?

obs.lc.cohort <- underlying.effects.lc.cohort$obs

source("inlabru_analyses.R")
runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.lc.cohort.kappa_high_variance_prior(obs.lc.cohort)})

###   ----   Posterior period effect   ----   

nt <- underlying.effects.lc.cohort$nt
samps = inla.posterior.sample(res.inlabru.lc.1, n = 1000)

phi.plus.kappa <- function(){
  t = 0:(nt-1)
  res = kappa + phi*t
  return(res)
}

posterior.phi.kappa <- inla.posterior.sample.eval(fun = phi.plus.kappa, samples=samps)

mean.post = apply(posterior.phi.kappa, 1, mean)
q1.post <- apply(posterior.phi.kappa, 1, quantile, 0.025)
q2.post <- apply(posterior.phi.kappa, 1, quantile, 0.975)

posterior.phi.kappa.df <- data.frame(t = 1:nt,
                                     mean = apply(posterior.phi.kappa, 1, mean),
                                     q1 = apply(posterior.phi.kappa, 1, quantile, 0.025),
                                     q2 = apply(posterior.phi.kappa, 1, quantile, 0.975)) %>%
  mutate(kappa = underlying.effects.lc.cohort$kappa.true[t]) %>%
  mutate(phi.t = underlying.effects.lc.cohort$phi.true*(t-1)) %>%
  mutate(kappa.phi = kappa + phi.t)

gg.posterior <- ggplot(data = posterior.phi.kappa.df) +
  geom_ribbon(aes(x = t, ymin = q1, ymax = q2, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_point(aes(x = t, y = mean, color = "Inlabru", fill = "Inlabru")) +
  geom_point(aes(x = t, y = kappa.phi, color = "True", fill = "True")) +
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) + 
  labs(title = "Phi*t + kappa", x = "t", y = "")

gg.posterior


###   ----   Plot results from inlabru inference   ----  

source("plot_inlabru_vs_underlying.R")

plot.period.posterior <- plot.posterior.period.effects(res.inlabru.lc.1, underlying.effects.lc.cohort)

# plotting results from run with cohort effects:
plots.lc.cohort <- plot.inlabru.vs.underlying.v7(res.inlabru.lc.1,
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
