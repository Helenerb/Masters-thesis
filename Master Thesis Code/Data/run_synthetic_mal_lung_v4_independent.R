# independent script for producing and running synthetic male lung cancer - v4

library(inlabru)
library(tidyverse)

synthetic.male.lung.v4 <- function(){
  obs <- read.csv("synthetic_male_lung_4.csv")
  
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

inlabru.rw2.lc.2 <- function(obs, max_iter=30){
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
  loggamma.prior <- list(prec = list(prior = 'loggamma', param = c(1,0.00005), initial = 4))
  loggamma.prior.high.variance <- list(prec = list(prior = 'loggamma', param = c(1,0.005), initial = 4))
  
  comp = ~ -1 +
    Int(1) +
    alpha(x, model = "rw1", values=unique(obs$x), hyper = loggamma.prior, constr = TRUE) +
    beta(x.c, model = "iid", extraconstr = list(A = A.beta, e = e.beta), hyper = loggamma.prior) +
    kappa(t, model = "rw2", values = unique(obs$t), constr = TRUE, hyper = loggamma.prior.high.variance) +
    epsilon(xt, model = "iid", hyper = loggamma.prior)
  
  formula = Y ~ Int + alpha + beta*kappa + epsilon
  
  likelihood = like(formula = formula, family = "poisson", data = obs, E = obs$E)
  
  c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)  # control.compute
  
  res.inlabru = bru(components = comp,
                    likelihood, 
                    options = list(verbose = F,
                                   bru_verbose = 1, 
                                   num.threads = "1:1",
                                   control.compute = c.c,
                                   bru_max_iter=max_iter,
                                   control.predictor = list(link = 1)
                    ))
  return(res.inlabru)
}

config_data <- synthetic.male.lung.v4()

underlying.effects.lc <- config_data$underlying.effects

obs.lc <- underlying.effects.lc$obs

runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.rw2.lc.2(obs.lc)})

