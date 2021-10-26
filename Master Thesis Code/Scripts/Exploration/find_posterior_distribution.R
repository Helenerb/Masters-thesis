# full script for finding posterior period distribution

# define underlying data: 

library(INLA)
library(inlabru)
library(tidyverse)
library(ggplot2)

palette <- c('#70A4D4', '#ECC64B', '#93AD80', '#da9124', '#696B8D',
             '#3290c1', '#5d8060', '#D7B36A', '#826133', '#A85150')


configuration.v10.2 <- function(){
  # version of v9 with more erratic beta
  
  #   ----   Setiting seed for reproductiveness   ----
  
  seed = 324
  set.seed(seed)
  
  #   ----   Defining data structure   ----
  
  nx = 16  # number of age groups
  nt = 20  # number of time steps
  
  #   ----   Definining underlying effects   ----
  
  alpha.true = 3.9*cos(((1:nx*5 + 30)* pi)/50)
  alpha.true = alpha.true - mean(alpha.true)
  
  age.intercept.true = - 4
  
  phi.true = -1 # used to be -1
  
  tau.beta.true = 1000
  
  beta.true = rnorm(nx, mean=0, sd=sqrt(1/tau.beta.true))
  beta.true = beta.true - mean(beta.true) + 1/nx
  
  tau.kappa.true = 0.5
  kappa.true.increments = rnorm(nt-1, mean=0, sd=sqrt(1/tau.kappa.true))
  kappa.true.increments = kappa.true.increments - mean(kappa.true.increments)
  kappa.true = rep(0, nt)
  for (idx in 2:nt){
    kappa.true[idx] = kappa.true[idx - 1] + kappa.true.increments[idx-1]
  }
  kappa.true[2:nt] = kappa.true[2:nt] - mean(kappa.true[2:nt])
  
  tau.epsilon.true = 1000
  
  at.risk = 10**6/nx
  
  obs <- data.frame(expand.grid(x=1:nx, t=1:nt)) %>%
    mutate(alpha = alpha.true[x]) %>%
    mutate(age.intercept = age.intercept.true) %>%
    mutate(beta = beta.true[x]) %>%
    mutate(kappa = kappa.true[t]) %>%
    mutate(phi.t = phi.true*(t-1)) %>%
    mutate(epsilon = rnorm(length(x), mean = 0, sd = sqrt(1/tau.epsilon.true))) %>%
    mutate(E = at.risk) %>%
    mutate(eta = age.intercept + alpha + beta*phi.t + beta*kappa + epsilon) %>%
    mutate(Y = rpois(length(x), E*exp(eta))) %>%
    mutate(xt = seq_along(x)) %>%
    mutate(x = x - 1, t = t - 1, xt = xt - 1) %>%
    mutate(x.c = x, t.c = t)
  
  underlying.effects <- list(
    obs = obs,
    alpha.true = alpha.true[unique(obs$x) + 1],
    age.intercept.true= age.intercept.true,
    beta.true = beta.true[unique(obs$x) + 1],
    kappa.true = kappa.true[unique(obs$t) + 1],
    phi.true = phi.true,
    at.risk = at.risk,
    tau.beta.true = tau.beta.true,
    tau.kappa.true = tau.kappa.true,
    tau.epsilon.true= tau.epsilon.true,
    config_name = "v10.2",
    nx=nx,
    nt=nt
  )
  return(underlying.effects)
}

inlabru.lc.kappa_high_variance_prior <- function(obs){
  # set-up of inlabru analysis of Poisson Lee-Carter cohort model
  # separate intercept
  # alpha_x modeled as rw1, summed to zero. 
  
  # returns: bru-object containing results of analysis
  
  nx = length(unique(obs$x))
  nt = length(unique(obs$t))
  nc = length(unique(obs$c))
  
  #   ----   Inlabru model fitting   ----
  
  A.mat = matrix(1, nrow = 1, ncol = nx)  #  helper values for constraining of beta
  e.vec = 1  #  helper values for constraining of beta
  
  pc.prior <- list(prec = list(prior = "pc.prec", param = c(1,0.05)))
  loggamma.prior <- list(prec = list(prior = 'loggamma', param = c(1,0.00005), initial = 4))  # Default values. Sufficiently uninformative?
  loggamma.prior.high.variance <- list(prec = list(prior = 'loggamma', param = c(1,0.005), initial = 4))
  
  comp = ~ -1 +
    Int(1) +
    alpha(x, model = "rw1", values=unique(obs$x), hyper = loggamma.prior, constr = TRUE) +
    phi(t, model = "linear", prec.linear = 1, mean.linear = 0) +
    beta(x.c, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = loggamma.prior) +
    kappa(t.c, model = "rw1", values = unique(obs$t), constr = TRUE, hyper = loggamma.prior.high.variance) +
    epsilon(xt, model = "iid", hyper = loggamma.prior)
  
  formula = Y ~ Int + alpha + beta*phi + beta*kappa + epsilon
  
  likelihood = like(formula = formula, family = "poisson", data = obs, E = obs$E)
  
  c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)  # control.compute
  
  res.inlabru = bru(components = comp,
                    likelihood, 
                    options = list(verbose = F,
                                   bru_verbose = 1, 
                                   num.threads = "1:1",
                                   control.compute = c.c,
                                   bru_max_iter=30,
                                   control.predictor = list(link = 1)
                    ))
  return(res.inlabru)
}

underlying.effects <- configuration.v10.2()

obs <- underlying.effects$obs

res.inlabru <- inlabru.lc.kappa_high_variance_prior(obs)

phi.plus.kappa <- function(){
  t = 0:(20-1)
  res = kappa + phi*t
  return(res)
}

plot.posterior.period.effects <- function(res.inlabru, underlying.effects, phi.plus.kappa.func){
  
  nt = underlying.effects$nt
  
  samps = inla.posterior.sample(res.inlabru, n = 1000)
  
  posterior.phi.kappa <- inla.posterior.sample.eval(fun = phi.plus.kappa.func, samples=samps)
  
  posterior.phi.kappa.df <- data.frame(t = 1:nt,
                                       mean = apply(posterior.phi.kappa, 1, mean),
                                       q1 = apply(posterior.phi.kappa, 1, quantile, 0.025),
                                       q2 = apply(posterior.phi.kappa, 1, quantile, 0.975)) %>%
    mutate(kappa = underlying.effects$kappa.true[t]) %>%
    mutate(phi.t = underlying.effects$phi.true*(t-1)) %>%
    mutate(kappa.phi = kappa + phi.t)
  
  gg.posterior <- ggplot(data = posterior.phi.kappa.df) +
    geom_ribbon(aes(x = t, ymin = q1, ymax = q2, fill = "Estimated"), alpha = 0.5) + 
    geom_point(aes(x = t, y = mean, color = "Estimated", fill = "Estimated"), size=0.5) +
    geom_point(aes(x = t, y = kappa.phi, color = "True", fill = "True"), size=0.5) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) + 
    labs(title = "Phi*t + kappa", x = "t", y = "")
  
  return(list(posterior.plot = gg.posterior, posterior.data = posterior.phi.kappa.df))
}

p.posterior.period <- plot.posterior.period.effects(res.inlabru, underlying.effects, phi.plus.kappa)

plot.posterior <- p.posterior.period$posterior.plot
plot.posterior