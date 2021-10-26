# run inlabru with cohort effects

library(INLA)
library(inlabru)
library(tidyverse)
library(ggplot2)

palette <- c('#70A4D4', '#ECC64B', '#93AD80', '#da9124', '#696B8D',
             '#3290c1', '#5d8060', '#D7B36A', '#826133', '#A85150')

#   ----   period effect centered at zero

configuration.v17.3 <- function(){
  # version of v17 with first element of kappa set exactly to zero. 
  # around zero with separate intercept
  
  #   ----   Setiting seed for reproductiveness   ----
  
  seed = 324
  set.seed(seed)
  
  #   ----   Defining data structure   ----
  
  nx = 16   # number of age groups - recall zero-indexing!
  nt = 20  # number of time steps
  
  cohorts <- (-15):19  # cohort indices, c = t-x
  nc = length(cohorts)  # number of cohorts
  
  #   ----   Definining underlying effects   ----
  
  alpha.true = 3.9*cos(((1:nx*5 + 30)* pi)/50)
  alpha.true = alpha.true - mean(alpha.true)  # center around zero
  
  age.intercept.true = -4
  
  phi.true = -1
  
  tau.beta.true = 30000
  
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
  
  kappa.true.drifted = kappa.true + phi.true*(0:(nt-1))
  kappa.true.drifted = kappa.true.drifted - mean(kappa.true.drifted)
  
  tau.gamma.true = 250
  gamma.true = rnorm(nc, mean = 0, sd=sqrt(1/tau.gamma.true))
  gamma.true = gamma.true - mean(gamma.true)
  for (idx in 2:nc) {
    gamma.true[idx] = gamma.true[idx - 1] + gamma.true[idx]
  }
  # increase risk for cohorts 40-50:
  gamma.true[8:10] = gamma.true[8:10] + rnorm(3, mean=0.4, sd = sqrt(1/(tau.gamma.true)))
  
  gamma.true = gamma.true - mean(gamma.true)  # center around zero
  
  tau.epsilon.true = 1000
  
  at.risk = 10**6/nx
  
  obs <- data.frame(expand.grid(x=1:nx, t=1:nt)) %>%
    mutate(alpha = alpha.true[x]) %>%
    mutate(age.intercept = age.intercept.true) %>%
    mutate(beta = beta.true[x]) %>%
    mutate(kappa = kappa.true[t]) %>%
    mutate(phi.t = phi.true*(t-1)) %>%
    mutate(kappa.drifted = kappa.true.drifted[t]) %>%
    mutate(xt = seq_along(x)) %>%
    mutate(x = x - 1, t = t - 1, xt = xt - 1) %>%
    mutate(c = t-x) %>%
    mutate(gamma = gamma.true[c + 15 + 1]) %>%
    mutate(epsilon = rnorm(length(x), mean = 0, sd = sqrt(1/tau.epsilon.true))) %>%
    mutate(E = at.risk) %>%
    mutate(eta = age.intercept + alpha + beta*kappa.drifted + gamma + epsilon) %>%
    mutate(Y = rpois(length(x), E*exp(eta))) %>%
    mutate(x.c = x, t.c = t)
  
  underlying.effects <- list(
    obs = obs,
    alpha.true = alpha.true,
    beta.true = beta.true,
    kappa.true = kappa.true,
    kappa.drifted = kappa.true.drifted,
    phi.true = phi.true,
    gamma.true = gamma.true,
    at.risk = at.risk,
    age.intercept.true = age.intercept.true,
    nt = nt,
    nx = nx,
    config_name = "v17"
  )
  return(underlying.effects)
}

# function for running inlabru without extra constraints on gamma
inlabru.undrifted.period.cohort.2 <- function(obs){
  # set up of inlabru where the period effect is modelled only as a random walk (without drift)
  
  nx = length(unique(obs$x))
  nt = length(unique(obs$t))
  nc = length(unique(obs$c))
  
  #   ----   Inlabru model fitting   ----
  
  A.mat = matrix(1, nrow = 1, ncol = nx)  #  helper values for constraining of beta
  e.vec = 1  #  helper values for constraining of beta
  
  # constraint of period effect - first element set to zero
  # A.kappa = matrix(0, nrow = 1, ncol = nt)  
  # A.kappa[,1] = 1
  # print(A.kappa)
  # e.kappa = 0  #  helper values for constraining of beta
  
  pc.prior <- list(prec = list(prior = "pc.prec", param = c(1,0.05)))
  loggamma.prior <- list(prec = list(prior = 'loggamma', param = c(1,0.00005), initial = 4))  # Default values. Sufficiently uninformative?
  loggamma.prior.high.variance <- list(prec = list(prior = 'loggamma', param = c(1,0.005), initial = 4))
  
  comp = ~ -1 +
    Int(1) +
    alpha(x, model = "rw1", values=unique(obs$x), hyper = loggamma.prior, constr = TRUE) +
    #phi(t, model = "linear", prec.linear = 1, mean.linear = 0) +
    beta(x.c, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = loggamma.prior) +
    kappa(t.c, model = "rw1", values = unique(obs$t), constr=TRUE, hyper = loggamma.prior.high.variance) +
    gamma(c, model = "rw1", values = unique(obs$c), constr = TRUE, hyper = loggamma.prior) +
    epsilon(xt, model = "iid", hyper = loggamma.prior)
  
  formula = Y ~ Int + alpha + beta*kappa + gamma + epsilon
  
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

# function for running inlabru with extra constraints on gamma
inlabru.undrifted.period.cohort.2.gamma.extraconstr <- function(obs){
  # set up of inlabru where the period effect is modelled only as a random walk (without drift)
  
  nx = length(unique(obs$x))
  nt = length(unique(obs$t))
  nc = length(unique(obs$c))
  
  cohorts = -(nx-1) : (nt-1)
  cohorts.2 = (0:nt) - (0:nx)
  print(cohorts)
  print(cohorts.2)
  
  #   ----   Inlabru model fitting   ----
  
  A.mat = matrix(1, nrow = 1, ncol = nx)  #  helper values for constraining of beta
  e.vec = 1  #  helper values for constraining of beta
  
  # constraint of period effect - first element set to zero
  # A.kappa = matrix(0, nrow = 1, ncol = nt)  
  # A.kappa[,1] = 1
  # print(A.kappa)
  # e.kappa = 0  #  helper values for constraining of beta
  
  # extraconstraining of gamma:
  A.gamma = matrix(1, nrow = 1, ncol = nc)
  e.gamma = 0
  A.gamma[1,] = cohorts - mean(cohorts)
  print(A.gamma)
  
  
  pc.prior <- list(prec = list(prior = "pc.prec", param = c(1,0.05)))
  loggamma.prior <- list(prec = list(prior = 'loggamma', param = c(1,0.00005), initial = 4))  # Default values. Sufficiently uninformative?
  loggamma.prior.high.variance <- list(prec = list(prior = 'loggamma', param = c(1,0.005), initial = 4))
  
  comp = ~ -1 +
    Int(1) +
    alpha(x, model = "rw1", values=unique(obs$x), hyper = loggamma.prior, constr = TRUE) +
    #phi(t, model = "linear", prec.linear = 1, mean.linear = 0) +
    beta(x.c, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = loggamma.prior) +
    kappa(t.c, model = "rw1", values = unique(obs$t), constr=TRUE, hyper = loggamma.prior.high.variance) +
    gamma(c, model = "rw1", values = unique(obs$c), constr = TRUE, extraconstr = list(A = A.gamma, e = e.gamma), hyper = loggamma.prior) +
    epsilon(xt, model = "iid", hyper = loggamma.prior)
  
  formula = Y ~ Int + alpha + beta*kappa + gamma + epsilon
  
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

plot.res <- function(res.inlabru, underlying.effects, path.to.storage="", save=FALSE){
  
  obs <- underlying.effects$obs
  nx <- underlying.effects$nx
  nt <- underlying.effects$nt
  nc <- underlying.effects$nc
  
  data.alpha = cbind(res.inlabru$summary.random$alpha,
                     alpha.true = underlying.effects$alpha.true[res.inlabru$summary.random$alpha$ID + 1])
  p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = alpha.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Alpha - inlabru", x = "x", y='')
  
  data.beta = cbind(res.inlabru$summary.random$beta,
                    beta.true = underlying.effects$beta.true[res.inlabru$summary.random$beta$ID + 1])
  p.beta <- ggplot(data = data.beta, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = beta.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "x", y = "beta", title = "Beta - inlabru")
  
  data.kappa = cbind(res.inlabru$summary.random$kappa,
                     kappa.drifted = underlying.effects$kappa.drifted[res.inlabru$summary.random$kappa$ID + 1])
  p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = kappa.drifted, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "t", y = "kappa", title = "Kappa - inlabru")
  
  data.fixed  <- data.frame(res.inlabru$marginals.fixed)
  
  p.intercept <- ggplot(data.fixed) + 
    geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Estimated", fill="Estimated")) + 
    geom_vline(aes(xintercept = underlying.effects$age.intercept.true, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of phi", y = " ", title = "Intercept - inlabru")
  
  data.gamma = cbind(res.inlabru$summary.random$gamma,
                     gamma.true = underlying.effects$gamma.true[res.inlabru$summary.random$gamma$ID + nx])
  p.gamma <- ggplot(data = data.gamma, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = gamma.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "t", y = "gamma", title = "Gamma - inlabru")
  
  p.random.effects <- (p.intercept | p.alpha ) / (p.beta | p.kappa | p.gamma) + 
    plot_layout(guides = "collect")
  
  if(save){
    save.figure(p.random.effects, name = "random_effects_inlabru", path = path.to.storage)
  }
  
  data.eta <- data.frame(eta.sim = res.inlabru$summary.linear.predictor$mean[1:length(obs$eta)]) %>%
    mutate(true.eta = obs$eta) %>%
    mutate(xt = obs$xt, x = obs$x, t = obs$t)
  
  p.eta <- ggplot(data = data.eta) +
    geom_point(aes(x = eta.sim, y = true.eta), color = palette[1]) + 
    labs(x="Estimated eta", y="True value for eta", title = "Eta")
  
  p.eta.2 <- ggplot(data = data.eta) +
    geom_line(aes(x=xt, y = eta.sim, color="Estimated")) +
    geom_line(aes(x=xt, y = true.eta, color="True")) +
    labs(x=" ", y="Eta", title="Eta- inlabru")
  
  p.eta.t <- ggplot(data = data.eta) + 
    geom_line(aes(x = x, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = x, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.eta.x <- ggplot(data = data.eta) + 
    geom_line(aes(x = t, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = t, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each age") + 
    facet_wrap(~x)
  
  p.eta.xt <- (p.eta | p.eta.2) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.xt, name="eta_xt_inlabru", path=path.to.storage)
  }
  
  p.eta.facet <- (p.eta.t | p.eta.x) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.facet, name="eta_facet_inlabru", path=path.to.storage)
  }
  
  plots <- list(p.alpha = p.alpha, p.beta = p.beta, p.kappa = p.kappa,
                p.eta = p.eta, p.gamma = p.gamma,
                p.eta.2 = p.eta.2, p.eta.x = p.eta.x, 
                p.eta.t = p.eta.t, p.intercept = p.intercept,
                p.random.effects = p.random.effects,
                p.eta.all = p.eta.xt,
                p.eta.facet = p.eta.facet)
  
  summaries <- list(data.alpha = data.alpha,
                    data.beta = data.beta,
                    data.kappa = data.kappa,
                    data.eta = data.eta,
                    data.fixed = data.fixed,
                    data.gamma = data.gamma,
                    intercept = res.inlabru$summary.fixed$mean[1],
                    phi = res.inlabru$summary.fixed$mean[2])
  
  return(list(plots = plots, summaries = summaries))
}

underlying.effects <- configuration.v17.3()
obs <- underlying.effects$obs

# run witout extra constraints on gamma
res.inlabru.1 <- inlabru.undrifted.period.cohort.2(obs)

plots.1 <- plot.res(res.inlabru.1, underlying.effects)
plots.1$plots$p.random.effects
plots.1$plots$p.eta.facet
plots.1$plots$p.eta.all

# run with extra constraints on gamma
res.inlabru.2 <- inlabru.undrifted.period.cohort.2.gamma.extraconstr(obs)

plots.2 <- plot.res(res.inlabru.2, underlying.effects)
plots.2$plots$p.random.effects
plots.2$plots$p.eta.facet
plots.2$plots$p.eta.all



