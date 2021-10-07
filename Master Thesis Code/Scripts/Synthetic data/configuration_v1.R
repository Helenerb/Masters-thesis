library(INLA)
library(inlabru)
library(rstan)
library(tidyverse)
library(ggplot2)

configuration.v1 <- function(){
  #   ----   Setiting seed for reproductiveness   ----
  
  seed = 324
  set.seed(seed)
  
  #   ----   Defining data structure   ----
  
  nx = 80   # number of age groups
  nt = 100  # number of time steps
  
  cohorts <- (-79):99  # cohort indices, c = t-x
  nc = length(cohorts)  # number of cohorts
  
  #   ----   Definining underlying effects   ----
  
  alpha.true = 3.9*cos(((1:nx + 30)* pi)/50) - 10
  
  phi.true = -0.25
  
  tau.beta.true = 150
  
  beta.true = rnorm(nx, mean=0, sd=sqrt(1/tau.beta.true))
  beta.true = beta.true - mean(beta.true) + 1/nx
  
  tau.kappa.true = 240
  kappa.true.increments = rnorm(nt-1, mean=0, sd=sqrt(1/tau.kappa.true))
  kappa.true.increments = kappa.true.increments - mean(kappa.true.increments)
  kappa.true = rep(0, nt)
  for (idx in 2:nt){
    kappa.true[idx] = kappa.true[idx - 1] + kappa.true.increments[idx-1]
  }
  kappa.true = kappa.true - mean(kappa.true)
  
  tau.gamma.true = 250
  gamma.true = rnorm(nc, mean = 0, sd=sqrt(1/tau.gamma.true))
  gamma.true = gamma.true - mean(gamma.true)
  for (idx in 2:nc) {
    gamma.true[idx] = gamma.true[idx - 1] + gamma.true[idx]
  }
  # increase risk for cohorts 40-50:
  gamma.true[40:50] = gamma.true[40:50] + rnorm(11, mean=0.4, sd = sqrt(1/(tau.gamma.true)))
  
  gamma.true = gamma.true - mean(gamma.true)  # center around zero
  
  tau.epsilon.true = 1000
  
  at.risk = 10**6/nx
  
  obs <- data.frame(expand.grid(x=1:nx, t=1:nt)) %>%
    mutate(c = t-x) %>%
    mutate(alpha = alpha.true[x]) %>%
    mutate(beta = beta.true[x]) %>%
    mutate(kappa = kappa.true[t]) %>%
    mutate(phi.t = phi.true*t) %>%
    mutate(gamma = gamma.true[c +80]) %>%
    mutate(epsilon = rnorm(length(x), mean = 0, sd = sqrt(1/tau.epsilon.true))) %>%
    mutate(E = at.risk) %>%
    mutate(eta = alpha + beta*phi.t + beta*kappa + gamma + epsilon) %>%
    mutate(Y = rpois(length(x), E*exp(eta))) %>%
    mutate(xt = seq_along(x)) %>%
    mutate(x.c = x, t.c = t)
  
  underlying.effects <- list(
    obs = obs,
    alpha.true = alpha.true,
    beta.true = beta.true,
    kappa.true = kappa.true,
    phi.true = phi.true,
    gamma.true = gamma.true,
    at.risk = at.risk
  )
  return(underlying.effects)
}

#   ----   Updating configuration to include better period effects in eta   ----

configuration.v2 <- function(){
  #   ----   Setiting seed for reproductiveness   ----
  
  seed = 324
  set.seed(seed)
  
  #   ----   Defining data structure   ----
  
  nx = 80   # number of age groups
  nt = 100  # number of time steps
  
  cohorts <- (-79):99  # cohort indices, c = t-x
  nc = length(cohorts)  # number of cohorts
  
  #   ----   Definining underlying effects   ----
  
  alpha.true = 3.9*cos(((1:nx + 30)* pi)/50) - 8
  
  phi.true = -1
  
  tau.beta.true = 30000
  
  beta.true = rnorm(nx, mean=0, sd=sqrt(1/tau.beta.true))
  beta.true = beta.true - mean(beta.true) + 1/nx
  
  tau.kappa.true = 240
  kappa.true.increments = rnorm(nt-1, mean=0, sd=sqrt(1/tau.kappa.true))
  kappa.true.increments = kappa.true.increments - mean(kappa.true.increments)
  kappa.true = rep(0, nt)
  for (idx in 2:nt){
    kappa.true[idx] = kappa.true[idx - 1] + kappa.true.increments[idx-1]
  }
  kappa.true = kappa.true - mean(kappa.true)
  
  tau.gamma.true = 250
  gamma.true = rnorm(nc, mean = 0, sd=sqrt(1/tau.gamma.true))
  gamma.true = gamma.true - mean(gamma.true)
  for (idx in 2:nc) {
    gamma.true[idx] = gamma.true[idx - 1] + gamma.true[idx]
  }
  # increase risk for cohorts 40-50:
  gamma.true[40:50] = gamma.true[40:50] + rnorm(11, mean=0.4, sd = sqrt(1/(tau.gamma.true)))
  
  gamma.true = gamma.true - mean(gamma.true)  # center around zero
  
  tau.epsilon.true = 1000
  
  at.risk = 10**6/nx
  
  obs <- data.frame(expand.grid(x=1:nx, t=1:nt)) %>%
    mutate(c = t-x) %>%
    mutate(alpha = alpha.true[x]) %>%
    mutate(beta = beta.true[x]) %>%
    mutate(kappa = kappa.true[t]) %>%
    mutate(phi.t = phi.true*t) %>%
    mutate(gamma = gamma.true[c +80]) %>%
    mutate(epsilon = rnorm(length(x), mean = 0, sd = sqrt(1/tau.epsilon.true))) %>%
    mutate(E = at.risk) %>%
    mutate(eta = alpha + beta*phi.t + beta*kappa + gamma + epsilon) %>%
    mutate(Y = rpois(length(x), E*exp(eta))) %>%
    mutate(xt = seq_along(x)) %>%
    mutate(x.c = x, t.c = t)
  
  underlying.effects <- list(
    obs = obs,
    alpha.true = alpha.true,
    beta.true = beta.true,
    kappa.true = kappa.true,
    phi.true = phi.true,
    gamma.true = gamma.true,
    at.risk = at.risk
  )
  return(underlying.effects)
}

configuration.v3 <- function(){
  # simpler version of v2 without cohort effects
  
  #   ----   Setiting seed for reproductiveness   ----
  
  seed = 324
  set.seed(seed)
  
  #   ----   Defining data structure   ----
  
  nx = 80   # number of age groups
  nt = 100  # number of time steps
  
  #   ----   Definining underlying effects   ----
  
  alpha.true = 3.9*cos(((1:nx + 30)* pi)/50) - 6
  
  phi.true = -1
  
  tau.beta.true = 30000
  
  beta.true = rnorm(nx, mean=0, sd=sqrt(1/tau.beta.true))
  beta.true = beta.true - mean(beta.true) + 1/nx
  
  tau.kappa.true = 240
  kappa.true.increments = rnorm(nt-1, mean=0, sd=sqrt(1/tau.kappa.true))
  kappa.true.increments = kappa.true.increments - mean(kappa.true.increments)
  kappa.true = rep(0, nt)
  for (idx in 2:nt){
    kappa.true[idx] = kappa.true[idx - 1] + kappa.true.increments[idx-1]
  }
  kappa.true = kappa.true - mean(kappa.true)
  
  tau.epsilon.true = 1000
  
  at.risk = 10**6/nx
  
  obs <- data.frame(expand.grid(x=1:nx, t=1:nt)) %>%
    mutate(c = t-x) %>%
    mutate(alpha = alpha.true[x]) %>%
    mutate(beta = beta.true[x]) %>%
    mutate(kappa = kappa.true[t]) %>%
    mutate(phi.t = phi.true*t) %>%
    mutate(epsilon = rnorm(length(x), mean = 0, sd = sqrt(1/tau.epsilon.true))) %>%
    mutate(E = at.risk) %>%
    mutate(eta = alpha + beta*phi.t + beta*kappa + epsilon) %>%
    mutate(Y = rpois(length(x), E*exp(eta))) %>%
    mutate(xt = seq_along(x)) %>%
    mutate(x.c = x, t.c = t)
  
  underlying.effects <- list(
    obs = obs,
    alpha.true = alpha.true,
    beta.true = beta.true,
    kappa.true = kappa.true,
    phi.true = phi.true,
    at.risk = at.risk
  )
  return(underlying.effects)
}

configuration.v4 <- function(){
  # version of v3 with more erratic kappa -> kappa taking more extreme values
  
  #   ----   Setiting seed for reproductiveness   ----
  
  seed = 324
  set.seed(seed)
  
  #   ----   Defining data structure   ----
  
  nx = 80   # number of age groups
  nt = 100  # number of time steps
  
  #   ----   Definining underlying effects   ----
  
  alpha.true = 3.9*cos(((1:nx + 30)* pi)/50) - 4
  
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
  kappa.true = kappa.true - mean(kappa.true)
  
  tau.epsilon.true = 1000
  
  at.risk = 10**6/nx
  
  obs <- data.frame(expand.grid(x=1:nx, t=1:nt)) %>%
    mutate(c = t-x) %>%
    mutate(alpha = alpha.true[x]) %>%
    mutate(beta = beta.true[x]) %>%
    mutate(kappa = kappa.true[t]) %>%
    mutate(phi.t = phi.true*(t)) %>%
    mutate(epsilon = rnorm(length(x), mean = 0, sd = sqrt(1/tau.epsilon.true))) %>%
    mutate(E = at.risk) %>%
    mutate(eta = alpha + beta*phi.t + beta*kappa + epsilon) %>%
    mutate(Y = rpois(length(x), E*exp(eta))) %>%
    mutate(xt = seq_along(x)) %>%
    mutate(x.c = x, t.c = t)
  
  underlying.effects <- list(
    obs = obs,
    alpha.true = alpha.true,
    beta.true = beta.true,
    kappa.true = kappa.true,
    phi.true = phi.true,
    at.risk = at.risk
  )
  return(underlying.effects)
}

configuration.v5 <- function(){
  # version of v4 with 0-indexed t and x
  # with alpha centered around zero and separate age intercept
  
  #   ----   Setiting seed for reproductiveness   ----
  
  seed = 324
  set.seed(seed)
  
  #   ----   Defining data structure   ----
  
  nx = 80   # number of age groups
  nt = 100  # number of time steps
  
  #   ----   Definining underlying effects   ----
  
  alpha.true = 3.9*cos(((1:nx + 30)* pi)/50)
  alpha.true = alpha.true - mean(alpha.true)
  
  age.intercept.true = - 4
  
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
  kappa.true = kappa.true - mean(kappa.true)
  
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
    alpha.true = alpha.true,
    age.intercept.true= age.intercept.true,
    beta.true = beta.true,
    kappa.true = kappa.true,
    phi.true = phi.true,
    at.risk = at.risk
  )
  return(underlying.effects)
}

configuration.v9 <- function(){
  # version of v7 with fewer period and age steps
  
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
  kappa.true = kappa.true - mean(kappa.true)
  
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
    tau.epsilon.true= tau.epsilon.true
  )
  return(underlying.effects)
}


configuration.v6 <- function(){
  # version of v5 with cohort effect, 0- indexed 
  
  #   ----   Setiting seed for reproductiveness   ----
  
  seed = 324
  set.seed(seed)
  
  #   ----   Defining data structure   ----
  
  nx = 80   # number of age groups
  nt = 100  # number of time steps
  
  cohorts <- (-79):99  # cohort indices, c = t-x
  nc = length(cohorts)  # number of cohorts
  
  #   ----   Definining underlying effects   ----
  
  alpha.true = 3.9*cos(((1:nx + 30)* pi)/50) - 4
  
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
  kappa.true = kappa.true - mean(kappa.true)
  
  tau.gamma.true = 250
  gamma.true = rnorm(nc, mean = 0, sd=sqrt(1/tau.gamma.true))
  gamma.true = gamma.true - mean(gamma.true)
  for (idx in 2:nc) {
    gamma.true[idx] = gamma.true[idx - 1] + gamma.true[idx]
  }
  # increase risk for cohorts 40-50:
  gamma.true[40:50] = gamma.true[40:50] + rnorm(11, mean=0.4, sd = sqrt(1/(tau.gamma.true)))
  
  gamma.true = gamma.true - mean(gamma.true)  # center around zero
  
  tau.epsilon.true = 1000
  
  at.risk = 10**6/nx
  
  obs <- data.frame(expand.grid(x=1:nx, t=1:nt)) %>%
    mutate(alpha = alpha.true[x]) %>%
    mutate(beta = beta.true[x]) %>%
    mutate(kappa = kappa.true[t]) %>%
    mutate(phi.t = phi.true*(t-1)) %>%
    #mutate(phi.t = phi.true*(t)) %>%
    mutate(xt = seq_along(x)) %>%
    mutate(x = x - 1, t = t - 1, xt = xt - 1) %>%
    mutate(c = t-x) %>%
    mutate(gamma = gamma.true[c + 80]) %>%
    mutate(epsilon = rnorm(length(x), mean = 0, sd = sqrt(1/tau.epsilon.true))) %>%
    mutate(E = at.risk) %>%
    mutate(eta = alpha + beta*phi.t + beta*kappa + gamma + epsilon) %>%
    mutate(Y = rpois(length(x), E*exp(eta))) %>%
    #mutate(c = c - 1) %>%
    mutate(x.c = x, t.c = t)
  
  underlying.effects <- list(
    obs = obs,
    alpha.true = alpha.true,
    beta.true = beta.true,
    kappa.true = kappa.true,
    phi.true = phi.true,
    gamma.true = gamma.true,
    at.risk = at.risk
  )
  return(underlying.effects)
}

configuration.v7 <- function(){
  # version of v6 with cohort effect, 0-indexed, redefining alpha as centered
  # around zero with separate intercept
  
  #   ----   Setiting seed for reproductiveness   ----
  
  seed = 324
  set.seed(seed)
  
  #   ----   Defining data structure   ----
  
  nx = 80   # number of age groups
  nt = 100  # number of time steps
  
  cohorts <- (-79):99  # cohort indices, c = t-x
  nc = length(cohorts)  # number of cohorts
  
  #   ----   Definining underlying effects   ----
  
  alpha.true = 3.9*cos(((1:nx + 30)* pi)/50)
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
  kappa.true = kappa.true - mean(kappa.true)
  
  tau.gamma.true = 250
  gamma.true = rnorm(nc, mean = 0, sd=sqrt(1/tau.gamma.true))
  gamma.true = gamma.true - mean(gamma.true)
  for (idx in 2:nc) {
    gamma.true[idx] = gamma.true[idx - 1] + gamma.true[idx]
  }
  # increase risk for cohorts 40-50:
  gamma.true[40:50] = gamma.true[40:50] + rnorm(11, mean=0.4, sd = sqrt(1/(tau.gamma.true)))
  
  gamma.true = gamma.true - mean(gamma.true)  # center around zero
  
  tau.epsilon.true = 1000
  
  at.risk = 10**6/nx
  
  obs <- data.frame(expand.grid(x=1:nx, t=1:nt)) %>%
    mutate(alpha = alpha.true[x]) %>%
    mutate(age.intercept = age.intercept.true) %>%
    mutate(beta = beta.true[x]) %>%
    mutate(kappa = kappa.true[t]) %>%
    mutate(phi.t = phi.true*(t-1)) %>%
    mutate(xt = seq_along(x)) %>%
    mutate(x = x - 1, t = t - 1, xt = xt - 1) %>%
    mutate(c = t-x) %>%
    mutate(gamma = gamma.true[c + 80]) %>%
    mutate(epsilon = rnorm(length(x), mean = 0, sd = sqrt(1/tau.epsilon.true))) %>%
    mutate(E = at.risk) %>%
    mutate(eta = age.intercept + alpha + beta*phi.t + beta*kappa + gamma + epsilon) %>%
    mutate(Y = rpois(length(x), E*exp(eta))) %>%
    mutate(x.c = x, t.c = t)
  
  underlying.effects <- list(
    obs = obs,
    alpha.true = alpha.true,
    beta.true = beta.true,
    kappa.true = kappa.true,
    phi.true = phi.true,
    gamma.true = gamma.true,
    at.risk = at.risk,
    age.intercept.true = age.intercept.true
  )
  return(underlying.effects)
}

configuration.v8 <- function(){
  # version of v7 with different settings for phi and age.intercept, to check if 
  # inlabru really does Int = age.intercept - phi
  
  #   ----   Setiting seed for reproductiveness   ----
  
  seed = 324
  set.seed(seed)
  
  #   ----   Defining data structure   ----
  
  nx = 80   # number of age groups
  nt = 100  # number of time steps
  
  cohorts <- (-79):99  # cohort indices, c = t-x
  nc = length(cohorts)  # number of cohorts
  
  #   ----   Definining underlying effects   ----
  
  alpha.true = 3.9*cos(((1:nx + 30)* pi)/50)
  alpha.true = alpha.true - mean(alpha.true)  # center around zero
  
  age.intercept.true = -4
  #age.intercept.true = 0
  
  phi.true = -1.8
  
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
  kappa.true = kappa.true - mean(kappa.true)
  
  tau.gamma.true = 250
  gamma.true = rnorm(nc, mean = 0, sd=sqrt(1/tau.gamma.true))
  gamma.true = gamma.true - mean(gamma.true)
  for (idx in 2:nc) {
    gamma.true[idx] = gamma.true[idx - 1] + gamma.true[idx]
  }
  # increase risk for cohorts 40-50:
  gamma.true[40:50] = gamma.true[40:50] + rnorm(11, mean=0.4, sd = sqrt(1/(tau.gamma.true)))
  
  gamma.true = gamma.true - mean(gamma.true)  # center around zero
  
  tau.epsilon.true = 1000
  
  at.risk = 10**6/nx
  
  obs <- data.frame(expand.grid(x=1:nx, t=1:nt)) %>%
    mutate(alpha = alpha.true[x]) %>%
    mutate(age.intercept = age.intercept.true) %>%
    mutate(beta = beta.true[x]) %>%
    mutate(kappa = kappa.true[t]) %>%
    mutate(phi.t = phi.true*(t-1)) %>%
    mutate(xt = seq_along(x)) %>%
    mutate(x = x - 1, t = t - 1, xt = xt - 1) %>%
    mutate(c = t-x) %>%
    mutate(gamma = gamma.true[c + 80]) %>%
    mutate(epsilon = rnorm(length(x), mean = 0, sd = sqrt(1/tau.epsilon.true))) %>%
    mutate(E = at.risk) %>%
    mutate(eta = age.intercept + alpha + beta*phi.t + beta*kappa + gamma + epsilon) %>%
    mutate(Y = rpois(length(x), E*exp(eta))) %>%
    mutate(x.c = x, t.c = t)
  
  underlying.effects <- list(
    obs = obs,
    alpha.true = alpha.true,
    beta.true = beta.true,
    kappa.true = kappa.true,
    phi.true = phi.true,
    gamma.true = gamma.true,
    at.risk = at.risk,
    age.intercept.true = age.intercept.true
  )
  return(underlying.effects)
}




plot.underlying.effects <- function(u.e){
  
  palette.basis <- c('#70A4D4', '#ECC64B', '#93AD80', '#da9124', '#696B8D',
                     '#3290c1',
                     '#5d8060', '#D7B36A', '#826133', '#A85150')
  
  p.1 <- ggplot(data = u.e$obs) + geom_line(aes(x=xt, y = eta))
  p.2 <- ggplot(data = u.e$obs) + geom_point(aes(x = xt, y = Y))  
  p.3 <- ggplot(data = u.e$obs) + geom_point(aes(x = xt, y = E*exp(eta)))
  p.4 <- ggplot(data = u.e$obs) + geom_line(aes(x = xt, y = exp(eta)))
  p.5 <- ggplot(data = u.e$obs) + geom_point(aes(x = xt, y = exp(eta)))
  p.6 <- ggplot(data = u.e$obs) + geom_point(aes(x = xt, y = eta))
  
  #   ----   plot by time   ----
  p.7 <- ggplot(data = u.e$obs) + geom_point(aes(x = x, y = eta)) + facet_wrap(~t)
  
  #   ----   plot by age   ----
  p.8 <- ggplot(data = u.e$obs) + geom_point(aes(x = t, y = eta)) + facet_wrap(~x) +
    labs(title = "Eta, for each age")
  
  p.beta <- ggplot(data = u.e$obs) + geom_point(aes(x = x, y = beta)) +
    labs(title = "Underlying beta", y = " ", x = "x")
  
  p.gamma <- ggplot(data = u.e$obs) + geom_point(aes(x = c, y = gamma)) +
    labs(title = "Underlying gamma", y = " ", x = "t-x")
  
  obs.stepwise.eta <- u.e$obs %>% 
    mutate(beta.phi = beta*phi.t) %>%
    mutate(beta.phi.kappa = beta.phi + beta*kappa) %>%
    mutate(alpha.beta.phi.kappa = alpha + beta.phi.kappa) %>%
    mutate(alpha.beta.phi.kappa.gamma = alpha.beta.phi.kappa + gamma)
  
  p.eta.stepwise <- ggplot(data = obs.stepwise.eta, aes(x=x)) +
    geom_point(aes(y = eta, color = "Eta"), size = 1, shape = 1) + 
    geom_point(aes(y = alpha.beta.phi.kappa.gamma, color = "Alpha + beta*phi.t + beta*kappa + gamma"), size = 1, shape = 2) +
    geom_point(aes(y = alpha.beta.phi.kappa, color = "Alpha + beta*phi.t + beta*kappa"), size = 1, shape = 3) +
    geom_point(aes(y = beta.phi.kappa, color = "Beta*phi.t + beta*kappa"), size = 1, shape = 4) + 
    geom_point(aes(y = beta.phi, color = "Beta*phi"), size = 1, shape = 5) + 
    #geom_point(aes(y = alpha, color = "Alpha"), shape = 1) +
    scale_color_manual(name = " ", values = palette.basis)+
    facet_wrap(~t)
  
  p.eta.stepwise.t <- ggplot(data = obs.stepwise.eta, aes(x=t)) +
    geom_point(aes(y = eta, color = "Eta"), size = 0.5, shape = 1) + 
    geom_point(aes(y = alpha.beta.phi.kappa.gamma, color = "Alpha + beta*phi + beta*kappa + gamma"), size = 0.5, shape = 2) +
    geom_point(aes(y = alpha.beta.phi.kappa, color = "Alpha + beta*phi + beta*kappa"), size = 0.5, shape = 3) +
    geom_point(aes(y = beta.phi.kappa, color = "Beta*phi.t + beta*kappa"), size = 0.5, shape = 4) + 
    geom_point(aes(y = beta.phi, color = "Beta*phi.t"), size = 0.5, shape = 5) + 
    #geom_point(aes(y = alpha, color = "Alpha"), shape = 2) +
    #geom_point(aes(y = kappa + phi.t, color = "kappa + phi*t"), shape = 4, size = 1) +
    #geom_point(aes(y = beta, color = "Beta"), shape = 5, size = 1.5) +
    scale_color_manual(name = " ", values = palette.basis)+
    facet_wrap(~x)
  
  plots <- list(
    p.1 = p.1, p.2 = p.2, p.3 = p.3, p.4 = p.4, p.5 = p.5, p.6 = p.6, p.7 = p.7,
    p.8 = p.8, p.beta = p.beta, p.eta.stepwise = p.eta.stepwise,
    p.eta.stepwise.t = p.eta.stepwise.t, p.gamma = p.gamma)
  
  return(plots)
}

plot.underlying.effects.age.period <- function(u.e){
  palette.basis <- c('#70A4D4', '#ECC64B', '#93AD80', '#da9124', '#696B8D',
                     '#3290c1',
                     '#5d8060', '#D7B36A', '#826133', '#A85150')
  
  p.1 <- ggplot(data = u.e$obs) + geom_line(aes(x=xt, y = eta))
  p.2 <- ggplot(data = u.e$obs) + geom_point(aes(x = xt, y = Y))  
  p.3 <- ggplot(data = u.e$obs) + geom_point(aes(x = xt, y = E*exp(eta)))
  p.4 <- ggplot(data = u.e$obs) + geom_line(aes(x = xt, y = exp(eta)))
  p.5 <- ggplot(data = u.e$obs) + geom_point(aes(x = xt, y = exp(eta)))
  p.6 <- ggplot(data = u.e$obs) + geom_point(aes(x = xt, y = eta))
  
  #   ----   plot by time   ----
  p.7 <- ggplot(data = u.e$obs) + geom_point(aes(x = x, y = eta)) + facet_wrap(~t)
  
  #   ----   plot by age   ----
  p.8 <- ggplot(data = u.e$obs) + geom_point(aes(x = t, y = eta)) + facet_wrap(~x) +
    labs(title = "Eta, for each age")
  
  p.beta <- ggplot(data = u.e$obs) + geom_point(aes(x = x, y = beta)) +
    labs(title = "Underlying beta", y = " ", x = "x")
  
  p.kappa <- ggplot(data = u.e$obs) + geom_point(aes(x = t, y = kappa)) +
    labs(title = "Kappa")
  
  p.kappa.phi <- ggplot(data = u.e$obs) + geom_point(aes(x = t, y = kappa + phi.t)) + 
    labs(title = "Kappa + phi*t")
  
  obs.stepwise.eta <- u.e$obs %>% 
    mutate(beta.phi = beta*phi.t) %>%
    mutate(beta.phi.kappa = beta.phi + beta*kappa) %>%
    mutate(alpha.beta.phi.kappa = alpha + beta.phi.kappa)
  
  p.eta.stepwise <- ggplot(data = obs.stepwise.eta, aes(x=x)) +
    geom_point(aes(y = eta, color = "Eta"), size = 1, shape = 1) + 
    geom_point(aes(y = alpha.beta.phi.kappa, color = "Alpha + beta*phi.t + beta*kappa"), size = 1, shape = 3) +
    geom_point(aes(y = beta.phi.kappa, color = "Beta*phi.t + beta*kappa"), size = 1, shape = 4) + 
    geom_point(aes(y = beta.phi, color = "Beta*phi"), size = 1, shape = 5) + 
    #geom_point(aes(y = alpha, color = "Alpha"), shape = 1) +
    scale_color_manual(name = " ", values = palette.basis)+
    facet_wrap(~t)
  
  p.eta.stepwise.t <- ggplot(data = obs.stepwise.eta, aes(x=t)) +
    geom_point(aes(y = eta, color = "Eta"), size = 0.5, shape = 1) + 
    geom_point(aes(y = alpha.beta.phi.kappa, color = "Alpha + beta*phi + beta*kappa"), size = 0.5, shape = 3) +
    geom_point(aes(y = beta.phi.kappa, color = "Beta*phi.t + beta*kappa"), size = 0.5, shape = 4) + 
    geom_point(aes(y = beta.phi, color = "Beta*phi.t"), size = 0.5, shape = 5) + 
    #geom_point(aes(y = alpha, color = "Alpha"), shape = 2) +
    #geom_point(aes(y = kappa + phi.t, color = "kappa + phi*t"), shape = 4, size = 1) +
    #geom_point(aes(y = beta, color = "Beta"), shape = 5, size = 1.5) +
    scale_color_manual(name = " ", values = palette.basis)+
    facet_wrap(~x)
  
  plots <- list(
    p.1 = p.1, p.2 = p.2, p.3 = p.3, p.4 = p.4, p.5 = p.5, p.6 = p.6, p.7 = p.7,
    p.8 = p.8, p.beta = p.beta, p.eta.stepwise = p.eta.stepwise,
    p.eta.stepwise.t = p.eta.stepwise.t, p.kappa.phi = p.kappa.phi,
    p.kappa = p.kappa)
  
  return(plots)
}

#u.e.v2 <- configuration.v2()

# u.e.v4 <- configuration.v4() 
# 

u.e.v7 <- configuration.v9()
plots.v2 <- plot.underlying.effects.age.period(u.e.v7)
plots.v2$p.8
plots.v2$p.beta
plots.v2$p.eta.stepwise
plots.v2$p.eta.stepwise.t
plots.v2$p.1
plots.v2$p.2
plots.v2$p.6
plots.v2$p.kappa.phi
plots.v2$p.kappa
plots.v2$p.gamma
