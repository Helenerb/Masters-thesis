# independent script for producing v4 version of synthetic lung cancer data:

library("tidyverse")
library("ggplot2")

load("population-germany.Rda")

load("lungCancer-germany.Rda")
lung.cancer <- cancer.data

lung.cancer <- lung.cancer %>%
  mutate(female.mr = female/female.t, male.mr = male/male.t)

lung.cancer.male <- lung.cancer %>%
  select(x, x.c, t, xt, cohort, c, age, age.int, year, birth.year, male.t, male, male.mr) %>%
  mutate(E = male.t, Y = male, `mortality rate` = male.mr)

# function running inlabru with Poisson lee-carter model:
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

res.lung.lc.m <- inlabru.rw2.lc.2(lung.cancer.male)

inlabru.tau.alpha <- res.lung.lc.m$summary.hyperpar$mean[1]
inlabru.tau.beta <- res.lung.lc.m$summary.hyperpar$mean[2]
inlabru.tau.kappa <- res.lung.lc.m$summary.hyperpar$mean[3]
inlabru.tau.epsilon <- res.lung.lc.m$summary.hyperpar$mean[4]

inlabru.intercept <- res.lung.lc.m$summary.fixed$mean[1]

rw2 <- function(tau, nx){
  #' samples a random walk of order two
  #' 
  #' @param tau (float): precision of random walk
  #' @param nx (int): length of random walk 
  
  increments <- rnorm(nx, mean = 0, sd = sqrt(1/tau))
  random.walk <- increments
  random.walk[2] <- random.walk[1] + increments[2]
  for(i in 3:nx){
    random.walk[i] <- 2*random.walk[i-1]  - random.walk[i-2] + increments[i]
  }
  return(random.walk)
}

rw1 <- function(tau, nx){
  #' samples a random walk of order one
  #' 
  #' @param tau (float): precision of random walk
  #' @param nx (int): length of random walk 
  
  increments <- rnorm(nx, mean = 0, sd = sqrt(1/tau))
  random.walk <- increments
  for(i in 2:nx){
    random.walk[i] <- random.walk[i-1] + increments[i]
  }
  return(random.walk)
}

#   ----   v4   ----

# attempt to find a configuration that is of somewhat comparable size to real data
set.seed(14) # quite good for alpha and eta, we try this. Note - still not very steep kappa
alpha.4 <- rw1(inlabru.tau.alpha, 18); alpha.4 = alpha.4 - mean(alpha.4)
beta.4 <- rnorm(18, sd = sqrt(1/inlabru.tau.beta)); beta.4 = beta.4 - mean(beta.4) + 1/18
kappa.4 <- rw2(inlabru.tau.kappa, 18); kappa.4 = kappa.4 - mean(kappa.4)
epsilon.4 <- rnorm(18*18, sd = sqrt(1/inlabru.tau.epsilon))

obs.4 <- data.frame(x = lung.cancer.male$x, t = lung.cancer.male$t, xt = lung.cancer.male$xt, E = lung.cancer.male$E) %>%
  mutate(age.int = lung.cancer.male$age.int, year = lung.cancer.male$year) %>%
  mutate(x.c = x) %>%
  mutate(alpha = alpha.4[x + 1]) %>%
  mutate(beta = beta.4[x + 1]) %>%
  mutate(kappa = kappa.4[t + 1]) %>%
  mutate(intercept = inlabru.intercept) %>%
  mutate(epsilon = epsilon.4[xt + 1]) %>%
  mutate(eta = intercept + alpha + beta*kappa + epsilon) %>%
  mutate(Y = rpois(length(x), E*exp(eta))) %>%
  mutate(mr = Y/E) %>%
  mutate(tau.alpha = inlabru.tau.alpha) %>%
  mutate(tau.beta = inlabru.tau.beta) %>%
  mutate(tau.kappa = inlabru.tau.kappa) %>%
  mutate(tau.epsilon = inlabru.tau.epsilon)

ggplot(obs.4) + geom_line(aes(x = x, y = eta, color = year))
ggplot(obs.4) + geom_line(aes(x = x, y = mr, color = year))
ggplot(obs.4) + geom_line(aes(x = x, y = alpha))
ggplot(obs.4) + geom_line(aes(x = x, y = beta))
ggplot(obs.4) + geom_line(aes(x = t, y = kappa))

#write.csv(obs.4, "Data/synthetic_male_lung_4.csv")
#obs.4 <- read.csv("Data/synthetic_male_lung_4.csv")


