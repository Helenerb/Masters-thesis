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

lung.cancer.male.above.45 <- lung.cancer.male %>%
  filter(age.int >= 45)

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

res.lung.lc.m.a45 <- inlabru.rw2.lc.2(lung.cancer.male.above.45)

inlabru.tau.alpha.a45 <- res.lung.lc.m.a45$summary.hyperpar$mean[1] 
inlabru.tau.beta.a45 <- res.lung.lc.m.a45$summary.hyperpar$mean[2]
inlabru.tau.kappa.a45 <- res.lung.lc.m.a45$summary.hyperpar$mean[3]
inlabru.tau.epsilon.a45 <- res.lung.lc.m.a45$summary.hyperpar$mean[4]

inlabru.intercept.a45 <- res.lung.lc.m.a45$summary.fixed$mean[1]

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

#   ----   v7   ----

# attempt to find a configuration that is of somewhat comparable size to real data
set.seed(2) 
alpha.9 <- rw1(inlabru.tau.alpha.a45, 9); alpha.9 = alpha.9 - mean(alpha.9)
beta.9 <- rnorm(9, sd = sqrt(1/inlabru.tau.beta.a45)); beta.9 = beta.9 - mean(beta.9) + 1/9
kappa.9 <- rw2(inlabru.tau.kappa.a45, 18); kappa.9 = kappa.9 - mean(kappa.9)
epsilon.9 <- rnorm(9*18, sd = sqrt(1/inlabru.tau.epsilon.a45))

obs.9 <- data.frame(x = lung.cancer.male.above.45$x, t = lung.cancer.male.above.45$t, E = lung.cancer.male.above.45$E) %>%
  mutate(age.int = lung.cancer.male.above.45$age.int, year = lung.cancer.male.above.45$year) %>%
  mutate(xt = seq_along(x)) %>%
  mutate(x.c = x) %>%
  mutate(alpha = alpha.9[x - 9 + 1]) %>%
  mutate(beta = beta.9[x - 9 + 1]) %>%
  mutate(kappa = kappa.9[t + 1]) %>%
  mutate(intercept = inlabru.intercept.a45) %>%
  mutate(epsilon = epsilon.9[xt]) %>%
  mutate(eta = intercept + alpha + beta*kappa + epsilon) %>%
  mutate(e.exp.eta = E*exp(eta)) %>%
  mutate(Y = rpois(length(x), E*exp(eta))) %>%
  mutate(mr = Y/E) %>%
  mutate(tau.alpha = inlabru.tau.alpha.a45) %>%
  mutate(tau.beta = inlabru.tau.beta.a45) %>%
  mutate(tau.kappa = inlabru.tau.kappa.a45) %>%
  mutate(tau.epsilon = inlabru.tau.epsilon.a45)

# write data to csv, and read data from csv.

#write.csv(obs.9, "synthetic_male_lung_9.csv")
# obs.9 <- read.csv("Data/synthetic_male_lung_9.csv")
# obs.9 <- obs.9 %>% mutate(year.str = as.character(year))



