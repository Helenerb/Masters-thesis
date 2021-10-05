# functions for setting up and running analysis with inlabru.
# All functions take an obs-dataframe containing data for analysis as argument
# and returns a bru-object. 

inlabru.lc.cohort.1 <- function(obs){
  # set-up of inlabru analysis of Poisson Lee-Carter model (no cohort effect)
  # Includes separate intercept
  # alpha x modeled as random walk and constrained to sum to zero. 
  
  nx = length(unique(obs$x))
  nt = length(unique(obs$t))
  nc = length(unique(obs$c))
  
  #   ----   Inlabru model fitting   ----
  
  A.mat = matrix(1, nrow = 1, ncol = nx)  #  helper values for constraining of beta
  e.vec = 1  #  helper values for constraining of beta
  
  # pc.prior <- list(prec = list(prior = "pc.prec", param = c(1,0.05)))
  loggamma.prior <- list(prec = list(prior = 'loggamma', param = c(1,0.00005), initial = 4))  # Default values. Sufficiently uninformative?
  
  comp = ~ -1 +
    Int(1) +
    alpha(x, model = "rw1", values=unique(obs$x), hyper = loggamma.prior, constr = TRUE) +
    phi(t, model = "linear", prec.linear = 1, mean.linear = 0) +
    beta(x.c, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = loggamma.prior) +
    kappa(t.c, model = "rw1", values = unique(obs$t), constr = TRUE, hyper = loggamma.prior) +
    gamma(c, model = "rw1", values = unique(obs$c), constr = TRUE, hyper = loggamma.prior) +
    epsilon(xt, model = "iid", hyper = loggamma.prior)
  
  formula = Y ~ Int + alpha + beta*phi + beta*kappa + gamma + epsilon
  
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

inlabru.lc.1 <- function(obs){
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
  
  comp = ~ -1 +
    Int(1) +
    alpha(x, model = "rw1", values=unique(obs$x), hyper = loggamma.prior, constr = TRUE) +
    phi(t, model = "linear", prec.linear = 1, mean.linear = 0) +
    beta(x.c, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = loggamma.prior) +
    kappa(t.c, model = "rw1", values = unique(obs$t), constr = TRUE, hyper = loggamma.prior) +
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