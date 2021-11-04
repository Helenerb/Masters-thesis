# script for creating "synthetic data" from real data
# we use the male stomach cancer data, when fitted to the 
# rw2 lc-model with inlabru

synthetic.male.stomach.lc <- function(){

  # save from working environment:
  #save(res.stomach.lc.m, file=file.path("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Real data/Output/Data/stomach_rw2_lc_male", "inlabru_res.RData"))
  
  # working directory at .../Synthetic data/Stan analyses
  load(file.path("../Real data/Output/Data/stomach_rw2_lc_male", "inlabru_res.RData"))
  
  load("../../Data/population-germany.Rda")
  
  load("../../Data/stomachCancer-germany.Rda")
  stomach.cancer <- cancer.data
  
  stomach.cancer.male <- stomach.cancer %>%
    mutate(female.mr = female/female.t, male.mr = male/male.t) %>%
    select(x, x.c, t, xt, cohort, c, age, age.int, year, birth.year, male.t, male, male.mr) %>%
    mutate(E = male.t, Y = male, `mortality rate` = male.mr)
  
  res.inlabru <- res.stomach.lc.m
    
  nx <- max(res.inlabru$summary.random$alpha$ID) + 1
  nt <- max(res.inlabru$summary.random$kappa$ID) + 1
  
  alpha.true = res.inlabru$summary.random$alpha$mean
  alpha.true = alpha.true - mean(alpha.true)
  
  beta.true = res.inlabru$summary.random$beta$mean
  beta.true = beta.true - mean(beta.true) + 1/nx
  
  kappa.true =  res.inlabru$summary.random$kappa$mean
  kappa.true = kappa.true - mean(kappa.true)
  
  epsilon.true = res.inlabru$summary.random$epsilon$mean
  
  intercept.true = res.inlabru$summary.fixed$mean[1]
  
  obs <- data.frame(x = stomach.cancer.male$x, t = stomach.cancer.male$t,
                    xt = stomach.cancer.male$xt, E = stomach.cancer.male$E) %>%
    mutate(alpha = alpha.true[x + 1]) %>%
    mutate(intercept = intercept.true) %>%
    mutate(beta = beta.true[x + 1]) %>%
    mutate(kappa = kappa.true[t + 1]) %>%
    mutate(epsilon = epsilon.true[xt + 1]) %>%
    mutate(eta = intercept + alpha + beta*kappa + epsilon) %>%
    mutate(Y = rpois(length(x), E*exp(eta))) %>%
    mutate(x.c = x, t.c = t)
  
  underlying.effects <- list(
    obs = obs,
    alpha.true = alpha.true,
    age.intercept.true= intercept.true,
    beta.true = beta.true,
    kappa.true = kappa.true,
    tau.beta.true = NA,
    tau.kappa.true = NA,
    tau.epsilon.true= NA,
    config_name = "male_stomach_rw2_lc",
    nx = nx, 
    nt = nt
  )
  
  return(underlying.effects)

}

