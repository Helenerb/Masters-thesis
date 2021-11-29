# test different variations of fits to the real data
library("tidyverse")

load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data/population-germany.Rda")

load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data/lungCancer-germany.Rda")
lung.cancer <- cancer.data

load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data/stomachCancer-germany.Rda")
stomach.cancer <- cancer.data

# working directory at <Master\ Thesis\ Code>
source("Scripts/Synthetic data/Inlabru analyses/inlabru_analyses.R")

# run analyisis for male and female cancer mortality separately:

#   ----   reformat data   ----

# calculate male and female mortality:

lung.cancer <- lung.cancer %>%
  mutate(female.mr = female/female.t, male.mr = male/male.t)

stomach.cancer <- stomach.cancer %>%
  mutate(female.mr = female/female.t, male.mr = male/male.t) 
# separate dfs for male and female mortality:

lung.cancer.female <- lung.cancer %>%
  select(x, x.c, t, xt, cohort, c, age, age.int, year, birth.year, female.t, female, female.mr) %>%
  mutate(E = female.t, Y = female, `mortality rate` = female.mr)

stomach.cancer.female <- stomach.cancer %>%
  select(x, x.c, t, xt, cohort, c, age, age.int, year, birth.year, female.t, female, female.mr) %>%
  mutate(E = female.t, Y = female, `mortality rate` = female.mr)

lung.cancer.male <- lung.cancer %>%
  select(x, x.c, t, xt, cohort, c, age, age.int, year, birth.year, male.t, male, male.mr) %>%
  mutate(E = male.t, Y = male, `mortality rate` = male.mr)

stomach.cancer.male <- stomach.cancer %>%
  select(x, x.c, t, xt, cohort, c, age, age.int, year, birth.year, male.t, male, male.mr) %>%
  mutate(E = male.t, Y = male, `mortality rate` = male.mr)

#   ----   perform inlabru analyses   ----
source("Scripts/Synthetic data/Inlabru analyses/inlabru_analyses.R")

# lc-models 

res.lung.lc.f <- inlabru.rw2.lc.2(lung.cancer.female)

res.stomach.lc.f <- inlabru.rw2.lc.2(stomach.cancer.female, max_iter=100)

res.lung.lc.m <- inlabru.rw2.lc.2(lung.cancer.male)

res.stomach.lc.m <- inlabru.rw2.lc.2(stomach.cancer.male)

# lcc-models

res.lung.cohort.f <- inlabru.rw2.cohort.2(lung.cancer.female)

res.stomach.cohort.f <- inlabru.rw2.cohort.2(stomach.cancer.female, max_iter = 100)

res.lung.cohort.m <- inlabru.rw2.cohort.2(lung.cancer.male, max_iter=100)

res.stomach.cohort.m <- inlabru.rw2.cohort.2(stomach.cancer.male, max_iter=100)


#   ----   plot and save results   ----

source("Scripts/Real data/plot_real_data.R")

# female lung
plots.lung.female <- plot.inlabru.real(
  res.lung.lc.f, lung.cancer.female, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_rw2_lc/female",
  cohort=FALSE)

plot.hypers.inlabru.real(
  res.lung.lc.f, lung.cancer.female, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_rw2_lc/female",
  cohort=FALSE)

# female stomach
plots.stomach.female <- plot.inlabru.real(
  res.stomach.lc.f, stomach.cancer.female, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2_lc/female",
  cohort=FALSE)

plot.hypers.inlabru.real(
  res.stomach.lc.f, stomach.cancer.female, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2_lc/female",
  cohort=FALSE)

# male lung
plots.lung.male <- plot.inlabru.real(
  res.lung.lc.m, lung.cancer.male, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_rw2_lc/male",
  cohort=FALSE)

plot.hypers.inlabru.real(
  res.lung.lc.m, lung.cancer.male, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_rw2_lc/male",
  cohort=FALSE)

# male stomach
plots.stomach.male <- plot.inlabru.real(
  res.stomach.lc.m, stomach.cancer.male, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2_lc/male",
  cohort=FALSE)

plot.hypers.inlabru.real(
  res.stomach.lc.m, stomach.cancer.male, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2_lc/male",
  cohort=FALSE)

# plot lcc:

# female lung
plots.lung.female <- plot.inlabru.real(
  res.lung.cohort.f, lung.cancer.female, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_rw2/female",
  cohort=TRUE)

plot.hypers.inlabru.real(
  res.lung.cohort.f, lung.cancer.female, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_rw2/female",
  cohort=TRUE)

# female stomach
plots.stomach.female <- plot.inlabru.real(
  res.stomach.cohort.f, stomach.cancer.female, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2/female",
  cohort=TRUE)

plot.hypers.inlabru.real(
  res.stomach.cohort.f, stomach.cancer.female, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2/female",
  cohort=T)

# male lung
plots.lung.male <- plot.inlabru.real(
  res.lung.cohort.m, lung.cancer.male, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_rw2/male",
  cohort=T)

plot.hypers.inlabru.real(
  res.lung.cohort.m, lung.cancer.male, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_rw2/male",
  cohort=T)

# male stomach
plots.stomach.male <- plot.inlabru.real(
  res.stomach.cohort.m, stomach.cancer.male, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2/male",
  cohort=T)

plot.hypers.inlabru.real(
  res.stomach.cohort.m, stomach.cancer.male, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2/male",
  cohort=T)

#   ----   Comparison plots   ----

# stomach male:
load("Scripts/Real data/Stan analyses/stomach_rw2_lc_male/stan_results/stan_stomach_rw2_lc_male.Rda") 
stan_stomach_rw2_lc_m <- stan_lc_df

path.to.stan.results = "Scripts/Real\ data/Stan analyses/stomach_rw2_lc_male/stan_results"

load(file=file.path(path.to.stan.results, "draws_intercept.RData"))
load(file=file.path(path.to.stan.results, "draws_tau_epsilon.RData"))
load(file.path(path.to.stan.results, "draws_tau_alpha.RData"))
load(file.path(path.to.stan.results, "draws_tau_beta.RData"))
load(file.path(path.to.stan.results, "draws_tau_kappa.RData"))
load(file.path(path.to.stan.results, "draws_alpha.RData"))
load(file.path(path.to.stan.results, "draws_beta.RData"))
load(file.path(path.to.stan.results, "draws_kappa.RData"))
load(file.path(path.to.stan.results, "draws_eta.RData"))

stan.marginals.stomach.m <- list(intercept_draws = intercept_draws,
                       tau_epsilon_draws = tau_epsilon_draws,
                       tau_alpha_draws = tau_alpha_draws,
                       tau_beta_draws = tau_beta_draws,
                       tau_kappa_draws = tau_kappa_draws,
                       alpha_draws = alpha_draws,
                       beta_draws = beta_draws,
                       kappa_draws = kappa_draws,
                       eta_draws = eta_draws)


plots.compared.stomach.male <- plot.comparison.real(
  res.stomach.lc.m, res.stan=stan_stomach_rw2_lc_m,
  stan.marginals=stan.marginals.stomach.m, cancer.data=stomach.cancer,
  path.to.storage="Scripts/Real data/Output/Figures/stomach_rw2_lc/male",
  cohort=FALSE, save=TRUE)

# stomach female:
load("Scripts/Real data/Stan analyses/stomach_rw2_lc_female/stan_results/stan_stomach_rw2_lc_female.Rda") 
stan_stomach_rw2_lc_f <- stan_lc_df

path.to.stan.results = "Scripts/Real\ data/Stan analyses/stomach_rw2_lc_female/stan_results"

load(file=file.path(path.to.stan.results, "draws_intercept.RData"))
load(file=file.path(path.to.stan.results, "draws_tau_epsilon.RData"))
load(file.path(path.to.stan.results, "draws_tau_alpha.RData"))
load(file.path(path.to.stan.results, "draws_tau_beta.RData"))
load(file.path(path.to.stan.results, "draws_tau_kappa.RData"))
load(file.path(path.to.stan.results, "draws_alpha.RData"))
load(file.path(path.to.stan.results, "draws_beta.RData"))
load(file.path(path.to.stan.results, "draws_kappa.RData"))
load(file.path(path.to.stan.results, "draws_eta.RData"))

stan.marginals.stomach.f <- list(intercept_draws = intercept_draws,
                                 tau_epsilon_draws = tau_epsilon_draws,
                                 tau_alpha_draws = tau_alpha_draws,
                                 tau_beta_draws = tau_beta_draws,
                                 tau_kappa_draws = tau_kappa_draws,
                                 alpha_draws = alpha_draws,
                                 beta_draws = beta_draws,
                                 kappa_draws = kappa_draws,
                                 eta_draws = eta_draws)


plots.compared.stomach.female <- plot.comparison.real(
  res.stomach.lc.f, res.stan=stan_stomach_rw2_lc_f,
  stan.marginals=stan.marginals.stomach.f, cancer.data=stomach.cancer,
  path.to.storage="Scripts/Real data/Output/Figures/stomach_rw2_lc/female",
  cohort=FALSE, save=TRUE)

# lung male:
load("Scripts/Real data/Stan analyses/lung_rw2_lc_male/stan_results/stan_lung_rw2_lc_male.Rda") 
stan_lung_rw2_lc_m <- stan_lc_df

path.to.stan.results = "Scripts/Real\ data/Stan analyses/lung_rw2_lc_male/stan_results"

load(file=file.path(path.to.stan.results, "draws_intercept.RData"))
load(file=file.path(path.to.stan.results, "draws_tau_epsilon.RData"))
load(file.path(path.to.stan.results, "draws_tau_alpha.RData"))
load(file.path(path.to.stan.results, "draws_tau_beta.RData"))
load(file.path(path.to.stan.results, "draws_tau_kappa.RData"))
load(file.path(path.to.stan.results, "draws_alpha.RData"))
load(file.path(path.to.stan.results, "draws_beta.RData"))
load(file.path(path.to.stan.results, "draws_kappa.RData"))
load(file.path(path.to.stan.results, "draws_eta.RData"))

stan.marginals.lung.m <- list(intercept_draws = intercept_draws,
                                 tau_epsilon_draws = tau_epsilon_draws,
                                 tau_alpha_draws = tau_alpha_draws,
                                 tau_beta_draws = tau_beta_draws,
                                 tau_kappa_draws = tau_kappa_draws,
                                 alpha_draws = alpha_draws,
                                 beta_draws = beta_draws,
                                 kappa_draws = kappa_draws,
                                 eta_draws = eta_draws)


plots.compared.lung.male <- plot.comparison.real(
  res.lung.lc.m, res.stan=stan_lung_rw2_lc_m,
  stan.marginals=stan.marginals.lung.m, cancer.data=lung.cancer,
  path.to.storage="Scripts/Real data/Output/Figures/lung_rw2_lc/male",
  cohort=FALSE, save=TRUE)

# lung female:
load("Scripts/Real data/Stan analyses/lung_rw2_lc_female/stan_results/stan_lung_rw2_lc_female.Rda") 
stan_lung_rw2_lc_f <- stan_lc_df

path.to.stan.results = "Scripts/Real\ data/Stan analyses/lung_rw2_lc_female/stan_results"

load(file=file.path(path.to.stan.results, "draws_intercept.RData"))
load(file=file.path(path.to.stan.results, "draws_tau_epsilon.RData"))
load(file.path(path.to.stan.results, "draws_tau_alpha.RData"))
load(file.path(path.to.stan.results, "draws_tau_beta.RData"))
load(file.path(path.to.stan.results, "draws_tau_kappa.RData"))
load(file.path(path.to.stan.results, "draws_alpha.RData"))
load(file.path(path.to.stan.results, "draws_beta.RData"))
load(file.path(path.to.stan.results, "draws_kappa.RData"))
load(file.path(path.to.stan.results, "draws_eta.RData"))

stan.marginals.lung.f <- list(intercept_draws = intercept_draws,
                              tau_epsilon_draws = tau_epsilon_draws,
                              tau_alpha_draws = tau_alpha_draws,
                              tau_beta_draws = tau_beta_draws,
                              tau_kappa_draws = tau_kappa_draws,
                              alpha_draws = alpha_draws,
                              beta_draws = beta_draws,
                              kappa_draws = kappa_draws,
                              eta_draws = eta_draws)


plots.compared.lung.female <- plot.comparison.real(
  res.lung.lc.f, res.stan=stan_lung_rw2_lc_f,
  stan.marginals=stan.marginals.lung.f, cancer.data=lung.cancer,
  path.to.storage="Scripts/Real data/Output/Figures/lung_rw2_lc/female",
  cohort=FALSE, save=TRUE)




