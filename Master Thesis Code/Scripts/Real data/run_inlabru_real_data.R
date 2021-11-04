# script for running inlabru analyses with the real data. 

load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data/population-germany.Rda")

load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data/lungCancer-germany.Rda")
lung.cancer <- cancer.data

load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data/stomachCancer-germany.Rda")
stomach.cancer <- cancer.data

# working directory at <Master\ Thesis\ Code>
source("Scripts/Synthetic data/Inlabru analyses/inlabru_analyses.R")

res.stomach <- inlabru.rw2.cohort.2(stomach.cancer)

# the lc-version:
res.stomach.lc <- inlabru.rw2.lc.2(stomach.cancer)

source("Scripts/Real data/plot_real_data.R")

plots.stomach <- plot.inlabru.real(
  res.stomach, stomach.cancer, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2",
  cohort=TRUE)

plot.hypers.inlabru.real(
  res.stomach, stomach.cancer, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2")

# run for lung cancer
res.lung <- inlabru.rw2.cohort.2(lung.cancer, max_iter = 100)

# lc-version:
res.lung.lc <- inlabru.rw2.lc.2(lung.cancer, max_iter = 100)

plots.lung <- plot.inlabru.real(
  res.lung, lung.cancer, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_rw2",
  cohort=TRUE)

plot.hypers.inlabru.real(
  res.lung, lung.cancer, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_rw2")

# run with extraconstraints on gamma, stomach
res.stomach.extraconstr <- inlabru.rw2.cohort.2.gamma.extraconstr(stomach.cancer, real_data = TRUE)

plot.hypers.inlabru.real(
  res.stomach.extraconstr, stomach.cancer, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2_extraconstr")

plots.stomach.extraconstr <- plot.inlabru.real(
  res.stomach.extraconstr, stomach.cancer, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2_extraconstr")

# run with extraconstraints on gamma, lung
res.lung.extraconstr <- inlabru.rw2.cohort.2.gamma.extraconstr(lung.cancer, real_data = TRUE, max_iter=100)

plot.hypers.inlabru.real(
  res.lung.extraconstr, lung.cancer, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_rw2_extraconstr")

plots.lung.extraconstr <- plot.inlabru.real(
  res.lung.extraconstr, lung.cancer, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_rw2_extraconstr")

save(res.lung.extraconstr, file="Scripts/Real data/Output/Data/lung_rw2_extraconstr/res_inlabru.RData")

#   ----    Run predictions   ----   
stomach.cancer.until2007 <- stomach.cancer %>% 
  mutate(Y_full = Y) %>%
  mutate(total = replace(total, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), NA)) %>%
  mutate(Y = replace(Y, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), NA)) %>%
  mutate(male = replace(male, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), NA)) %>%
  mutate(female = replace(female, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), NA)) %>%
  mutate(predict = "observed") %>% mutate(predict = replace(predict, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), "predicted"))

res.stomach.predict <- inlabru.rw2.cohort.2(stomach.cancer.until2007, max_iter = 100)

plots.stomach <- plot.inlabru.real.predicted(
  res.stomach.predict, stomach.cancer.until2007, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2_predict")

plot.hypers.inlabru.real(
  res.stomach.predict, stomach.cancer, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2_predict")

# prediction with ar1c:

res.stomach.predict.ar1c <- inlabru.ar1c.cohort.2(stomach.cancer.until2007)

plots.stomach <- plot.inlabru.real.predicted(
  res.stomach.predict.ar1c, stomach.cancer.until2007, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_ar1c_predict")

plot.hypers.inlabru.real(
  res.stomach.predict.ar1c, stomach.cancer, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_ar1c_predict")

# same for lung cancer:

lung.cancer.until2007 <- lung.cancer %>% 
  mutate(Y_full = Y) %>%
  mutate(total = replace(total, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), NA)) %>%
  mutate(Y = replace(Y, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), NA)) %>%
  mutate(male = replace(male, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), NA)) %>%
  mutate(female = replace(female, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), NA)) %>%
  mutate(predict = "observed") %>% mutate(predict = replace(predict, year %in% c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016"), "predicted"))

res.lung.predict <- inlabru.rw2.cohort.2(lung.cancer.until2007, max_iter = 100)

plots.lung <- plot.inlabru.real.predicted(
  res.lung.predict, lung.cancer.until2007, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_rw2_predict")

plot.hypers.inlabru.real(
  res.lung.predict, lung.cancer, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_rw2_predict")

# prediction with ar1c:

res.lung.predict.ar1c <- inlabru.ar1c.cohort.2(lung.cancer.until2007)

plots.lung <- plot.inlabru.real.predicted(
  res.lung.predict.ar1c, lung.cancer.until2007, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_ar1c_predict")

plot.hypers.inlabru.real(
  res.lung.predict.ar1c, lung.cancer, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_ar1c_predict")


#   ----   Compare results of STAN and inlabru:   ----

# load("Scripts/Real data/Stan analyses/lung_rw2_lc/stan_results/stan_lung_rw2_lc.Rda")  # TODO: Check that this is the correct path!
# stan_lung_rw2_lc <- stan_lc_df

load("Scripts/Real data/Stan analyses/stomach_rw2_lc/stan_results/stan_stomach_rw2_lc.Rda")  # TODO: Check that this is the correct path!
stan_stomach_rw2_lc <- stan_lc_df

path.to.stan.results = "Scripts/Real\ data/Stan analyses/stomach_rw2_lc/stan_results"

load(file=file.path(path.to.stan.results, "draws_intercept.RData"))
load(file=file.path(path.to.stan.results, "draws_tau_epsilon.RData"))
load(file.path(path.to.stan.results, "draws_tau_alpha.RData"))
load(file.path(path.to.stan.results, "draws_tau_beta.RData"))
load(file.path(path.to.stan.results, "draws_tau_kappa.RData"))
load(file.path(path.to.stan.results, "draws_alpha.RData"))
load(file.path(path.to.stan.results, "draws_beta.RData"))
load(file.path(path.to.stan.results, "draws_kappa.RData"))
load(file.path(path.to.stan.results, "draws_eta.RData"))

stan.marginals <- list(intercept_draws = intercept_draws,
                       tau_epsilon_draws = tau_epsilon_draws,
                       tau_alpha_draws = tau_alpha_draws,
                       tau_beta_draws = tau_beta_draws,
                       tau_kappa_draws = tau_kappa_draws,
                       alpha_draws = alpha_draws,
                       beta_draws = beta_draws,
                       kappa_draws = kappa_draws)#,
                      # eta_draws = eta_draws)


# plots.compared.lung <- plot.comparison.real(
#   res.lung.lc, res.stan=stan_lung_rw2_lc, 
#   stan.marginals=stan.marginals, cancer.data=lung.cancer, 
#   path.to.storage="Scripts/Real data/Output/Figures/lung_rw2_lc",
#   cohort=FALSE, save=TRUE)

plots.compared.stomach <- plot.comparison.real(
  res.stomach.lc, res.stan=stan_stomach_rw2_lc,
  stan.marginals=stan.marginals, cancer.data=stomach.cancer,
  path.to.storage="Scripts/Real data/Output/Figures/stomach_rw2_lc",
  cohort=FALSE, save=TRUE)






