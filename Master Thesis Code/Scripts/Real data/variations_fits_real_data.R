# test different variations of fits to the real data
library("tidyverse")

load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data/population-germany.Rda")

load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data/lungCancer-germany.Rda")
lung.cancer <- cancer.data

load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data/stomachCancer-germany.Rda")
stomach.cancer <- cancer.data

# working directory at <Master\ Thesis\ Code>
source("Scripts/Functions/inlabru_analyses.R")

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

# save trace plots:
source("Scripts/Functions/plotters.R")

# save trace of intercept:
trace.intercept <- trace_plot(intercept_draws, chains = 6, iterations = 200000, warmup = 20000, title= "Trace plot intercept - stomach cancer, male")
save.figure(trace.intercept, name="trace_intercept", path="Scripts/Real data/Output/Figures/stomach_rw2_lc/male", pdf=F)

trace.alpha <- trace_plot_matrix(alpha_draws, chains = 6, iterations = 200000, warmup = 20000, title = "Trace plot alpha - stomach cancer, male")
save.figure(trace.alpha, name="trace_alpha", path = "Scripts/Real data/Output/Figures/stomach_rw2_lc/male", pdf=F)

trace.beta <- trace_plot_matrix(beta_draws, chains = 6, iterations = 200000, warmup = 20000, title = "Trace plot beta - stomach cancer, male")
save.figure(trace.beta, name="trace_beta", path = "Scripts/Real data/Output/Figures/stomach_rw2_lc/male", pdf=F)

trace.kappa <- trace_plot_matrix(kappa_draws, chains = 6, iterations = 200000, warmup = 20000, title = "Trace plot kappa - stomach cancer, male")
save.figure(trace.kappa, name="trace_kappa", path = "Scripts/Real data/Output/Figures/stomach_rw2_lc/male", pdf=F)

trace.tau.alpha <- trace_plot(tau_alpha_draws, chains = 6, iterations = 200000, warmup = 20000, title = "Trace plot Precision of Alpha - stomach cancer, male")
save.figure(trace.tau.alpha, name = "trace_tau_alpha", path = "Scripts/Real data/Output/Figures/stomach_rw2_lc/male", pdf=F)

trace.tau.beta <- trace_plot(tau_beta_draws, chains = 6, iterations = 200000, warmup = 20000, title = "Trace plot Precision of Beta - stomach cancer, male")
save.figure(trace.tau.beta, name = "trace_tau_beta", path = "Scripts/Real data/Output/Figures/stomach_rw2_lc/male", pdf=F)

trace.tau.kappa <- trace_plot(tau_kappa_draws, chains = 6, iterations = 200000, warmup = 20000, title = "Trace plot Precision of Kappa - stomach cancer, male")
save.figure(trace.tau.kappa, name = "trace_tau_kappa", path = "Scripts/Real data/Output/Figures/stomach_rw2_lc/male", pdf=F)


stan.marginals.stomach.m <- list(intercept_draws = intercept_draws,
                       tau_epsilon_draws = tau_epsilon_draws,
                       tau_alpha_draws = tau_alpha_draws,
                       tau_beta_draws = tau_beta_draws,
                       tau_kappa_draws = tau_kappa_draws,
                       alpha_draws = alpha_draws,
                       beta_draws = beta_draws,
                       kappa_draws = kappa_draws)#,
                       #eta_draws = eta_draws)


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
                                 kappa_draws = kappa_draws)#,
                                 #eta_draws = eta_draws)

# save trace plots:
source("Scripts/Functions/plotters.R")

# save trace of intercept:
trace.intercept <- trace_plot(intercept_draws, chains = 6, iterations = 200000, warmup = 20000, title= "Trace plot intercept - stomach cancer, female")
save.figure(trace.intercept, name="trace_intercept", path="Scripts/Real data/Output/Figures/stomach_rw2_lc/female", pdf=F)

trace.alpha <- trace_plot_matrix(alpha_draws, chains = 6, iterations = 200000, warmup = 20000, title = "Trace plot alpha - stomach cancer, female")
save.figure(trace.alpha, name="trace_alpha", path = "Scripts/Real data/Output/Figures/stomach_rw2_lc/female", pdf=F)

trace.beta <- trace_plot_matrix(beta_draws, chains = 6, iterations = 200000, warmup = 20000, title = "Trace plot beta - stomach cancer, female")
save.figure(trace.beta, name="trace_beta", path = "Scripts/Real data/Output/Figures/stomach_rw2_lc/female", pdf=F)

trace.kappa <- trace_plot_matrix(kappa_draws, chains = 6, iterations = 200000, warmup = 20000, title = "Trace plot kappa - stomach cancer, female")
save.figure(trace.kappa, name="trace_kappa", path = "Scripts/Real data/Output/Figures/stomach_rw2_lc/female", pdf=F)

trace.tau.alpha <- trace_plot(tau_alpha_draws, chains = 6, iterations = 200000, warmup = 20000, title = "Trace plot Precision of Alpha - stomach cancer, female")
save.figure(trace.tau.alpha, name = "trace_tau_alpha", path = "Scripts/Real data/Output/Figures/stomach_rw2_lc/female")

trace.tau.beta <- trace_plot(tau_beta_draws, chains = 6, iterations = 200000, warmup = 20000, title = "Trace plot Precision of Beta - stomach cancer, female")
save.figure(trace.tau.beta, name = "trace_tau_beta", path = "Scripts/Real data/Output/Figures/stomach_rw2_lc/female")

trace.tau.kappa <- trace_plot(tau_kappa_draws, chains = 6, iterations = 200000, warmup = 20000, title = "Trace plot Precision of Kappa - stomach cancer, female")
save.figure(trace.tau.kappa, name = "trace_tau_kappa", path = "Scripts/Real data/Output/Figures/stomach_rw2_lc/female")


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

trace.tau.alpha <- trace_plot(tau_alpha_draws, chains = 5, iterations = 200000, warmup = 20000, title = "Trace plot Precision of Alpha - lung cancer, male")
save.figure(trace.tau.alpha, name = "trace_tau_alpha", path = "Scripts/Real data/Output/Figures/lung_rw2_lc/male")

trace.tau.beta <- trace_plot(tau_beta_draws, chains = 5, iterations = 200000, warmup = 20000, title = "Trace plot Precision of Beta - lung cancer, male")
save.figure(trace.tau.beta, name = "trace_tau_beta", path = "Scripts/Real data/Output/Figures/lung_rw2_lc/male")

trace.tau.kappa <- trace_plot(tau_kappa_draws, chains = 5, iterations = 200000, warmup = 20000, title = "Trace plot Precision of Kappa - lung cancer, male")
save.figure(trace.tau.kappa, name = "trace_tau_kappa", path = "Scripts/Real data/Output/Figures/lung_rw2_lc/male")


stan.marginals.lung.m <- list(intercept_draws = intercept_draws,
                                 tau_epsilon_draws = tau_epsilon_draws,
                                 tau_alpha_draws = tau_alpha_draws,
                                 tau_beta_draws = tau_beta_draws,
                                 tau_kappa_draws = tau_kappa_draws,
                                 alpha_draws = alpha_draws,
                                 beta_draws = beta_draws,
                                 kappa_draws = kappa_draws)#,
                                 #eta_draws = eta_draws)


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



#   ----   Omitting ages < 45, x < 9   ----

# male stomach

stomach.cancer.male.above.45 <- stomach.cancer.male %>%
  filter(age.int >= 45)

res.stomach.lc.m.a45 <- inlabru.rw2.lc.2(stomach.cancer.male.above.45)

plots.stomach.male <- plot.inlabru.real(
  res.stomach.lc.m.a45, stomach.cancer.male.above.45, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2_lc_a45/male",
  cohort=FALSE)

plot.hypers.inlabru.real(
  res.stomach.lc.m.a45, stomach.cancer.male.above.45, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2_lc_a45/male",
  cohort=FALSE)

# lcc-model
res.stomach.cohort.m.a45 <- inlabru.rw2.cohort.2(stomach.cancer.male.above.45)

plots.stomach.male.cohort <- plot.inlabru.real(
  res.stomach.cohort.m.a45, stomach.cancer.male.above.45, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2_a45/male",
  cohort=TRUE)

plot.hypers.inlabru.real(
  res.stomach.cohort.m.a45, stomach.cancer.male.above.45, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2_a45/male",
  cohort=TRUE)

# female stomach

stomach.cancer.female.above.45 <- stomach.cancer.female %>%
  filter(age.int >= 45)

res.stomach.lc.f.a45 <- inlabru.rw2.lc.2(stomach.cancer.female.above.45)

plots.stomach.female <- plot.inlabru.real(
  res.stomach.lc.f.a45, stomach.cancer.female.above.45, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2_lc_a45/female",
  cohort=FALSE)

plot.hypers.inlabru.real(
  res.stomach.lc.f.a45, stomach.cancer.female.above.45, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2_lc_a45/female",
  cohort=FALSE)







#   ----   Prediction    ----

lung.cancer.female.until2011 <- lung.cancer.female %>% 
  mutate(Y_full = Y) %>%
  mutate(Y = replace(Y, year %in% c("2012", "2013", "2014", "2015", "2016"), NA)) %>%
  mutate(female = replace(female, year %in% c("2012", "2013", "2014", "2015", "2016"), NA)) %>%
  mutate(female.mr = replace(female.mr, year %in% c("2012", "2013", "2014", "2015", "2016"), NA)) %>%
  mutate(predict = "observed") %>% mutate(predict = replace(predict, year %in% c("2012", "2013", "2014", "2015", "2016"), "predicted"))

res.lung.f.u11.predict <- inlabru.rw2.cohort.2(lung.cancer.female.until2011, max_iter = 100)

# produce posterior samples
lung.cancer.female.until2011.samps <- generate(
  res.lung.f.u11.predict,
  data = data.frame(
    x = lung.cancer.female.until2011$x, t = lung.cancer.female.until2011$t,
    x.c = lung.cancer.female.until2011$x.c, c = lung.cancer.female.until2011$c,
    xt = lung.cancer.female.until2011$xt),
  formula = ~ Int + alpha + beta*kappa + gamma + epsilon,
  n.samples = 1000)

lung.cancer.female.until2011.lambda <- lung.cancer.female.until2011$E * exp(lung.cancer.female.until2011.samps)

lung.cancer.female.until2011.Y <- matrix(rpois(324*1000, lambda = lung.cancer.female.until2011.lambda), nrow = 324, ncol = 1000)

lung.cancer.female.until2011.Y.df <- data.frame(x = lung.cancer.female.until2011$x, t = lung.cancer.female.until2011$t,
                           xt = lung.cancer.female.until2011$xt, c = lung.cancer.female.until2011$c,
                           mean = apply(lung.cancer.female.until2011.Y, 1, mean),
                           X0.975 = apply(lung.cancer.female.until2011.Y, 1, quantile, 0.975),
                           X0.025 = apply(lung.cancer.female.until2011.Y, 1, quantile, 0.025))


source("Scripts/Real data/plot_real_data.R")

plots.lung <- plot.inlabru.real.predicted(
  res.lung.f.u11.predict, lung.cancer.female.until2011, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_rw2_predict/female")

plot.hypers.inlabru.real(
  res.lung.f.u11.predict, lung.cancer.female, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_rw2_predict/female")

plot.counts.inlabru.real.predicted(
  lung.cancer.female.until2011.Y.df, lung.cancer.female.until2011, save = TRUE,
  path.to.storage = "Scripts/Real data/Output/Figures/lung_rw2_predict/female"
)




