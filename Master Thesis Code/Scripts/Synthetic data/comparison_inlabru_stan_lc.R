library(patchwork)

setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data")
source("configurations_synthetic_data.R")

# configuration without cohort effect
#underlying.effects.lc <- configuration.v5()
#underlying.effects.lc <- configuration.v9()  #  config with coarser grid
#underlying.effects.lc <- configuration.v10()
#underlying.effects.lc <- configuration.v10.1()# config with coarser grid and larger variance in beta
underlying.effects.lc <- configuration.v10.2()
#underlying.effects.lc <- configuration.v11()  # config with coarser grid, larger variance in beta and steeper phi (compared to alpha)
#underlying.effects.lc <- configuration.v11.1()
#underlying.effects.lc <- configuration.v12()  #  config with coarser grid, smaller variance in beta, steeper phi (compared to alpha)
#underlying.effects.lc <- configuration.v13()  # config with coarser grid, even smaller variance in beta, a bit less steep phi, higher variance in kappa
#underlying.effects.lc <- configuration.v14()
#underlying.effects.lc <- configuration.v15()  # half the grid of the original v5

figures.folder = "/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Output/Figures"

#storage_path = file.path(figures.folder, "v10_2")
storage_path = file.path(figures.folder, "v10d")
# storage_path = file.path(figures.folder, "v10dh")

obs.lc <- underlying.effects.lc$obs

source("Inlabru\ analyses/inlabru_analyses.R")
runtime.inlabru <- system.time({res.inlabru.lc.1 <- inlabru.lc.kappa_high_variance_prior(obs.lc)})

source("plot_inlabru_vs_underlying.R")

plots.lc <- plot.inlabru.vs.underlying.lc(
  res.inlabru.lc.1, 
  underlying.effects.lc,
  path.to.storage = storage_path,
  save=TRUE,
  phi.plus.kappa.func = phi.plus.kappa.v17)

plots.lc$p.alpha
plots.lc$p.beta
plots.lc$p.phi
plots.lc$p.intercept
plots.lc$p.kappa
plots.lc$p.eta
plots.lc$p.eta.2
plots.lc$p.eta.t
plots.lc$p.eta.x
plots.lc$p.gamma


print("Runtime for inlabru: ")
print(runtime.inlabru)

#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v10/stan_results/stan_v10.Rda")
#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v10d/stan_results/stan_v10d.Rda")
load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Stan analyses/v10dh/stan_results/stan_v10dh.Rda")

source("plot_stan_vs_underlying.R")

#undrifted version of stan program:
# stan.results <- produce.stan.plots(stan_df=stan_lc_df,
#                    underlying.effects=underlying.effects.lc,
#                    plot.func=plot.stan.vs.underlying.lc.undrifted,
#                    save.func=save.stan.plots.lc,
#                    path.to.storage=storage_path,
#                    summaries.func=produce.summaries.stan.lc)

# drifted version, soft constraints
# stan.res <- produce.stan.plots(stan_df=stan_lc_df,
#                                underlying.effects=underlying.effects.lc,
#                                plot.func=plot.stan.vs.underlying.lc.drifted,
#                                save.func=save.stan.plots.lc,
#                                path.to.storage=storage_path,
#                                summaries.func=produce.summaries.stan.lc)

# drifted version of stan program, hard constraints
stan.res <- produce.stan.plots(stan_df=stan_lc_df,
                               underlying.effects=underlying.effects.lc,
                               plot.func=plot.stan.vs.underlying.lc.drifted,
                               save.func=save.stan.plots.lc,
                               path.to.storage=storage_path,
                               summaries.func=produce.summaries.stan.hard.lc)


# samps = inla.posterior.sample(res.inlabru.lc.1, n = 1000)
# 
# phi.plus.kappa <- function(t_max){
#   t = 0:t_max
#   res = kappa + phi*t
#   return(res)
# }
# 
# phi.plus.kappa.99 <- function(){
#   t = 0:99
#   res = kappa + phi*t
#   return(res)
# }

# posterior <- inla.posterior.sample.eval(fun = phi.plus.kappa.99, samples=samps)
# 
# posterior.df <- data.frame(t = 1:100,
#                            mean = apply(posterior, 1, mean),
#                            q1 = apply(posterior, 1, quantile, 0.025),
#                            q2 = apply(posterior, 1, quantile, 0.975)) %>%
#   mutate(kappa = underlying.effects.lc$kappa.true[t]) %>%
#   mutate(phi.t = underlying.effects.lc$phi.true*(t-1)) %>%
#   mutate(kappa.phi = kappa + phi.t)
# 
# gg.posterior <- ggplot(data = posterior.df) +
#   geom_ribbon(aes(x = t, ymin = q1, ymax = q2, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
#   geom_point(aes(x = t, y = mean, color = "Inlabru", fill = "Inlabru")) +
#   geom_point(aes(x = t, y = kappa.phi, color = "True", fill = "True")) +
#   scale_color_manual(name = " ", values = palette.basis) + 
#   scale_fill_manual(name = " ", values = palette.basis) + 
#   labs(title = "Phi*t + kappa", x = "t", y = "")
# 
# gg.posterior

# phi.kappa.samp <- inla.posterior.sample.eval(samples=samps, c("phi", "kappa"))
# 
# samp.matrix <- matrix(rep(1:20, 100), byrow = FALSE, nrow= 20)
# phi.matrix <- matrix(rep(phi.kappa.samp[1,], 20), nrow = 20, byrow=TRUE)
# samp.matrix <- samp.matrix * phi.matrix + phi.kappa.samp[-1,]
# 
# data.frame(t = 1:20, mean = apply(samp.matrix, 1, mean),
#            q1 = apply(samp.matrix, 1, quantile, 0.025),
#            q2 = apply(samp.matrix, 1, quantile, 0.975)) %>%
#   ggplot() + 
#   geom_line(aes(x = t, y= mean)) + 
#   geom_line(aes(x = t, y = q1), alpha = 0.5) +
#   geom_line(aes(x = t, y = q2), alpha = 0.5)
# 
# samples[1]

# ggplot(data = posterior.phi) + geom_point(aes(x = t, y = mean)) + geom_line(aes(x = t, y =q0.025)) + geom_line(aes(x = t, y = q0.975))
# 

#res.inlabru.lc.1 <- inlabru.lc.pc_priors(obs.lc)
#res.inlabru.lc.1 <- inlabru.lc.kappa_pc_prior(obs.lc)
#res.inlabru.lc.1 <- inlabru.lc.kappa_high_variance_prior(obs.lc)
# res.inlabru.lc.1 <- inlabru.lc.alpha.iid(obs.lc)

#res.inlabru.lc.1 <- inlabru.lc.1(obs.lc)

# #   ----   plotting results of inlabru fit   ----
# source("plot_inlabru_vs_underlying.R")
# 
# # plotting results from run with cohort effects:
# plots <- plot.inlabru.vs.underlying.v5(res.inlabru.lc.1, underlying.effects.lc)
# plots$p.alpha
# plots$p.beta
# plots$p.phi
# #plots$p.intercept
# plots$p.kappa
# plots$p.eta
# plots$p.eta.2
# plots$p.eta.t
# plots$p.eta.x


save.comparison.plots(plots$p.alpha, 'alpha_inlabru', path.to.folder)
save.comparison.plots(plots$p.beta, 'beta_inlabru', path.to.folder)
save.comparison.plots(plots$p.phi, 'phi_inlabru', path.to.folder)
save.comparison.plots(plots$p.intercept, 'intercept_inlabru', path.to.folder)
save.comparison.plots(plots$p.kappa, 'kappa_inlabru', path.to.folder)
save.comparison.plots(plots$p.eta, 'eta_inlabru', path.to.folder)
save.comparison.plots(plots$p.eta.2, 'eta_2_inlabru', path.to.folder)
save.comparison.plots(plots$p.eta.t, 'eta_t_inlabru', path.to.folder)
save.comparison.plots(plots$p.eta.x, 'eta_x_inlabru', path.to.folder)

#   ----   Perform similar analysis in STAN   ----   
library("rstan")

input_stan.lc <- list(
  X=length(unique(obs.lc$x)),
  T=length(unique(obs.lc$t)),
  x=(obs.lc$x + 1),
  t=(obs.lc$t + 1),
  ts = obs.lc$t,
  E = obs.lc$E,
  Y = obs.lc$Y,
  nx = length(unique(obs.lc$x)),
  nt = length(unique(obs.lc$t))
)

# save workspace image 
#save.image("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Workspaces/stan_analysis_lc.RData")

#save.image("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Workspaces/stan_analysis_lc_v9.RData")
#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Workspaces/stan_analysis_lc_v9.RData")  # configuration v9 uninformative priors

#save.image("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Workspaces/stan_analysis_lc_v9_informative_priors.RData")

#save.image("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Workspaces/stan_analysis_lc_v15.RData")
#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Workspaces/stan_analysis_lc_v15.RData")

#save.image("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Workspaces/stan_analysis_lc_v10.RData")

#save.image("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Workspaces/stan_analysis_lc_v5.RData")

#save.image("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Workspaces/stan_analysis_lc_v15_iid.RData")
#save.image("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Workspaces/stan_analysis_lc_v10_iid.RData")

#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Workspaces/stan_analysis_lc_v3.RData")
# load workspace image
#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Workspaces/stan_analysis_lc.RData")

# information variables for workspace:
workspace_information <- list(
  stan.file = "stan_analysis_lc_v4.stan",
  chains = "4",
  warmup = "1000",
  iter = "10000",
  configuration = "v9",
  type_of_prior = "loggamma, less informative on kappa"
)

# fit <- stan(
#   file="stan_analysis_lc_v3.stan",
#   data = input_stan.lc,
#   chains=4,
#   warmup = 2000,
#   iter = 15000,
#   refresh = 100,
#   seed=123
# )

# with more informative prior on kappa:

fit <- stan(
  file="stan_analysis_lc_v4.stan",
  data = input_stan.lc,
  chains=4,
  warmup = 1000,
  iter = 10000,
  refresh = 1000,
  seed=123
)

pairs(fit, pars=c("tau_alpha", "phi", "intercept", "eta[1]"))

#   ----   plot stan results   ---- 

save.comparison.plots(traceplot(fit, pars=c("eta[1]", "eta[2]", "eta[3]", "eta[4]", "eta[5]",
                                            "eta[6]", "eta[7]", "eta[8]", "eta[9]", "eta[10]" )),
                      'trace_eta_alpha_iid', path.to.folder)

save.comparison.plots(traceplot(fit, pars=c("kappa[1]", "kappa[2]", "kappa[3]", "kappa[4]",
                                            "kappa[5]", "kappa[6]", "kappa[7]", "kappa[8]",
                                            "kappa[9]", "kappa[10]" )),
                      'trace_kappa_alpha_iid', path.to.folder)

save.comparison.plots(traceplot(fit, pars=c("beta[1]", "beta[2]", "beta[3]", "beta[4]",
                                            "beta[5]", "beta[6]", "beta[7]", "beta[8]",
                                            "beta[9]", "beta[10]" )),
                      'trace_beta_alpha_iid', path.to.folder)


save.comparison.plots(traceplot(fit, pars=c("alpha[1]", "alpha[3]", "alpha[11]", "alpha[13]",
                                            "alpha[5]", "alpha[15]", "alpha[7]", "alpha[8]",
                                            "alpha[12]", "alpha[10]" )),
                      'trace_alpha_alpha_iid', path.to.folder)

save.comparison.plots(traceplot(fit, pars=c("alpha_raw[1]", "alpha_raw[3]", "alpha_raw[11]", "alpha_raw[13]",
                                            "alpha_raw[5]", "alpha_raw[15]", "alpha_raw[7]", "alpha_raw[8]",
                                            "alpha_raw[12]", "alpha_raw[10]" )),
                      'trace_alpha_raw_alpha_iid', path.to.folder)

save.comparison.plots(traceplot(fit, pars=c("beta_raw[1]", "beta_raw[3]", "beta_raw[11]", "beta_raw[13]",
                                            "beta_raw[5]", "beta_raw[15]", "beta_raw[7]", "beta_raw[8]",
                                            "beta_raw[12]", "beta_raw[10]" )),
                      'trace_beta_raw_alpha_iid', path.to.folder)

save.comparison.plots(traceplot(fit),
                      'general_trace_alpha_iid', path.to.folder)



stan_plot(fit, pars = c("intercept", "phi"))
stan_plot(fit, pars = c("tau_alpha", "tau_beta", "tau_kappa", "tau_epsilon"))

# autocorrelation plot:
library("bayesplot")
posterior_stan_lc <- as.array(fit)
acf_plot <- mcmc_acf(posterior_stan_lc, pars = "eta[1]", lags = 3)
acf_plot

print(fit)
#plot(fit)
fit_summary.stan.lc <- summary(fit)
print(fit_summary.stan.lc$summary)
stan_lc_df <- as.data.frame(fit_summary.stan.lc$summary)

# look at alpha
summary_alpha <- stan_lc_df %>%
  rownames_to_column("parameter") %>%
  filter(grepl('alpha', parameter)) %>%
  filter(!grepl('tau_alpha', parameter)) %>%
  filter(!grepl('alpha_raw', parameter)) %>%
  mutate(index = parse_number(parameter)) %>%
  mutate(true_alpha = obs.lc$alpha[index])

plot_alpa <- ggplot(data=summary_alpha) +
  geom_point(aes(x=index, y=mean, color="estimated")) + 
  geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
  geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha=0.5) +
  geom_point(aes(x=index, y=true_alpha, color="true value")) +
  ggtitle("Alpha")
plot_alpa  

# look at beta
summary_beta <- stan_lc_df %>%
  rownames_to_column("parameter") %>%
  filter(grepl('beta', parameter)) %>%
  filter(!grepl('tau_beta', parameter)) %>%
  filter(!grepl('beta_raw', parameter)) %>%
  mutate(index = parse_number(parameter)) %>%
  mutate(true_beta = underlying.effects.lc$beta.true[index])

summary_beta_raw <- stan_lc_df %>%
  rownames_to_column("parameter") %>%
  filter(grepl('beta_raw', parameter)) %>%
  mutate(index = parse_number(parameter)) %>%
  mutate(true_beta = underlying.effects.lc$beta.true[index])

plot_beta <- ggplot(data=summary_beta) +
  geom_point(aes(x=index, y=mean, color="estimated")) + 
  geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
  geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha = 0.5) +
  geom_point(aes(x=index, y=true_beta, color="true value")) +
  ggtitle("Beta -STAN")
plot_beta

# look at kappa
summary_kappa <- stan_lc_df %>%
  rownames_to_column("parameter") %>%
  filter(grepl('kappa', parameter)) %>%
  filter(!grepl('tau_kappa', parameter)) %>%
  filter(!grepl('kappa_raw', parameter)) %>%
  mutate(index = parse_number(parameter)) %>%
  mutate(true_kappa = underlying.effects.lc$kappa.true[index]) %>%
  mutate(kappa_drifted = true_kappa + underlying.effects.lc$phi.true*(index-1))

plot_kappa <- ggplot(data=summary_kappa) +
  geom_point(aes(x=index, y=mean, color="estimated")) + 
  geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
  geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha=0.5) +
  geom_point(aes(x=index, y=true_kappa, color="true value")) +
  ggtitle("Kappa - STAN")
plot_kappa

# look at eta
summary_eta <- stan_lc_df %>%
  rownames_to_column("parameter") %>%
  filter(grepl('eta', parameter)) %>%
  filter(!grepl('beta', parameter)) %>%
  mutate(index = parse_number(parameter)) %>%
  mutate(true_eta = obs.lc$eta)

plot_eta <- ggplot(data=summary_eta) +
  geom_line(aes(x=index, y=mean, color="estimated")) + 
  geom_line(aes(x=index, y=`2.5%`, color="estimated"), alpha=0.5) + 
  geom_line(aes(x=index, y=`97.5%`, color="estimated"), alpha=0.5) +
  geom_line(aes(x=index, y=true_eta, color="true value"), alpha = 0.5) +
  ggtitle("Eta - STAN")
plot_eta

###    ----   Comparing results from inlabru and STAN   ----

plot.comparison <- function(underlying.effects.lc, res.inlabru,
                            fit_summary.stan.lc, summary_alpha,
                            summary_beta, summary_kappa, summary_eta){
  palette.basis <- c('#70A4D4', '#ECC64B', '#93AD80', '#da9124', '#696B8D',
                     '#3290c1', '#5d8060', '#D7B36A', '#826133', '#A85150')
  
  obs <- underlying.effects.lc$obs
  alpha_true <- underlying.effects.lc$alpha.true
  beta_true <- underlying.effects.lc$beta.true
  kappa_true <- underlying.effects.lc$kappa.true
  phi_true <- underlying.effects.lc$phi.true
  intercept_true <- underlying.effects.lc$age.intercept.true
  tau.beta <- underlying.effects.lc$tau.beta.true
  tau.epsilon <- underlying.effects.lc$tau.epsilon.true

  p.prec.alpha <- ggplot(data.frame(res.inlabru$marginals.hyperpar) %>%
                           filter(Precision.for.alpha.x < 35)) + 
    geom_area(aes(x = Precision.for.alpha.x, y = Precision.for.alpha.y),fill = palette.basis[1], alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.hyperpar, aes(xintercept = mean[1]), color = palette.basis[1]) + 
    geom_vline(aes(xintercept = 1/sqrt(fit_summary.stan.lc$mean[1])), color=palette.basis[2]) + 
    labs(x = "Value of precision", y = " ", title = "Precision for alpha - comparison")
  
  p.prec.beta <- ggplot(data.frame(res.inlabru$marginals.hyperpar) %>%
                          filter(Precision.for.beta.x < 200)) + 
    geom_area(aes(x = Precision.for.beta.x, y = Precision.for.beta.y, fill = "Inlabru"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.hyperpar, aes(xintercept = mean[2], color = "Inlabru", fill = "Inlabru")) + 
    geom_vline(aes(xintercept = tau.beta, color = "True value", fill = "True value")) + 
    geom_vline(aes(xintercept = 1/sqrt(fit_summary.stan.lc$mean[2]), color="STAN")) + 
    scale_color_manual(name = " ", values = palette.basis) + 
    scale_fill_manual(name = " ", values = palette.basis) +
    labs(x = "Value of precision", y = " ", title = "Precision for beta - comparsion")
  
  p.prec.kappa <- ggplot(data.frame(res.inlabru$marginals.hyperpar) %>%
                           filter(Precision.for.kappa.x < 175)) + 
    geom_area(aes(x = Precision.for.kappa.x, y = Precision.for.kappa.y), alpha = 0.4, fill = palette.basis[1]) + 
    geom_vline(data = res.inlabru$summary.hyperpar, aes(xintercept = mean[3]), color = palette.basis[1]) + 
    geom_vline(aes(xintercept = 1/sqrt(fit_summary.stan.lc$mean[3])), color=palette.basis[2]) + 
    labs(x = "Value of precision", y = " ", title = "Precision for kappa - comparison")
  
  p.prec.epsilon <- ggplot(data.frame(res.inlabru$marginals.hyperpar) %>%
                             filter(Precision.for.epsilon.x < 150000)) + 
    geom_area(aes(x = Precision.for.epsilon.x, y = Precision.for.epsilon.y, fill = "Inlabru"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.hyperpar, aes(xintercept = mean[4], color = "Inlabru", fill = "Inlabru")) + 
    geom_vline(aes(xintercept = 1/sqrt(fit_summary.stan.lc$mean[4]), color="STAN")) +
    geom_vline(aes(xintercept = tau.epsilon, color = "True value", fill = "True value")) + 
    scale_color_manual(name = " ", values = palette.basis) + 
    scale_fill_manual(name = " ", values = palette.basis) +
    labs(x = "Value of precision", y = " ", title = "Precision for epsilon - comparison")
  
  #   ----   plot results together   ---- 
  
  p.int.c <- ggplot(data.frame(res.inlabru$marginals.fixed)) + 
    geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Inlabru")) + 
    geom_vline(data = fit_summary.stan.lc, aes(xintercept = mean[6], color = "STAN")) +
    geom_vline(data = obs, aes(xintercept = intercept_true, color = "True value")) + 
    scale_color_manual(name = " ", values = palette.basis) + 
    scale_fill_manual(name = " ", values = palette.basis) +
    labs(x = "Value of intercept", y = " ", title = "Intercept - comparison")
  
  
  data.alpha = cbind(res.inlabru$summary.random$alpha, alpha.true = alpha_true[res.inlabru$summary.random$alpha$ID + 1])
  p.alpha.c <- ggplot(data = data.alpha, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
    geom_point(aes(y = mean, color = "Inlabru", fill = "Inlabru")) + 
    geom_point(data = summary_alpha, aes(x = index - 1, y = mean, color = "STAN", fill = "STAN")) +
    geom_line(data = summary_alpha, aes(x = index - 1, y = `2.5%`, color = "STAN", fill = "STAN"), alpha = 0.6) +
    geom_line(data = summary_alpha, aes(x = index - 1, y = `97.5%`, color = "STAN", fill = "STAN"), alpha = 0.6) +
    geom_point(aes(y = alpha.true, color = "True value", fill = "True value")) + 
    scale_color_manual(name = "",
                       values = palette.basis ) +
    scale_fill_manual(name = "",
                      values = palette.basis ) +
    labs(title="Alpha - comparison", x = "x", y='')
  
  data.beta = cbind(res.inlabru$summary.random$beta, beta.true = beta_true[res.inlabru$summary.random$beta$ID + 1])
  p.beta.c <- ggplot(data = data.beta, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
    geom_point(data = summary_beta, aes(x = index - 1, y = mean, color = "STAN", fill = "STAN")) +
    geom_line(data = summary_beta, aes(x = index - 1, y = `2.5%`, color = "STAN", fill = "STAN"), alpha = 0.6) +
    geom_line(data = summary_beta, aes(x = index - 1, y = `97.5%`, color = "STAN", fill = "STAN"), alpha = 0.6) +
    geom_point(aes(y = beta.true, color = "True value", fill = "True value")) + 
    geom_point(aes(y = mean, color = "Inlabru", fill = "Inlabru")) + 
    scale_color_manual(name = "",
                       values = palette.basis ) +
    scale_fill_manual(name = "",
                      values = palette.basis ) +
    labs(x = "x", y = "beta", title = "Beta - comparison")
  
  data.kappa = cbind(res.inlabru$summary.random$kappa, kappa.true = kappa_true[res.inlabru$summary.random$kappa$ID + 1])
  p.kappa.c <- ggplot(data = data.kappa, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
    geom_point(aes(y = kappa.true, color = "True undrifted", fill = "True undrifted")) + 
    geom_point(aes(y = mean, color = "Inlabru", fill = "Inlabru")) + 
    geom_point(data = summary_kappa, aes(x = index - 1, y = mean, color = "STAN", fill = "STAN")) +
    geom_line(data = summary_kappa, aes(x = index - 1, y = `2.5%`, color = "STAN", fill = "STAN"), alpha = 0.6) +
    geom_line(data = summary_kappa, aes(x = index - 1, y = `97.5%`, color = "STAN", fill = "STAN"), alpha = 0.6) +
    scale_color_manual(name = "",
                       values = palette.basis ) +
    scale_fill_manual(name = "",
                      values = palette.basis ) +
    labs(x = "t", y = "kappa", title = "Kappa - comparison")
  
  p.phi.c <- ggplot(data.frame(res.inlabru$marginals.fixed)) + 
    geom_area(aes(x = phi.x, y = phi.y, fill = "Inlabru"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[2], color = "Inlabru", fill="Inlabru")) + 
    geom_vline(aes(xintercept = phi_true, color="True", fill="True")) +
    geom_vline(data=fit_summary.stan.lc, aes(xintercept=mean[5], color="STAN", fill="STAN")) + 
    geom_vline(data=fit_summary.stan.lc, aes(xintercept=`2.5%`[5], color="STAN", fill="STAN"), alpha=0.5) +
    geom_vline(data=fit_summary.stan.lc, aes(xintercept=`97.5%`[5], color="STAN", fill="STAN"), alpha=0.5) +
    scale_color_manual(name = " ", values = palette.basis) + 
    scale_fill_manual(name = " ", values = palette.basis) +
    labs(x = "Value of phi", y = " ", title = "Phi - comparison")
  
  data.eta <- data.frame(eta.sim = res.inlabru$summary.linear.predictor$mean[1:length(obs$eta)]) %>%
    mutate(inlabru.975 = res.inlabru$summary.linear.predictor$`0.975quant`[1:length(obs$eta)]) %>%
    mutate(inlabru.025 = res.inlabru$summary.linear.predictor$`0.025quant`[1:length(obs$eta)]) %>%
    mutate(true.eta = obs$eta) %>%
    mutate(eta.stan = summary_eta$mean) %>%
    mutate(stan.975 = summary_eta$`2.5%`) %>%
    mutate(stan.025 = summary_eta$`97.5%`) %>%
    mutate(x = obs$x, t = obs$t, xt = obs$xt)
  
  print(res.inlabru$summary.linear.predictor)
  print(res.inlabru$summary.linear.predictor$`0.025quant`[1:length(obs$eta)])
  print(res.inlabru$summary.linear.predictor$`0.975quant`[1:length(obs$eta)])
  print(data.eta)
  
  p.eta.c <- ggplot(data = data.eta) +
    geom_line(aes(x=obs$xt, y = eta.sim, color="Inlabru")) +
    geom_line(data=summary_eta, aes(x=index - 1, y=mean, color="STAN")) +
    geom_line(data=obs, aes(x=xt, y = eta, color="True value")) +
    scale_color_manual(name = " ", values = palette.basis) + 
    labs(x=" ", y="Eta", title="Eta- comparison")
  
  p.eta.t <- ggplot(data = data.eta) + 
    geom_line(aes(x = x, y = eta.sim, color = "Inlabru")) +
    geom_line(aes(x = x, y = inlabru.975, color = "Inlabru"), alpha = 0.5) +
    geom_line(aes(x = x, y = inlabru.025, color = "Inlabru"), alpha = 0.5) +
    geom_line(aes(x = x, y = eta.stan, color = "STAN")) +
    geom_line(aes(x = x, y = stan.975, color = "STAN"), alpha = 0.5) +
    geom_line(aes(x = x, y = stan.025, color = "STAN"), alpha = 0.5) +
    geom_line(aes(x = x, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.eta.x <- ggplot(data = data.eta) + 
    geom_line(aes(x = t, y = eta.sim, color = "Inlabru")) +
    geom_line(aes(x = t, y = inlabru.975, color = "Inlabru"), alpha = 0.5) +
    geom_line(aes(x = t, y = inlabru.025, color = "Inlabru"), alpha = 0.5) +
    geom_line(aes(x = t, y = eta.stan, color = "STAN")) +
    geom_line(aes(x = t, y = stan.975, color = "STAN"), alpha = 0.5) +
    geom_line(aes(x = t, y = stan.025, color = "STAN"), alpha = 0.5) +
    geom_line(aes(x = t, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each age") + 
    facet_wrap(~x)
  
  plots <- list(p.prec.alpha = p.prec.alpha, 
                p.prec.beta.compared = p.prec.beta, 
                p.prec.kappa.compared = p.prec.kappa,
                p.prec.epsilon.compared = p.prec.epsilon,
                p.alpha.compared = p.alpha.c,
                p.beta.compared = p.beta.c,
                p.kappa.compared = p.kappa.c,
                p.intercept.compared = p.int.c,
                p.phi.compared = p.phi.c,
                p.eta.compared = p.eta.c,
                p.eta.compared.t = p.eta.t,
                p.eta.compared.x = p.eta.x)
}

plots.comparison <- plot.comparison(underlying.effects.lc, res.inlabru.lc.1, 
                                    stan_lc_df, summary_alpha,
                                    summary_beta, summary_kappa, summary_eta)

plots.comparison$p.prec.alpha
plots.comparison$p.prec.beta.compared
plots.comparison$p.prec.kappa.compared
plots.comparison$p.prec.epsilon.compared

plots.comparison$p.alpha.compared
plots.comparison$p.beta.compared
plots.comparison$p.kappa.compared
plots.comparison$p.intercept.compared
plots.comparison$p.phi.compared
plots.comparison$p.eta.compared
plots.comparison$p.eta.compared.t
plots.comparison$p.eta.compared.x

save.comparison.plots <- function(plot, name, path){
  #'plot: <gg object>
  #'name: <string>  on the format '<name>.png'
  #'
  ggsave(paste(name, '.png', sep=""),
         plot = plot,
         device = "png",
         path = path,
         height = 5, width = 8, 
         dpi = "retina"
  )
  ggsave(paste(name, '.pdf', sep=""),
         plot = plot,
         device = "pdf",
         path = path,
         height = 5, width = 8, 
         dpi = "retina"
  )
}

save.comparison.plots(plots.comparison$p.prec.alpha, 'prec_alpha', path.to.folder)
save.comparison.plots(plots.comparison$p.prec.beta.compared, 'prec_beta', path.to.folder)
save.comparison.plots(plots.comparison$p.prec.kappa.compared, 'prec_kappa', path.to.folder)
save.comparison.plots(plots.comparison$p.prec.epsilon.compared, 'prec_epsilon', path.to.folder)

save.comparison.plots(plots.comparison$p.alpha.compared, 'alpha', path.to.folder)
save.comparison.plots(plots.comparison$p.beta.compared, 'beta', path.to.folder)
save.comparison.plots(plots.comparison$p.kappa.compared, 'kappa', path.to.folder)
save.comparison.plots(plots.comparison$p.intercept.compared, 'intercept', path.to.folder)
save.comparison.plots(plots.comparison$p.phi.compared, 'phi', path.to.folder)
save.comparison.plots(plots.comparison$p.eta.compared, 'eta', path.to.folder)
save.comparison.plots(plots.comparison$p.eta.compared.t, 'eta_t', path.to.folder)
save.comparison.plots(plots.comparison$p.eta.compared.x, 'eta_x', path.to.folder)


#    ----   Old code: don´t know if you still need it yet.   ----

obs <- underlying.effects$obs

obs.pred = data.frame(x = obs$x, t = obs$t, xt = obs$xt, c = obs$c,
                      x.c = obs$x.c, t.c = obs$t.c, E = obs$E, Y = obs$Y,
                      age.intercept.obs = obs$age.intercept,
                      alpha.obs = obs$alpha, beta.obs = obs$beta,
                      phi.t.obs = obs$phi.t, kappa.obs = obs$kappa,
                      gamma.obs = obs$gamma, eta.obs = obs$eta,
                      epsilon.obs = obs$epsilon
                      )

# posterior distribution of beta*phi*t
#beta.phi.posterior <- predict(res.inlabru, data = NULL, formula = ~ beta*phi*t)
beta.phi.posterior <- predict(res.inlabru, data = obs.pred , formula = ~ beta*phi + beta*kappa)

ggplot(data = beta.phi.posterior) +
  geom_line(aes(x = xt, y = mean, color = "Predict")) + 
  geom_line(aes(x = xt, y = beta.obs*phi.t.obs + beta.obs*kappa.obs, color = "Observed"))

alpha.beta.kappa.posterior <- predict(res.inlabru, data = obs.pred , formula = ~ alpha + beta*phi + beta*kappa)

ggplot(data = alpha.beta.kappa.posterior) +
  geom_point(aes(x = xt, y = mean, color = "Predict"), size = 1) + 
  geom_point(aes(x = xt, y = alpha.obs + beta.obs*phi.t.obs + beta.obs*kappa.obs, color = "Observed"), size = 1)

alpha.beta.kappa.gamma.posterior <- predict(res.inlabru, data = obs.pred , formula = ~ alpha + beta*phi + beta*kappa + gamma)

ggplot(data = alpha.beta.kappa.gamma.posterior) +
  geom_point(aes(x = xt, y = mean, color = "Predict")) + 
  geom_point(aes(x = xt, y = alpha.obs + beta.obs*phi.t.obs + beta.obs*kappa.obs + gamma.obs, color = "Observed"))

# observe shioft for this configuration:
int.alpha.beta.kappa.gamma.posterior <- predict(res.inlabru, data = obs.pred , formula = ~ Int + alpha + beta*phi + beta*kappa + gamma)

ggplot(data = int.alpha.beta.kappa.gamma.posterior) +
  geom_point(aes(x = xt, y = mean, color = "Predict")) + 
  geom_point(aes(x = xt, y = age.intercept.obs + alpha.obs + beta.obs*phi.t.obs + beta.obs*kappa.obs + gamma.obs, color = "Observed"))

int.alpha.beta.kappa.gamma.epsilon.posterior <- predict(res.inlabru, data = obs.pred , formula = ~ Int + alpha + beta*phi + beta*kappa + gamma + epsilon)

ggplot(data = int.alpha.beta.kappa.gamma.epsilon.posterior) +
  geom_point(aes(x = xt, y = mean, color = "Predict")) + 
  geom_point(aes(x = xt, y = age.intercept.obs + alpha.obs + beta.obs*phi.t.obs + beta.obs*kappa.obs + gamma.obs  + epsilon.obs, color = "Observed"))

# attempt removing intercept from calculation
alpha.beta.kappa.gamma.epsilon.posterior <- predict(res.inlabru, data = obs.pred , formula = ~ alpha + beta*phi + beta*kappa + gamma + epsilon)

ggplot(data = alpha.beta.kappa.gamma.epsilon.posterior) +
  #geom_line(aes(x = xt, y = mean, color = "Predict wo Int")) + 
  geom_line(data = int.alpha.beta.kappa.gamma.epsilon.posterior,
             aes(x = xt, y = mean, color = "Predict w Int")) +
  # geom_line(aes(x = xt,
  #                y = age.intercept.obs + alpha.obs + beta.obs*phi.t.obs +
  #                  beta.obs*kappa.obs + gamma.obs  + epsilon.obs, color = "Observed")) +
  # geom_line(aes(x = xt,
  #                y = alpha.obs + beta.obs*phi.t.obs + beta.obs*kappa.obs +
  #                  gamma.obs  + epsilon.obs, color = "Observed wo age intercept")) +
  geom_line(aes(x = xt,
                 y = age.intercept.obs + alpha.obs + beta.obs*phi.t.obs +
                   beta.obs*kappa.obs + gamma.obs  + epsilon.obs -
                   mean(beta.obs*phi.t.obs),
                 color = "Observed mean(period effect extracted)")) +
  geom_line(aes(x = xt,
                y = age.intercept.obs + alpha.obs + beta.obs*phi.t.obs +
                  beta.obs*kappa.obs + gamma.obs  + epsilon.obs + 0.6296,
                color = "Observed mean(mean(period effect extracted))"))
  
all.posteroops <- predict(res.inlabru, data = obs.pred , formula = NULL)

ggplot() +
  geom_point(data = data.frame(mean = all.posteroops$Predictor$mean[1:8000]),
             aes(x = obs.pred$xt, y = mean, color = "Predict")) + 
  geom_point(data = obs.pred, aes(x = xt, y = exp(eta.obs), color = "Observed"))


