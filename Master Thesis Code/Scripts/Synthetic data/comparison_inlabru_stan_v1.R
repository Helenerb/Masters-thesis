setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data")

source("configuration_v1.R")

# configuration without cohort effect
#underlying.effects.lc <- configuration.v5()
underlying.effects.lc <- configuration.v9()

obs.lc <- underlying.effects.lc$obs

# cohfiguration with cohort effect
underlying.effects.lc.cohort <- configuration.v7()

obs.lc.cohort <- underlying.effects.lc.cohort$obs

source("inlabru_analyses.R")

res.inlabru.lc.1 <- inlabru.lc.1(obs.lc)

res.inlabru.lc.cohort.1 <- inlabru.lc.cohort.1(obs.lc.cohort)

#   ----   plotting results of inlabru fit   ----
source("plot_inlabru_vs_underlying.R")

# plotting results from run with cohort effects:
plots.lc.cohort <- plot.inlabru.vs.underlying.v1(res.inlabru.lc.cohort.1,
                                                 underlying.effects.lc.cohort)
plots.lc.cohort$p.alpha
plots.lc.cohort$p.beta
plots.lc.cohort$p.phi
plots.lc.cohort$p.intercept
plots.lc.cohort$p.kappa
plots.lc.cohort$p.eta
plots.lc.cohort$p.eta.2
plots.lc.cohort$p.eta.t
plots.lc.cohort$p.eta.x
plots.lc.cohort$p.gamma

# plotting results from run with cohort effects:
plots <- plot.inlabru.vs.underlying.v5(res.inlabru.lc.1, underlying.effects.lc)
plots$p.alpha
plots$p.beta
plots$p.phi
plots$p.intercept
plots$p.kappa
plots$p.eta
plots$p.eta.2
plots$p.eta.t
plots$p.eta.x

#   ----   Perform similar analysis in STAN   ----   
library("rstan")

input_stan.lc <- list(
  X=length(unique(obs.lc$x)),
  T=length(unique(obs.lc$t)),
  x=(obs.lc$x + 1),
  t=(obs.lc$t + 1),
  ts = obs.lc$t,
  E = obs.lc$E,
  Y = obs.lc$Y
)

# fit <- stan(
#   file="stan_analysis_lc.stan",
#   data = input_stan.lc,
#   chains=4,
#   warmup = 1000,
#   iter = 4000,
#   refresh = 200,
#   seed=123
# )

# save workspace image 
#save.image("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Workspaces/stan_analysis_lc.RData")
# load workspace image
#load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Workspaces/stan_analysis_lc.RData")

# attempt with testin-stan file with fewer runs - try to speed computations up!

fit <- stan(
  file="stan_analysis_lc_v2.stan",
  data = input_stan.lc,
  chains=4,
  warmup = 2000,
  iter = 6000,
  refresh = 100,
  seed=123
)

pairs(fit, pars=c("tau_alpha", "phi", "intercept", "eta[1]"))

#   ----   plot stan results   ---- 

traceplot(fit, pars=c("eta[1]", "eta[2]", "eta[3]", "eta[4]", "eta[5]",
                      "eta[6]", "eta[7]", "eta[8]", "eta[9]", "eta[10]" ))

traceplot(fit, pars=c("kappa[1]", "kappa[2]", "kappa[3]", "kappa[4]",
                      "kappa[5]", "kappa[6]", "kappa[7]", "kappa[8]",
                      "kappa[9]", "kappa[10]" ))
traceplot(fit, pars=c("beta[1]", "beta[2]", "beta[3]", "beta[4]",
                      "beta[5]", "beta[6]", "beta[7]", "beta[8]",
                      "beta[9]", "beta[10]" ))

traceplot(fit)

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
  geom_point(aes(x = index, y = kappa_drifted, color = "true drifted")) +
  ggtitle("Kappa; time effect with drift - STAN")
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


