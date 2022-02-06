#' This program investigates inference with:
#' - The traditional_lc_log_prec version - Gaussian lc-model. 
#' All random effects are modelled as iid, in inlabru and in stan
#' Constraints are imposed as usual in inlabru, and with soft constraints in stan
#' They hyperparameters are fixed. 
#' 

#   ----   Load libraries and set workspace   ----
library("tidyverse")
library("inlabru")
library("ggplot2")
library("INLA")
library("patchwork")
library("rstan")

#TODO: Change to the location of the folder
setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic\ data/Step_by_step_results_2/poiss_gp_rw1_stand_alones")

# setting paths, some are deprecated
investigation.path <- file.path("v4") #TODO: "v4" should be changed to the name of the folder where you store all files for the v4 data. 
output.path <- investigation.path
stan.output <- output.path  # deprecated

#   ----    Retrieve the data   ----

synthetic.male.lung.v4 <- function(){
  obs <- read.csv(file.path(investigation.path, "synthetic_male_lung_4.csv"))
  
  obs.trad <- obs %>% 
    select(c(x, t, xt, age.int, year, x.c, alpha, beta, kappa, intercept, epsilon,
             eta, tau.alpha, tau.beta, tau.kappa, tau.epsilon, E, Y)) %>%
    mutate(x = x+1, x.c = x.c + 1, t = t + 1, xt = xt + 1)
  
  underlying.effects <- list(obs = obs.trad, nx = 18, nt = 18,
                             alpha.true = {obs %>% filter(t == 0)}$alpha,
                             beta.true = {obs %>% filter(t == 0)}$beta,
                             kappa.true = {obs %>% filter(x == 0)}$kappa,
                             intercept = unique(obs$intercept),
                             age.intercept.true = unique(obs$intercept),
                             tau.alpha.true = unique(obs$tau.alpha),
                             tau.beta.true = unique(obs$tau.beta),
                             tau.kappa.true = unique(obs$tau.kappa),
                             tau.epsilon.true = unique(obs$tau.epsilon))
  
  return(list(obs = obs.trad, underlying.effects = underlying.effects))
}

# We use this data for both inlabru and stan
config.data <- synthetic.male.lung.v4()
obs <- config.data$obs
underlying.effects <- config.data$underlying.effects

####   ----    Plot the data   ----

obs.zeroes <- obs %>% filter(Y == 0)
print("Number of zero counts of data v4: ")
print(length(obs.zeroes$Y))

p.Y <- ggplot(obs) + 
  geom_tile(aes(x = x, y = t, fill = Y)) + 
  labs(title = "Observed deaths Y of v4", x = "x", y = "t")
p.Y; ggsave("data_Y.pdf", p.Y, path = output.path)

p.Y.2 <- ggplot(obs) + 
  geom_line(aes(x = xt, y = Y)) + 
  labs(title = "Observed deaths Y of v4", x = "x,t", y = "")
p.Y.2; ggsave("data_Y_2.pdf", p.Y.2, path = output.path)

#   ----   Run STAN analysis   ----
# Running traditional lc version of stan, 
# implemented with log-precisions, random effects as iid, and no constraints

run_stan <- function(stan_program, obs, chains, warmup, iter, output.path){
  
  input_stan <- list(
    X=length(unique(obs$x)),
    T=length(unique(obs$t)),
    x = obs$x,
    t = obs$t,
    E = as.integer(obs$E),
    Y = as.integer(obs$Y),
    nx = length(unique(obs$x)),
    nt = length(unique(obs$t))
  )
  
  stan.run.time <- system.time({
    stan_fit <- stan(
      file = stan_program,
      data = input_stan,
      chains = chains,
      iter = iter,
      warmup = warmup,
      refresh = iter%/%10,
      seed = 123
    )
  }
  )
  
  write.table(input_stan, file = file.path(output.path, "input_stan.txt"))
  write.table(list(program = stan_program, chains = chains, iter = iter, warmup = warmup), file = file.path(output.path, "info_stan.txt"))
  write.table(summary(stan.run.time), file= file.path(output.path, "stan_run_time.txt"))
  
  return(stan_fit)
}

# TODO: Update to the desired amount of cainns, warmup and iterations
stan_fit <- run_stan(
  stan_program = file.path(output.path, "stan_poiss_gp_rw1_4.stan"),
  obs = obs, chains=2, warmup = 1000, iter = 10000, output.path = stan.output)

ggsave("trace.pdf", traceplot(stan_fit, pars = c("tau_alpha", "tau_beta", "tau_kappa",
                                                 "tau_epsilon", "eta[1]", "alpha[1]",
                                                 "beta[1]", "kappa[1]", "eta[54]",
                                                 "alpha[9]", "beta[9]", "kappa[9]")),
       path = output.path, dpi = "retina", width = 8, height = 5)

stan.summary <- data.frame(summary(stan_fit))
save(stan.summary, file = file.path(output.path, "stan_summary.Rda"))
#load(file = file.path(output.path, "stan_summary.Rda"))
list_of_draws <- rstan::extract(stan_fit)
save(list_of_draws, file = file.path(output.path, "list_of_draws.RData"))
#load(file = file.path(output.path, "list_of_draws.RData"))


inlabru.poiss.fh.rw1 <- function(obs, output.path, max_iter=30, write = T){
  #'Implements inlabru analysis for lc model, fixing the precisions and modelling all random effects as iid
  #'
  #'@param obs: Contains the observed data and the real underlying random effects
  #'@param max_iter (int): maximum number of iterations in inlabru
  
  nx = length(unique(obs$x))
  nt = length(unique(obs$t))
  
  loggamma.prior <- list(prec = list(prior = 'loggamma', param = c(1,0.00005), initial = log(1)))
  loggamma.prior.high.variance <- list(prec = list(prior = 'loggamma', param = c(1,0.005), initial = log(1)))
  
  # For constraining of beta:
  A.beta = matrix(1, nrow = 1, ncol = nx)
  e.beta = 1
  
  comp = ~ -1 +
    alpha(x, model = "rw1", hyper = loggamma.prior, constr = FALSE, scale.model = F) +
    beta(x.c, model = "iid", hyper = loggamma.prior, extraconstr = list(A = A.beta, e = e.beta)) +
    kappa(t, model = "rw1", hyper = loggamma.prior.high.variance, constr = TRUE, scale.model = F) +
    epsilon(xt, model = "iid", hyper = loggamma.prior, constr = FALSE)
  
  formula = Y ~ alpha + beta*kappa + epsilon
  
  likelihood = like(formula = formula, family = "poisson", data = obs, E = obs$E)
  
  c.compute <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE, return.marginals.predictor = TRUE)  # control.compute
  
  inlabru.run.time <- system.time(
    {
      res.inlabru = bru(components = comp,
                        likelihood,
                        options = list(verbose = F,
                                       bru_verbose = 4,
                                       num.threads = "1:1",
                                       control.compute = c.compute,
                                       bru_max_iter=max_iter,
                                       control.predictor = list(link = 1)
                        ))
    }
  )
  
  if(write){
    write.table(summary(inlabru.run.time), file = file.path(output.path, "runtime_inlabru.txt"))
  }
  
  return(res.inlabru)
}


res.inlabru <- inlabru.poiss.fh.rw1(obs, output.path = output.path, max_iter = 100, write = T)

source(file.path(investigation.path, "plotters.R"))

#   ----   Produce stan summaries    ----

summary_alpha <- stan.summary %>%
  rownames_to_column("parameter") %>%
  filter(grepl('alpha', parameter)) %>%
  filter(!grepl('tau_alpha', parameter)) %>%
  filter(!grepl('theta_alpha', parameter)) %>%
  mutate(index = parse_number(parameter))

summary_beta <- stan.summary %>%
  rownames_to_column("parameter") %>%
  filter(grepl('beta', parameter)) %>%
  filter(!grepl('tau_beta', parameter)) %>%
  filter(!grepl('theta_beta', parameter)) %>%
  mutate(index = parse_number(parameter))

summary_kappa <- stan.summary %>%
  rownames_to_column("parameter") %>%
  filter(grepl('kappa', parameter)) %>%
  filter(!grepl('tau_kappa', parameter)) %>%
  filter(!grepl('theta_kappa', parameter)) %>%
  mutate(index = parse_number(parameter))

summary_eta <- stan.summary %>%
  rownames_to_column("parameter") %>%
  filter(grepl('eta', parameter)) %>%
  filter(!grepl('beta', parameter)) %>%
  filter(!grepl('theta', parameter)) %>%
  mutate(index = parse_number(parameter))

#   ----   Plot Random effects and summary of predictor   ----

palette <- c('#70A4D4', '#ECC64B', '#93AD80', '#da9124', '#696B8D',
             '#3290c1', '#5d8060', '#D7B36A', '#826133', '#A85150')

p.predictor.summary <- ggplot(data = data.frame(xt = obs$xt,
                                                mean = res.inlabru$summary.linear.predictor$mean[1:324],
                                                X0.025quant = res.inlabru$summary.linear.predictor$`0.025quant`[1:324],
                                                X0.975quant = res.inlabru$summary.linear.predictor$`0.975quant`[1:324])) + 
  geom_point(data = underlying.effects$obs, aes(x = xt, y = eta, color = "True", fill = "True"), alpha = 0.5) + 
  geom_ribbon(aes(x = xt, ymin = X0.025quant, ymax = X0.975quant, color = "Inlabru", fill = "Inlabru"), alpha = 0.3) + 
  geom_ribbon(data = summary_eta, aes(x = index, ymin = summary.2.5., ymax = summary.97.5., color = "Stan", fill = "Stan"), alpha = 0.3) + 
  geom_point(aes(x = xt, y = mean, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
  geom_point(data = summary_eta, aes(x = index, y = summary.mean, fill = "Stan", color = "Stan"), alpha = 0.5) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Predictor", x = "x, t", y = "")

ggsave("predictor.pdf", p.predictor.summary, path = output.path, dpi = "retina", height = 5, width = 8)

p.alpha <- ggplot() + 
  geom_ribbon(data = data.frame(res.inlabru$summary.random$alpha), aes(x = ID, ymin = X0.025quant, ymax = X0.975quant, fill = "Inlabru", color = "Inlabru"), alpha = 0.3) + 
  geom_ribbon(data = summary_alpha, aes(x = index, ymin = summary.2.5., ymax = summary.97.5., fill = "Stan", color = "Stan"), alpha = 0.3) +
  geom_point(data = obs %>% filter(t == 1), aes(x = x, y = alpha + intercept, color = "True", fill = "True"), alpha = 0.7) + 
  geom_point(data = data.frame(res.inlabru$summary.random$alpha), aes(x = ID, y = mean, fill = "Inlabru", color = "Inlabru"), alpha = 0.7) + 
  geom_point(data = summary_alpha, aes(x = index, y = summary.mean, color = "Stan", fill = "Stan"), alpha = 0.7) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Alpha", x = "x", y = "")

p.beta <- ggplot() +
  geom_errorbar(data = data.frame(res.inlabru$summary.random$beta), aes(x = ID, ymin = X0.025quant, ymax = X0.975quant, fill = "Inlabru", color = "Inlabru"), alpha = 0.7) +
  geom_errorbar(data = summary_beta, aes(x = index, ymin = summary.2.5., ymax = summary.97.5., fill = "Stan", color = "Stan"), alpha = 0.7) +
  geom_point(data = obs %>% filter(t == 1), aes(x = x, y = beta, color = "True", fill = "True"), alpha = 0.7) +
  geom_point(data = data.frame(res.inlabru$summary.random$beta), aes(x = ID, y = mean, fill = "Inlabru", color = "Inlabru"), alpha = 0.7) +
  geom_point(data = summary_beta, aes(x = index, y = summary.mean, color = "Stan", fill = "Stan"), alpha = 0.7) +
  theme_classic() +
  scale_color_manual(name = "", values = palette) +
  scale_fill_manual(name = "", values = palette) +
  labs(title = "Beta", x = "x", y = "")


p.kappa <- ggplot() +
  geom_ribbon(data = data.frame(res.inlabru$summary.random$kappa), aes(x = ID, ymin = X0.025quant, ymax = X0.975quant, fill = "Inlabru", color = "Inlabru"), alpha = 0.3) + 
  geom_ribbon(data = summary_kappa, aes(x = index, ymin = summary.2.5., ymax = summary.97.5., fill = "Stan", color = "Stan"), alpha = 0.3) +
  geom_point(data = obs %>% filter(x == 1), aes(x = t, y = kappa, color = "True", fill = "True"), alpha = 0.7) + 
  geom_point(data = data.frame(res.inlabru$summary.random$kappa), aes(x = ID, y = mean, fill = "Inlabru", color = "Inlabru"), alpha = 0.7) + 
  geom_point(data = summary_kappa, aes(x = index, y = summary.mean, color = "Stan", fill = "Stan"), alpha = 0.7) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Kappa", x = "t", y = "")

p.random <- (p.alpha | p.beta) / (p.kappa) + plot_layout(guides = "collect")

ggsave("random.pdf", p.random, path = output.path, dpi = "retina", height = 4, width = 6.4)


#   ----   Sample predictor   ----

stan.predictor.df <- data.frame(list_of_draws$eta)
plot.predictor.inlabru.stan.compared(res.inlabru, stan.predictor.df, path.to.storage = output.path, a45=F)

#   ----   Plot marginals of random effects   ----

stan.kappa.df <- data.frame(list_of_draws$kappa)
plot.kappa.inlabru.stan.compared(res.inlabru, stan.kappa.df, path.to.storage = output.path)

stan.alpha.df <- data.frame(list_of_draws$alpha)
plot.alpha.inlabru.stan.compared(res.inlabru, stan.alpha.df, path.to.storage = output.path)

stan.epsilon.df <- data.frame(list_of_draws$epsilon)
plot.epsilon.inlabru.stan.compared(res.inlabru, stan.epsilon.df, path.to.storage = output.path)

stan.beta.df <- data.frame(list_of_draws$beta)
plot.beta.inlabru.stan.compared(res.inlabru, stan.beta.df, path.to.storage = output.path)

#   ----   Specifically check the predictors at xt = 54:   ----

pred.54.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.054)

p.pred.54 <- ggplot(pred.54.inlabru) + 
  geom_area(aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
  geom_density(data = stan.predictor.df, aes(x = X54, fill = "Stan", color = "Stan"), alpha = 0.5) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Predictor at xt=54", x = "", y = "")
p.pred.54

save.figure(p.pred.54, name = "predictor_54", path = output.path, png= F)

pred.36.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.036)

p.pred.36 <- ggplot(pred.36.inlabru) + 
  geom_area(aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
  geom_density(data = stan.predictor.df, aes(x = X36, fill = "Stan", color = "Stan"), alpha = 0.5) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Predictor at xt=36", x = "", y = "")
p.pred.36

save.figure(p.pred.36, name = "predictor_36", path = output.path, png= F)

pred.1.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.001)

p.pred.1 <- ggplot(pred.1.inlabru) + 
  geom_area(aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
  geom_density(data = stan.predictor.df, aes(x = X1, fill = "Stan", color = "Stan"), alpha = 0.5) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  labs(title = "Predictor at xt=1", x = "", y = "")
p.pred.1

save.figure(p.pred.1, name = "predictor_1", path = output.path, png= F)

#   ----   Plot densities of hyperparameters   ----

p.tau.alpha <- ggplot() +
  geom_density(data=data.frame("x" = list_of_draws$tau_alpha), aes(x = x, color = "Stan", fill = "Stan"), alpha = 0.2) +
  geom_area(data = data.frame(res.inlabru$marginals.hyperpar$`Precision for alpha`) %>% filter(x < 7), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) +
  theme_classic() +
  scale_color_manual(name = "", values = palette) +
  scale_fill_manual(name = "", values = palette) +
  labs(title = "Tau alpha", x = "", y = "")

p.theta.alpha <- ggplot() +
  geom_density(data=data.frame("x" = list_of_draws$theta_alpha), aes(x = x, color = "Stan", fill = "Stan"), alpha = 0.2) +
  geom_area(data = data.frame(res.inlabru$internal.marginals.hyperpar$`Log precision for alpha`), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) +
  theme_classic() +
  scale_color_manual(name = "", values = palette) +
  scale_fill_manual(name = "", values = palette) +
  labs(title = "Theta alpha", x = "", y = "")

p.tau.beta <- ggplot() +
  geom_density(data=data.frame("x" = list_of_draws$tau_beta), aes(x = x, color = "Stan", fill = "Stan"), alpha = 0.2) +
  geom_area(data = data.frame(res.inlabru$marginals.hyperpar$`Precision for beta`) %>% filter(x < 500), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) +
  theme_classic() +
  scale_color_manual(name = "", values = palette) +
  scale_fill_manual(name = "", values = palette) +
  labs(title = "Tau beta", x = "", y = "")

p.theta.beta <- ggplot() +
  geom_density(data=data.frame("x" = list_of_draws$theta_beta), aes(x = x, color = "Stan", fill = "Stan"), alpha = 0.2) +
  geom_area(data = data.frame(res.inlabru$internal.marginals.hyperpar$`Log precision for beta`), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) +
  theme_classic() +
  scale_color_manual(name = "", values = palette) +
  scale_fill_manual(name = "", values = palette) +
  labs(title = "Theta beta", x = "", y = "")

p.tau.kappa <- ggplot() +
  geom_density(data=data.frame("x" = list_of_draws$tau_kappa) %>% filter(x < 750), aes(x = x, color = "Stan", fill = "Stan"), alpha = 0.2) +
  geom_area(data = data.frame(res.inlabru$marginals.hyperpar$`Precision for kappa`) %>% filter(x < 750), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) +
  theme_classic() +
  scale_color_manual(name = "", values = palette) +
  scale_fill_manual(name = "", values = palette) +
  labs(title = "Tau kappa", x = "", y = "")

p.theta.kappa <- ggplot() +
  geom_density(data=data.frame("x" = list_of_draws$theta_kappa), aes(x = x, color = "Stan", fill = "Stan"), alpha = 0.2) +
  geom_area(data = data.frame(res.inlabru$internal.marginals.hyperpar$`Log precision for kappa`), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) +
  theme_classic() +
  scale_color_manual(name = "", values = palette) +
  scale_fill_manual(name = "", values = palette) +
  labs(title = "Theta kappa", x = "", y = "")

p.tau.eta <- ggplot() +
  geom_density(data=data.frame("x" = list_of_draws$tau_epsilon), aes(x = x, color = "Stan", fill = "Stan"), alpha = 0.2) +
  geom_area(data = data.frame(res.inlabru$marginals.hyperpar$`Precision for epsilon`) %>% filter(x < 600), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) +
  theme_classic() +
  scale_color_manual(name = "", values = palette) +
  scale_fill_manual(name = "", values = palette) +
  labs(title = "Tau epsilon", x = "", y = "")

p.theta.eta <- ggplot() +
  geom_density(data=data.frame("x" = list_of_draws$theta_epsilon), aes(x = x, color = "Stan", fill = "Stan"), alpha = 0.2) +
  geom_area(data = data.frame(res.inlabru$internal.marginals.hyperpar$`Log precision for epsilon`), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.2) +
  theme_classic() +
  scale_color_manual(name = "", values = palette) +
  scale_fill_manual(name = "", values = palette) +
  labs(title = "Theta epsilon", x = "", y = "")

p.tau <- (p.tau.alpha | p.tau.beta)/(p.tau.kappa | p.tau.eta) + plot_layout(guides = "collect")

ggsave("marginals_tau.pdf", p.tau, path=output.path, device = "pdf", dpi="retina", height = 5, width = 8)

