# synthetic configuration based on male lung data

# assuming workspace at .../Master Thesis Code

library("tidyverse")
library("ggplot2")

load("Data/population-germany.Rda")

load("Data/lungCancer-germany.Rda")
lung.cancer <- cancer.data

lung.cancer <- lung.cancer %>%
  mutate(female.mr = female/female.t, male.mr = male/male.t)

lung.cancer.male <- lung.cancer %>%
  select(x, x.c, t, xt, cohort, c, age, age.int, year, birth.year, male.t, male, male.mr) %>%
  mutate(E = male.t, Y = male, `mortality rate` = male.mr)

#   ----   perform inlabru analyses   ----
source("Scripts/Functions/inlabru_analyses.R")

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
obs.4 <- read.csv("Data/synthetic_male_lung_4.csv")

underlying.effects.4 <- list(obs = obs.4, nx = 18, nt = 18,
                             alpha.true = {obs.4 %>% filter(t == 0)}$alpha,
                             beta.true = {obs.4 %>% filter(t == 0)}$beta,
                             kappa.true = {obs.4 %>% filter(x == 0)}$kappa,
                             intercept = unique(obs.4$intercept),
                             tau.alpha.true = unique(obs.4$tau.alpha),
                             tau.beta.true = unique(obs.4$tau.beta),
                             tau.kappa.true = unique(obs.4$tau.kappa),
                             tau.epsilon.true = unique(obs.4$tau.epsilon))

source("Scripts/Functions/inlabru_analyses.R")

inlabru.synthetic.male.lung.lc.4 <- inlabru.rw2.lc.2(obs.4, max_iter = 100)

source("Scripts/Synthetic data/plot_inlabru_vs_underlying.R")

plots.summaries.inlabru <- plot.inlabru.vs.underlying.synthetic.cancer(
  inlabru.synthetic.male.lung.lc.4,
  underlying.effects.4,
  path.to.storage = "Scripts/Synthetic\ data/Output/Figures/synthetic_male_lung_lc/v4",
  save=TRUE)

#   ----   v5   ----
# attempt to find a configuration that is of somewhat comparable size to real data
# v4 was ok, try to find configuration with more realistic span for kappa and beta. 
set.seed(26)
alpha.5 <- rw1(inlabru.tau.alpha, 18); alpha.5 = alpha.5 - mean(alpha.5)
beta.5 <- rnorm(18, sd = sqrt(1/inlabru.tau.beta)); beta.5 = beta.5 - mean(beta.5) + 1/18
kappa.5 <- rw2(inlabru.tau.kappa, 18); kappa.5 = kappa.5 - mean(kappa.5)
epsilon.5 <- rnorm(18*18, sd = sqrt(1/inlabru.tau.epsilon))

obs.5 <- data.frame(x = lung.cancer.male$x, t = lung.cancer.male$t, xt = lung.cancer.male$xt, E = lung.cancer.male$E) %>%
  mutate(age.int = lung.cancer.male$age.int, year = lung.cancer.male$year) %>%
  mutate(x.c = x) %>%
  mutate(alpha = alpha.5[x + 1]) %>%
  mutate(beta = beta.5[x + 1]) %>%
  mutate(kappa = kappa.5[t + 1]) %>%
  mutate(intercept = inlabru.intercept) %>%
  mutate(epsilon = epsilon.4[xt + 1]) %>%
  mutate(eta = intercept + alpha + beta*kappa + epsilon) %>%
  mutate(Y = rpois(length(x), E*exp(eta))) %>%
  mutate(mr = Y/E) %>%
  mutate(tau.alpha = inlabru.tau.alpha) %>%
  mutate(tau.beta = inlabru.tau.beta) %>%
  mutate(tau.kappa = inlabru.tau.kappa) %>%
  mutate(tau.epsilon = inlabru.tau.epsilon)

ggplot(obs.5) + geom_line(aes(x = x, y = mr, color = year))
ggplot(obs.5) + geom_line(aes(x = x, y = eta, color = year))
ggplot(obs.5) + geom_line(aes(x = x, y = alpha))
ggplot(obs.5) + geom_line(aes(x = x, y = beta))
ggplot(obs.5) + geom_line(aes(x = t, y = kappa))

write.csv(obs.5, "Data/synthetic_male_lung_5.csv")

underlying.effects.5 <- list(obs = obs.5, nx = 18, nt = 18,
                             alpha.true = {obs.5 %>% filter(t == 0)}$alpha,
                             beta.true = {obs.5 %>% filter(t == 0)}$beta,
                             kappa.true = {obs.5 %>% filter(x == 0)}$kappa,
                             intercept = unique(obs.5$intercept))

source("Scripts/Functions/inlabru_analyses.R")

inlabru.synthetic.male.lung.lc.5 <- inlabru.rw2.lc.2(obs.5, max_iter = 100)

source("Scripts/Synthetic data/plot_inlabru_vs_underlying.R")

plots.summaries.inlabru <- plot.inlabru.vs.underlying.synthetic.cancer(
  inlabru.synthetic.male.lung.lc.5,
  underlying.effects.5,
  path.to.storage = "Scripts/Synthetic\ data/Output/Figures/synthetic_male_lung_lc/v5",
  save=TRUE)
  
###    ----    v6   ----

# attempt 6 - with adjusted stan hyperpars instead

stan.tau.alpha = 2
stan.tau.beta = 1000
stan.tau.kappa = 192

adjusted.tau.alpha = stan.tau.alpha - 1
adjusted.tau.beta = stan.tau.beta
adjusted.tau.kappa = stan.tau.kappa - 150

set.seed(17) # could work, somewhat ok mr


alpha.6 <- rw1(adjusted.tau.alpha, 18); alpha.6 = alpha.6 - mean(alpha.6)
beta.6 <- rnorm(18, sd = sqrt(1/adjusted.tau.beta)); beta.6 = beta.6 - mean(beta.6) + 1/18
kappa.6 <- rw2(adjusted.tau.kappa, 18); kappa.6 = kappa.6 - mean(kappa.6)
epsilon.6 <- rnorm(18*18, sd = sqrt(1/inlabru.tau.epsilon))

obs.6 <- data.frame(x = lung.cancer.male$x, t = lung.cancer.male$t, xt = lung.cancer.male$xt, E = lung.cancer.male$E) %>%
  mutate(age.int = lung.cancer.male$age.int, year = lung.cancer.male$year) %>%
  mutate(x.c = x) %>%
  mutate(alpha = alpha.6[x + 1]) %>%
  mutate(beta = beta.6[x + 1]) %>%
  mutate(kappa = kappa.6[t + 1]) %>%
  mutate(intercept = inlabru.intercept) %>%
  mutate(epsilon = epsilon.6[xt + 1]) %>%
  mutate(eta = intercept + alpha + beta*kappa + epsilon) %>%
  mutate(Y = rpois(length(x), E*exp(eta))) %>%
  mutate(mr = Y/E) %>%
  mutate(tau.alpha = adjusted.tau.alpha) %>%
  mutate(tau.beta = adjusted.tau.beta) %>%
  mutate(tau.kappa = adjusted.tau.kappa) %>%
  mutate(tau.epsilon = inlabru.tau.epsilon)

ggplot(obs.6) + geom_line(aes(x = x, y = mr, color = year))
ggplot(obs.6) + geom_line(aes(x = x, y = eta, color = year))
ggplot(obs.6) + geom_line(aes(x = x, y = alpha))
ggplot(obs.6) + geom_line(aes(x = x, y = beta))
ggplot(obs.6) + geom_line(aes(x = t, y = kappa))

#write.csv(obs.6, "Data/synthetic_male_lung_6.csv")
obs.6 <- read.csv("Data/synthetic_male_lung_6.csv")

underlying.effects.6 <- list(obs = obs.6, nx = 18, nt = 18,
                             alpha.true = {obs.6 %>% filter(t == 0)}$alpha,
                             beta.true = {obs.6 %>% filter(t == 0)}$beta,
                             kappa.true = {obs.6 %>% filter(x == 0)}$kappa,
                             intercept = unique(obs.6$intercept),
                             tau.alpha.true = unique(obs.6$tau.alpha),
                             tau.beta.true = unique(obs.6$tau.beta),
                             tau.kappa.true = unique(obs.6$tau.kappa),
                             tau.epsilon.true = unique(obs.6$tau.epsilon)
                             )

source("Scripts/Functions/inlabru_analyses.R")

inlabru.synthetic.male.lung.lc.6<- inlabru.rw2.lc.2(obs.6, max_iter = 100)

source("Scripts/Synthetic data/plot_inlabru_vs_underlying.R")

plots.summaries.inlabru <- plot.inlabru.vs.underlying.synthetic.cancer(
  inlabru.synthetic.male.lung.lc.6,
  underlying.effects.6,
  path.to.storage = "Scripts/Synthetic\ data/Output/Figures/synthetic_male_lung_lc/v6",
  save=TRUE)

###   ----   v7   ----

# seems like inlabru struggles for low absolute values. 
# first, to show results for non-zero mortalities, we run inlabru for male lung cancer with only ages above 45:
lung.cancer.male.above.45 <- lung.cancer.male %>%
  filter(age.int >= 45)

res.lung.lc.m.a45 <- inlabru.rw2.lc.2(lung.cancer.male.above.45)

source("Scripts/Real data/plot_real_data.R")

plots.stomach.male <- plot.inlabru.real(
  res.lung.lc.m.a45, lung.cancer.male.above.45, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_rw2_lc_a45/male",
  cohort=FALSE)

plot.hypers.inlabru.real(
  res.lung.lc.m.a45, lung.cancer.male.above.45, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_rw2_lc_a45/male",
  cohort=FALSE)

inlabru.tau.alpha.a45 <- res.lung.lc.m.a45$summary.hyperpar$mean[1] 
inlabru.tau.beta.a45 <- res.lung.lc.m.a45$summary.hyperpar$mean[2]
inlabru.tau.kappa.a45 <- res.lung.lc.m.a45$summary.hyperpar$mean[3]
inlabru.tau.epsilon.a45 <- res.lung.lc.m.a45$summary.hyperpar$mean[4]

inlabru.intercept.a45 <- res.lung.lc.m.a45$summary.fixed$mean[1]

set.seed(1)

alpha.7 <- rw1(inlabru.tau.alpha.a45, 9); alpha.7 = alpha.7 - mean(alpha.7)
beta.7 <- rnorm(9, sd = sqrt(1/inlabru.tau.beta.a45)); beta.7 = beta.7 - mean(beta.7) + 1/9
kappa.7 <- rw2(inlabru.tau.kappa.a45, 18); kappa.7 = kappa.7 - mean(kappa.7)
epsilon.7 <- rnorm(9*18, sd = sqrt(1/inlabru.tau.epsilon.a45))

obs.7 <- data.frame(x = lung.cancer.male.above.45$x, t = lung.cancer.male.above.45$t, E = lung.cancer.male.above.45$E) %>%
  mutate(age.int = lung.cancer.male.above.45$age.int, year = lung.cancer.male.above.45$year) %>%
  mutate(xt = seq_along(x)) %>%
  mutate(x.c = x) %>%
  mutate(alpha = alpha.7[x - 9 + 1]) %>%
  mutate(beta = beta.7[x - 9 + 1]) %>%
  mutate(kappa = kappa.7[t + 1]) %>%
  mutate(intercept = inlabru.intercept.a45) %>%
  mutate(epsilon = epsilon.7[xt]) %>%
  mutate(eta = intercept + alpha + beta*kappa + epsilon) %>%
  mutate(e.exp.eta = E*exp(eta)) %>%
  mutate(Y = rpois(length(x), E*exp(eta))) %>%
  mutate(mr = Y/E) %>%
  mutate(tau.alpha = inlabru.tau.alpha.a45) %>%
  mutate(tau.beta = inlabru.tau.beta.a45) %>%
  mutate(tau.kappa = inlabru.tau.kappa.a45) %>%
  mutate(tau.epsilon = inlabru.tau.epsilon.a45)

#write.csv(obs.7, "Data/synthetic_male_lung_7.csv")
obs.7 <- read.csv("Data/synthetic_male_lung_7.csv")
obs.7 <- obs.7 %>% mutate(year.str = as.character(year))

ggplot(obs.7) + geom_line(aes(x = x, y = mr, color = year.str))
ggplot(obs.7) + geom_line(aes(x = x, y = eta, color = year.str))
ggplot(obs.7) + geom_line(aes(x = x, y = alpha))
ggplot(obs.7) + geom_line(aes(x = x, y = beta))
ggplot(obs.7) + geom_line(aes(x = t, y = kappa))

underlying.effects.7 <- list(obs = obs.7, nx = 9, nt = 18,
                             alpha.true = c(rep(0,9), {obs.7 %>% filter(t == 0)}$alpha),
                             beta.true = c(rep(0,9), {obs.7 %>% filter(t == 0)}$beta),
                             kappa.true = {obs.7 %>% filter(x == 9)}$kappa,
                             intercept = unique(obs.7$intercept),
                             tau.alpha.true = unique(obs.7$tau.alpha),
                             tau.beta.true = unique(obs.7$tau.beta),
                             tau.kappa.true = unique(obs.7$tau.kappa),
                             tau.epsilon.true = unique(obs.7$tau.epsilon))

inlabru.synthetic.male.lung.lc.7<- inlabru.rw2.lc.2(obs.7, max_iter = 100)

plots.summaries.inlabru <- plot.inlabru.vs.underlying.synthetic.cancer(
  inlabru.synthetic.male.lung.lc.7,
  underlying.effects.7,
  path.to.storage = "Scripts/Synthetic\ data/Output/Figures/synthetic_male_lung_lc/v7",
  save=TRUE,
  cutoff_alpha = 75, cutoff_beta = 200, cutoff_kappa = 1000, cutoff_epsilon = 1000)



###   ----   v8   ----

# v8: use estimated lung cancer effects as input

inlabru.alpha <- res.lung.lc.m$summary.random$alpha$mean
inlabru.beta <- res.lung.lc.m$summary.random$beta$mean
inlabru.kappa <- res.lung.lc.m$summary.random$kappa$mean
inlabru.epsilon <- res.lung.lc.m$summary.random$epsilon$mean

alpha.8 <- inlabru.alpha - mean(inlabru.alpha)
beta.8 <- inlabru.beta - mean(inlabru.beta) + 1/18
kappa.8 <- inlabru.kappa - mean(inlabru.kappa) 
epsilon.8 <- inlabru.epsilon

set.seed(1)

obs.8 <- data.frame(x = lung.cancer.male$x, t = lung.cancer.male$t, xt = lung.cancer.male$xt, E = lung.cancer.male$E) %>%
  mutate(age.int = lung.cancer.male$age.int, year = lung.cancer.male$year) %>%
  mutate(x.c = x) %>%
  mutate(alpha = alpha.8[x + 1]) %>%
  mutate(beta = beta.8[x + 1]) %>%
  mutate(kappa = kappa.8[t + 1]) %>%
  mutate(intercept = inlabru.intercept) %>%
  mutate(epsilon = epsilon.8[xt + 1]) %>%
  mutate(eta = intercept + alpha + beta*kappa + epsilon) %>%
  mutate(Y = rpois(length(x), E*exp(eta))) %>%
  mutate(mr = Y/E) %>%
  mutate(tau.alpha = inlabru.tau.alpha) %>%
  mutate(tau.beta = inlabru.tau.beta) %>%
  mutate(tau.kappa = inlabru.tau.kappa) %>%
  mutate(tau.epsilon = inlabru.tau.epsilon)

#write.csv(obs.8, "Data/synthetic_male_lung_8.csv")
obs.8 <- read.csv("Data/synthetic_male_lung_8.csv")

underlying.effects.8 <- list(obs = obs.8, nx = 18, nt = 18,
                             alpha.true = {obs.8 %>% filter(t == 0)}$alpha,
                             beta.true = {obs.8 %>% filter(t == 0)}$beta,
                             kappa.true = {obs.8 %>% filter(x == 0)}$kappa,
                             intercept = unique(obs.8$intercept),
                             tau.alpha.true = unique(obs.8$tau.alpha),
                             tau.beta.true = unique(obs.8$tau.beta),
                             tau.kappa.true = unique(obs.8$tau.kappa),
                             tau.epsilon.true = unique(obs.8$tau.epsilon)
)

inlabru.synthetic.male.lung.lc.8<- inlabru.rw2.lc.2(obs.8, max_iter = 100)

source("Scripts/Synthetic data/plot_inlabru_vs_underlying.R")

plots.summaries.inlabru <- plot.inlabru.vs.underlying.synthetic.cancer(
  inlabru.synthetic.male.lung.lc.8,
  underlying.effects.8,
  path.to.storage = "Scripts/Synthetic\ data/Output/Figures/synthetic_male_lung_lc/v8",
  save=TRUE, cutoff_alpha = 50)
