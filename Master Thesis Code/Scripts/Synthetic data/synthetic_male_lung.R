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

# sample tree different versions - set seeds accordingly
set.seed(123)
alpha.1 <- rw1(inlabru.tau.alpha, 18); alpha.1 = alpha.1 - mean(alpha.1)
beta.1 <- rnorm(18, sd = sqrt(1/inlabru.tau.beta)); beta.1 = beta.1 - mean(beta.1) + 1/18
kappa.1 <- rw2(inlabru.tau.kappa, 18); kappa.1 = kappa.1 - mean(kappa.1)
epsilon.1 <- rnorm(18*18, sd = sqrt(1/inlabru.tau.epsilon))

obs.1 <- data.frame(x = lung.cancer.male$x, t = lung.cancer.male$t, xt = lung.cancer.male$xt, E = lung.cancer.male$E) %>%
  mutate(age.int = lung.cancer.male$age.int, year = lung.cancer.male$year) %>%
  mutate(x.c = x) %>%
  mutate(alpha = alpha.1[x + 1]) %>%
  mutate(beta = beta.1[x + 1]) %>%
  mutate(kappa = kappa.1[t + 1]) %>%
  mutate(intercept = inlabru.intercept) %>%
  mutate(epsilon = epsilon.1[xt + 1]) %>%
  mutate(eta = intercept + alpha + beta*kappa + epsilon) %>%
  mutate(Y = rpois(length(x), E*exp(eta))) %>%
  mutate(mr = Y/E) %>%
  mutate(tau.alpha = inlabru.tau.alpha) %>%
  mutate(tau.beta = inlabru.tau.beta) %>%
  mutate(tau.kappa = inlabru.tau.kappa) %>%
  mutate(tau.epsilon = inlabru.tau.epsilon)

#write.csv(obs.1, "Data/synthetic_male_lung_1.csv", row.names = FALSE)
obs.1 <- read.csv("Data/synthetic_male_lung_1.csv")
underlying.effects.1 <- list(obs = obs.1, nx = 18, nt = 18,
                             alpha.true = {obs.1 %>% filter(t == 0)}$alpha,
                             beta.true = {obs.1 %>% filter(t == 0)}$beta,
                             kappa.true = {obs.1 %>% filter(x == 0)}$kappa,
                             intercept = unique(obs.1$intercept))

source("Scripts/Functions/inlabru_analyses.R")

inlabru.synthetic.male.lung.lc <- inlabru.rw2.lc.2(obs.1, max_iter = 100)

source("Scripts/Synthetic data/plot_inlabru_vs_underlying.R")

plots.summaries.inlabru <- plot.inlabru.vs.underlying.synthetic.cancer(
  inlabru.synthetic.male.lung.lc,
  underlying.effects.1,
  path.to.storage = "Scripts/Synthetic\ data/Output/Figures/synthetic_male_lung_lc/v1",
  save=TRUE)

set.seed(124)
alpha.2 <- rw1(inlabru.tau.alpha, 18); alpha.2 = alpha.2 - mean(alpha.2)
beta.2 <- rnorm(18, sd = sqrt(1/inlabru.tau.beta)); beta.2 = beta.2 - mean(beta.2) + 1/18
kappa.2 <- rw2(inlabru.tau.kappa, 18); kappa.2 = kappa.2 - mean(kappa.2)
epsilon.2 <- rnorm(18*18, sd = sqrt(1/inlabru.tau.epsilon))

obs.2 <- data.frame(x = lung.cancer.male$x, t = lung.cancer.male$t, xt = lung.cancer.male$xt, E = lung.cancer.male$E) %>%
  mutate(age.int = lung.cancer.male$age.int, year = lung.cancer.male$year) %>%
  mutate(x.c = x) %>%
  mutate(alpha = alpha.2[x + 1]) %>%
  mutate(beta = beta.2[x + 1]) %>%
  mutate(kappa = kappa.2[t + 1]) %>%
  mutate(intercept = inlabru.intercept) %>%
  mutate(epsilon = epsilon.2[xt + 1]) %>%
  mutate(eta = intercept + alpha + beta*kappa + epsilon) %>%
  mutate(Y = rpois(length(x), E*exp(eta))) %>%
  mutate(mr = Y/E) %>%
  mutate(tau.alpha = inlabru.tau.alpha) %>%
  mutate(tau.beta = inlabru.tau.beta) %>%
  mutate(tau.kappa = inlabru.tau.kappa) %>%
  mutate(tau.epsilon = inlabru.tau.epsilon)

#write.csv(obs.2, "Data/synthetic_male_lung_2.csv")
#obs.2 <- read.csv("Data/synthetic_male_lung_2.csv")

underlying.effects.2 <- list(obs = obs.1, nx = 18, nt = 18,
                             alpha.true = {obs.2 %>% filter(t == 0)}$alpha,
                             beta.true = {obs.2 %>% filter(t == 0)}$beta,
                             kappa.true = {obs.2 %>% filter(x == 0)}$kappa,
                             intercept = unique(obs.2$intercept))

source("Scripts/Functions/inlabru_analyses.R")

inlabru.synthetic.male.lung.lc.2 <- inlabru.rw2.lc.2(obs.2, max_iter = 100)

source("Scripts/Synthetic data/plot_inlabru_vs_underlying.R")

plots.summaries.inlabru <- plot.inlabru.vs.underlying.synthetic.cancer(
  inlabru.synthetic.male.lung.lc.2,
  underlying.effects.2,
  path.to.storage = "Scripts/Synthetic\ data/Output/Figures/synthetic_male_lung_lc/v2",
  save=TRUE)

set.seed(125)
alpha.3 <- rw1(inlabru.tau.alpha, 18); alpha.3 = alpha.3 - mean(alpha.3)
beta.3 <- rnorm(18, sd = sqrt(1/inlabru.tau.beta)); beta.3 = beta.3 - mean(beta.3) + 1/18
kappa.3 <- rw2(inlabru.tau.kappa, 18); kappa.3 = kappa.3 - mean(kappa.3)
epsilon.3 <- rnorm(18*18, sd = sqrt(1/inlabru.tau.epsilon))

obs.3 <- data.frame(x = lung.cancer.male$x, t = lung.cancer.male$t, xt = lung.cancer.male$xt, E = lung.cancer.male$E) %>%
  mutate(age.int = lung.cancer.male$age.int, year = lung.cancer.male$year) %>%
  mutate(x.c = x) %>%
  mutate(alpha = alpha.3[x + 1]) %>%
  mutate(beta = beta.3[x + 1]) %>%
  mutate(kappa = kappa.3[t + 1]) %>%
  mutate(intercept = inlabru.intercept) %>%
  mutate(epsilon = epsilon.3[xt + 1]) %>%
  mutate(eta = intercept + alpha + beta*kappa + epsilon) %>%
  mutate(Y = rpois(length(x), E*exp(eta))) %>%
  mutate(mr = Y/E) %>%
  mutate(tau.alpha = inlabru.tau.alpha) %>%
  mutate(tau.beta = inlabru.tau.beta) %>%
  mutate(tau.kappa = inlabru.tau.kappa) %>%
  mutate(tau.epsilon = inlabru.tau.epsilon)

write.csv(obs.3, "Data/synthetic_male_lung_3.csv")

# attempt to find a configuration that is of somewhat comparable size to real data
#set.seed(4) # OK, try to find better
# set.seed(7) # Ok for alpha and eta, worse for kappa
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
ggplot(obs.4) + geom_line(aes(x = x, y = alpha))
ggplot(obs.4) + geom_line(aes(x = x, y = beta))
ggplot(obs.4) + geom_line(aes(x = t, y = kappa))

write.csv(obs.4, "Data/synthetic_male_lung_4.csv")

underlying.effects.4 <- list(obs = obs.4, nx = 18, nt = 18,
                             alpha.true = {obs.4 %>% filter(t == 0)}$alpha,
                             beta.true = {obs.4 %>% filter(t == 0)}$beta,
                             kappa.true = {obs.4 %>% filter(x == 0)}$kappa,
                             intercept = unique(obs.4$intercept))

source("Scripts/Functions/inlabru_analyses.R")

inlabru.synthetic.male.lung.lc.4 <- inlabru.rw2.lc.2(obs.4, max_iter = 100)

source("Scripts/Synthetic data/plot_inlabru_vs_underlying.R")

plots.summaries.inlabru <- plot.inlabru.vs.underlying.synthetic.cancer(
  inlabru.synthetic.male.lung.lc.4,
  underlying.effects.4,
  path.to.storage = "Scripts/Synthetic\ data/Output/Figures/synthetic_male_lung_lc/v4",
  save=TRUE)

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
  

# attempt with stan hyperpars instead:

stan.tau.alpha = 2
stan.tau.beta = 1000
stan.tau.kappa = 192

#set.seed(14) # at least hight enough mortality rates, very little period effect
set.seed(17) # could work, somewhat ok mr


alpha.6 <- rw1(stan.tau.alpha - 1, 18); alpha.6 = alpha.6 - mean(alpha.6)
beta.6 <- rnorm(18, sd = sqrt(1/stan.tau.beta)); beta.6 = beta.6 - mean(beta.6) + 1/18
kappa.6 <- rw2(stan.tau.kappa - 150, 18); kappa.6 = kappa.6 - mean(kappa.6)
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
  mutate(tau.alpha = inlabru.tau.alpha) %>%
  mutate(tau.beta = inlabru.tau.beta) %>%
  mutate(tau.kappa = inlabru.tau.kappa) %>%
  mutate(tau.epsilon = inlabru.tau.epsilon)

ggplot(obs.6) + geom_line(aes(x = x, y = mr, color = year))
ggplot(obs.6) + geom_line(aes(x = x, y = eta, color = year))
ggplot(obs.6) + geom_line(aes(x = x, y = alpha))
ggplot(obs.6) + geom_line(aes(x = x, y = beta))
ggplot(obs.6) + geom_line(aes(x = t, y = kappa))

write.csv(obs.6, "Data/synthetic_male_lung_6.csv")

underlying.effects.6 <- list(obs = obs.6, nx = 18, nt = 18,
                             alpha.true = {obs.6 %>% filter(t == 0)}$alpha,
                             beta.true = {obs.6 %>% filter(t == 0)}$beta,
                             kappa.true = {obs.6 %>% filter(x == 0)}$kappa,
                             intercept = unique(obs.6$intercept))

source("Scripts/Functions/inlabru_analyses.R")

inlabru.synthetic.male.lung.lc.6<- inlabru.rw2.lc.2(obs.6, max_iter = 100)

source("Scripts/Synthetic data/plot_inlabru_vs_underlying.R")

plots.summaries.inlabru <- plot.inlabru.vs.underlying.synthetic.cancer(
  inlabru.synthetic.male.lung.lc.6,
  underlying.effects.6,
  path.to.storage = "Scripts/Synthetic\ data/Output/Figures/synthetic_male_lung_lc/v6",
  save=TRUE)


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

ggplot(obs.7) + geom_line(aes(x = x, y = mr, color = year))
ggplot(obs.7) + geom_line(aes(x = x, y = eta, color = year))
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

