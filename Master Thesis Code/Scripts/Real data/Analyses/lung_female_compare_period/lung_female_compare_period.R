# Compare different model choices for period effects in inlabru.

# Running full inlabru analysis on female lung cancer data

#   ----   Load libraries and set workspace   ----

setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code")

output.path = file.path("Scripts/Real\ Data/Analyses", "lung_female_compare_period")

library("tidyverse")
library("inlabru")
library("ggplot2")
library("INLA")
library("patchwork")
library("rstan")

#    ----    Load data   ----
load("Data/population-germany.Rda")
load("Data/lungCancer-germany.Rda")

#   ----   Format data   ----
female.lung.cancer <- cancer.data %>% select(c(age, year, female, t, age.int, x, x.c, xt, cohort, c, birth.year, female.t)) %>%
  mutate(Y = replace(female, year %in% 2011:2016, NA), E = female.t) %>%
  mutate(pred = "In data") %>% mutate(pred = replace(pred, year %in% 2011:2016, "Out of data")) %>%
  mutate(year = parse_integer(year))
observed <- female.lung.cancer

#   ----   kappa ~ rw2   ----

#   ----   Define inlabru analysis   ----

run.inlabru.rw2 <- function(obs, max_iter = 100){
  #' Defines the model in inlabru and runs it on the observed data in obs
  #' 
  #' @param obs <data.frame>: The observed data
  pc.prior <- list(prec = list(hyper = "pc.prec", param=c(1, 0.01)))
  loggamma.prior <- list(prec = list(hyper = "loggamma", param = c(1, 0.00005)))
  
  nx = length(unique(obs$x))
  nt = length(unique(obs$t))
  
  # constraints for the age effect beta
  A.beta = matrix(1, nrow = 1, ncol = nx)
  e.beta = 1
  
  comp = ~ -1 + 
    alpha(x, model = "rw1", hyper = pc.prior, constr = FALSE, scale.model = T) + 
    beta(x.c, model = "iid", hyper = loggamma.prior, extraconstr = list(A = A.beta, e = e.beta)) + 
    kappa(t, model = "rw2", hyper = pc.prior, constr = TRUE, scale.model = T) + 
    gamma(c, model = "rw1", hyper = pc.prior, constr = TRUE, scale.model = T) + 
    epsilon(xt, model = "iid", hyper = loggamma.prior, constr = FALSE)
  
  formula = Y ~ alpha + beta*kappa + gamma +  epsilon
  
  likelihood = like(formula = formula, family = "poisson", data = obs, E = obs$E)
  
  c.compute <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE, return.marginals.predictor = TRUE)  # control.compute
  
  result = bru(components = comp,
               likelihood, 
               options = list(verbose = F,
                              bru_verbose = 4, 
                              num.threads = "1:1",
                              control.compute = c.compute,
                              bru_max_iter=max_iter,
                              control.predictor = list(link = 1)
               ))
  
  return(result)
}

res.inlabru.rw2 <- run.inlabru.rw2(female.lung.cancer)

source("Scripts/Misc/palette.R")

#   ----   Generate samples for Y   ----

Y.samples.rw2 <- generate(res.inlabru.rw2, female.lung.cancer, ~ E*exp(alpha + beta*kappa + gamma + epsilon), n.samples = 10000)
Y.samples.df.rw2 <- data.frame(Y.samples.rw2) 

Y.inlabru.rw2 <- female.lung.cancer %>%
  mutate(Y.mean = apply(Y.samples.df.rw2, 1, mean)) %>%
  mutate(Y.0.025 = apply(Y.samples.df.rw2, 1, quantile, 0.025)) %>%
  mutate(Y.0.975 = apply(Y.samples.df.rw2, 1, quantile, 0.975)) %>%
  mutate(Y.sd = apply(Y.samples.df.rw2, 1, sd)) %>%
  mutate(DSS = ((female - Y.mean)/Y.sd)^2 + 2*log(Y.sd))

MDSS.all.rw2 <- mean(Y.inlabru.rw2$DSS)
MDSS.x.above.5.rw2 <- mean({Y.inlabru.rw2 %>% filter(x > 5)}$DSS)

MDSS.all.in.data.rw2 <- mean({Y.inlabru.rw2 %>% filter(year %in% 1999:2010)}$DSS)
MDSS.x.above.5.in.data.rw2 <- mean({Y.inlabru.rw2 %>% filter(x > 5) %>% filter(year %in% 1999:2010)}$DSS)

MDSS.all.out.data.rw2 <- mean({Y.inlabru.rw2 %>% filter(year %in% 2011:2016)}$DSS)
MDSS.x.above.5.out.data.rw2 <- mean({Y.inlabru.rw2 %>% filter(x > 5) %>% filter(year %in% 2011:2016)}$DSS)

write.table(list(MDSS.all = MDSS.all.rw2,
                 MDSS.x.above.5 = MDSS.x.above.5.rw2,
                 all.in.data = MDSS.all.in.data.rw2,
                 above.5.in.data = MDSS.x.above.5.in.data.rw2,
                 all.out.data = MDSS.all.out.data.rw2,
                 above.5.out.data = MDSS.x.above.5.out.data.rw2), file = file.path(output.path, "DSS_rw2.txt"))


p.Y.age.rw2 <- ggplot(Y.inlabru.rw2 %>% filter(year %in% 2011:2016)) + 
  geom_ribbon(aes(x = age.int, ymin = Y.0.025, ymax = Y.0.975, color = "Inlabru", fill = "Inlabru", shape = "Inlabru"), alpha = 0.2, size = 0.5) + 
  geom_point(aes(x = age.int, y = Y.mean, color = "Inlabru", fill = "Inlabru", shape = "Inlabru")) + 
  geom_point(aes(x = age.int, y = female, color = "Observed", fill = "Observed", shape = "Observed"), size = 2) + 
  facet_wrap(~ as.factor(year)) + 
  scale_color_manual(name="", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  scale_shape_manual(name = "", values = c(16,4)) + 
  theme_classic() + 
  labs(x = "Age", y = "")

ggsave("Y_by_age_rw2.pdf", p.Y.age.rw2, path = output.path, dpi = "retina", height = 5, width = 8)  

p.Y.age.in.data.rw2 <- ggplot(Y.inlabru.rw2 %>% filter(year %in% 1999:2010)) + 
  geom_ribbon(aes(x = age.int, ymin = Y.0.025, ymax = Y.0.975, color = "Inlabru", fill = "Inlabru", shape = "Inlabru"), alpha = 0.2, size = 0.5) + 
  geom_point(aes(x = age.int, y = Y.mean, color = "Inlabru", fill = "Inlabru", shape = "Inlabru")) + 
  geom_point(aes(x = age.int, y = female, color = "Observed", fill = "Observed", shape = "Observed"), size = 2) + 
  facet_wrap(~ as.factor(year)) + 
  scale_color_manual(name="", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  scale_shape_manual(name = "", values = c(16,4)) + 
  theme_classic() + 
  labs(x = "Age", y = "")

ggsave("Y_by_age_in_data_rw2.pdf", p.Y.age.in.data.rw2, path = output.path, dpi = "retina", height = 5, width = 8)  


p.Y.year.rw2 <- ggplot(Y.inlabru.rw2 %>% filter(age.int >= 30)) + 
  geom_ribbon(aes(x = year, ymin = Y.0.025, ymax = Y.0.975, color = "Inlabru", fill = "Inlabru", shape = "Inlabru"), alpha = 0.2, size = 0.5) + 
  geom_point(aes(x = year, y = Y.mean, color = "Inlabru", fill = "Inlabru", shape = "Inlabru")) + 
  geom_point(aes(x = year, y = female, color = "Observed", fill = "Observed", shape = "Observed"), size = 2) + 
  geom_vline(aes(xintercept = 2011, color="Predicted period", fill = "Predicted period", shape = "Predicted period")) + 
  facet_wrap(~ as.factor(age)) + 
  scale_color_manual(name="", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  scale_shape_manual(name = "", values = c(16,4,3)) + 
  theme_classic() + 
  labs(x = "Year", y = "")

ggsave("Y_by_year_rw2.pdf", p.Y.year.rw2, path = output.path, dpi = "retina", height = 5, width = 8)  

#   ----   Plot random effects   ----

p.alpha.rw2 <- ggplot(data.frame(res.inlabru.rw2$summary.random$alpha), aes(x = ID)) + 
  geom_ribbon(aes(ymin = X0.025quant, ymax = X0.975quant), color = palette[1], fill = palette[1], alpha = 0.3) + 
  geom_point(aes(y = mean), color = palette[1]) + 
  theme_classic() + 
  labs(x = "x", y = "", title = "Alpha")

p.beta.rw2 <- ggplot(data.frame(res.inlabru.rw2$summary.random$beta), aes(x = ID)) + 
  geom_errorbar(aes(ymin = X0.025quant, ymax = X0.975quant), color = palette[1], fill = palette[1], alpha = 0.7) + 
  geom_point(aes(y = mean), color = palette[1]) + 
  theme_classic() + 
  labs(x = "x", y = "", title = "Beta")

p.kappa.rw2 <- ggplot(data.frame(res.inlabru.rw2$summary.random$kappa), aes(x = ID)) + 
  geom_vline(aes(xintercept = 12), color = palette[3]) + 
  geom_ribbon(aes(ymin = X0.025quant, ymax = X0.975quant), color = palette[1], fill = palette[1], alpha = 0.3) + 
  geom_point(aes(y = mean), color = palette[1]) + 
  theme_classic() + 
  #scale_color_manual(name="", values = c(palette[3])) + 
  labs(x = "t", y = "", title = "Kappa")

p.gamma.rw2 <- ggplot(data.frame(res.inlabru.rw2$summary.random$gamma), aes(x = ID)) + 
  geom_vline(aes(xintercept = 11), color = palette[3]) + 
  geom_vline(aes(xintercept  = 97), color = palette[4]) + 
  geom_ribbon(aes(ymin = X0.025quant, ymax = X0.975quant), color = palette[1], fill = palette[1], alpha = 0.3) + 
  geom_point(aes(y = mean), color = palette[1]) + 
  theme_classic() + 
  labs(x = "c", y = "", title = "Gamma")

p.random.rw2 <- (p.alpha.rw2 | p.beta.rw2)/(p.kappa.rw2 | p.gamma.rw2) + plot_layout(guides = "collect")
ggsave("random_rw2.pdf", p.random.rw2, path = output.path, dpi = "retina", height = 5, width = 8)

#   ----   Heatmap of epsilon   ----

p.epsilon.rw2 <- ggplot(data.frame(res.inlabru.rw2$summary.random$epsilon) %>% left_join(observed %>% select(c("x","t","xt")), by = c("ID" = "xt"))) + 
  geom_tile(aes(x = x, y = t, fill = mean)) + 
  scale_fill_gradient(low = palette[2], high = palette[6]) + 
  theme_classic() + 
  labs(x = "x", y = "t")

ggsave("epsilon_rw2.pdf", p.epsilon.rw2, path = output.path, dpi = "retina", height = 3, width = 4)

#   ----   Plot precisions of random effects   ----

p.alpha.prec.rw2 <- ggplot(data.frame(res.inlabru.rw2$marginals.hyperpar$`Precision for alpha`) %>% filter(x < 4)) + 
  geom_area(aes(x=x, y=y), color = palette[1], fill = palette[1], alpha = 0.5) + 
  theme_classic() + 
  labs(title = "Alpha", x="", y = "")

p.beta.prec.rw2 <- ggplot(data.frame(res.inlabru.rw2$marginals.hyperpar$`Precision for beta`) %>% filter(x < 750)) + 
  geom_area(aes(x=x, y=y), color = palette[1], fill = palette[1], alpha = 0.5) + 
  theme_classic() + 
  labs(title = "Beta", x="", y = "")

p.kappa.prec.rw2 <- ggplot(data.frame(res.inlabru.rw2$marginals.hyperpar$`Precision for kappa`) %>% filter(x < 600)) + 
  geom_area(aes(x=x, y=y), color = palette[1], fill = palette[1], alpha = 0.5) + 
  theme_classic() + 
  labs(title = "Kappa", x="", y = "")

p.gamma.prec.rw2 <- ggplot(data.frame(res.inlabru.rw2$marginals.hyperpar$`Precision for gamma`) %>% filter(x < 100)) + 
  geom_area(aes(x=x, y=y), color = palette[1], fill = palette[1], alpha = 0.5) + 
  theme_classic() + 
  labs(title = "Gamma", x="", y = "")

p.epsilon.prec.rw2 <- ggplot(data.frame(res.inlabru.rw2$marginals.hyperpar$`Precision for epsilon`) %>% filter(x < 150000)) + 
  geom_area(aes(x=x, y=y), color = palette[1], fill = palette[1], alpha = 0.5) + 
  theme_classic() + 
  labs(title = "Epsilon", x="", y = "")

p.prec.rw2 <- (p.alpha.prec.rw2 | p.beta.prec.rw2)/(p.kappa.prec.rw2 | p.gamma.prec.rw2 |p.epsilon.prec.rw2) + plot_layout(guides = "collect")
ggsave("precisions_rw2.pdf", p.prec.rw2, path = output.path, dpi = "retina", height = 5, width = 8)

#   ----   kappa ~ decomposed random walk with drift   ----

#   ----   Define inlabru analysis   ----

run.inlabru.drift <- function(obs, max_iter = 100){
  #' Defines the model in inlabru and runs it on the observed data in obs
  #' 
  #' @param obs <data.frame>: The observed data
  pc.prior <- list(prec = list(hyper = "pc.prec", param=c(1, 0.01)))
  loggamma.prior <- list(prec = list(hyper = "loggamma", param = c(1, 0.00005)))
  
  nx = length(unique(obs$x))
  nt = length(unique(obs$t))
  
  # constraints for the age effect beta
  A.beta = matrix(1, nrow = 1, ncol = nx)
  e.beta = 1
  
  comp = ~ -1 + 
    alpha(x, model = "rw1", hyper = pc.prior, constr = FALSE, scale.model = T) + 
    beta(x.c, model = "iid", hyper = loggamma.prior, extraconstr = list(A = A.beta, e = e.beta)) + 
    phi(t, model = "linear", prec.linear = 1, mean.linear = 0) + 
    kappa(t, model = "rw1", hyper = pc.prior, constr = TRUE, scale.model = T) + 
    gamma(c, model = "rw1", hyper = pc.prior, constr = TRUE, scale.model = T) + 
    epsilon(xt, model = "iid", hyper = loggamma.prior, constr = FALSE)
  
  formula = Y ~ alpha + beta*phi + beta*kappa + gamma +  epsilon
  
  likelihood = like(formula = formula, family = "poisson", data = obs, E = obs$E)
  
  c.compute <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE, return.marginals.predictor = TRUE)  # control.compute
  
  result = bru(components = comp,
               likelihood, 
               options = list(verbose = F,
                              bru_verbose = 4, 
                              num.threads = "1:1",
                              control.compute = c.compute,
                              bru_max_iter=max_iter,
                              control.predictor = list(link = 1)
               ))
  
  return(result)
}

res.inlabru.drift <- run.inlabru.drift(female.lung.cancer)

#   ----   Generate samples for Y   ----

Y.samples.drift <- generate(res.inlabru.drift, female.lung.cancer, ~ E*exp(alpha + beta*phi + beta*kappa + gamma + epsilon), n.samples = 10000)
Y.samples.df.drift <- data.frame(Y.samples.drift) 

Y.inlabru.drift <- female.lung.cancer %>%
  mutate(Y.mean = apply(Y.samples.df.drift, 1, mean)) %>%
  mutate(Y.0.025 = apply(Y.samples.df.drift, 1, quantile, 0.025)) %>%
  mutate(Y.0.975 = apply(Y.samples.df.drift, 1, quantile, 0.975)) %>%
  mutate(Y.sd = apply(Y.samples.df.drift, 1, sd)) %>%
  mutate(DSS = ((female - Y.mean)/Y.sd)^2 + 2*log(Y.sd))

MDSS.all.drift <- mean(Y.inlabru.drift$DSS)
MDSS.x.above.5.drift <- mean({Y.inlabru.drift %>% filter(x > 5)}$DSS)

MDSS.all.in.data.drift <- mean({Y.inlabru.drift %>% filter(year %in% 1999:2010)}$DSS)
MDSS.x.above.5.in.data.drift <- mean({Y.inlabru.drift %>% filter(x > 5) %>% filter(year %in% 1999:2010)}$DSS)

MDSS.all.out.data.drift <- mean({Y.inlabru.drift %>% filter(year %in% 2011:2016)}$DSS)
MDSS.x.above.5.out.data.drift <- mean({Y.inlabru.drift %>% filter(x > 5) %>% filter(year %in% 2011:2016)}$DSS)

write.table(list(MDSS.all = MDSS.all.drift,
                 MDSS.x.above.5 = MDSS.x.above.5.drift,
                 all.in.data = MDSS.all.in.data.drift,
                 above.5.in.data = MDSS.x.above.5.in.data.drift,
                 all.out.data = MDSS.all.out.data.drift,
                 above.5.out.data = MDSS.x.above.5.out.data.drift), file = file.path(output.path, "DSS_drift.txt"))


p.Y.age.drift <- ggplot(Y.inlabru.drift %>% filter(year %in% 2011:2016)) + 
  geom_ribbon(aes(x = age.int, ymin = Y.0.025, ymax = Y.0.975, color = "Inlabru", fill = "Inlabru", shape = "Inlabru"), alpha = 0.2, size = 0.5) + 
  geom_point(aes(x = age.int, y = Y.mean, color = "Inlabru", fill = "Inlabru", shape = "Inlabru")) + 
  geom_point(aes(x = age.int, y = female, color = "Observed", fill = "Observed", shape = "Observed"), size = 2) + 
  facet_wrap(~ as.factor(year)) + 
  scale_color_manual(name="", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  scale_shape_manual(name = "", values = c(16,4)) + 
  theme_classic() + 
  labs(x = "Age", y = "")

ggsave("Y_by_age_drift.pdf", p.Y.age.drift, path = output.path, dpi = "retina", height = 5, width = 8)  

p.Y.age.in.data.drift <- ggplot(Y.inlabru.drift %>% filter(year %in% 1999:2010)) + 
  geom_ribbon(aes(x = age.int, ymin = Y.0.025, ymax = Y.0.975, color = "Inlabru", fill = "Inlabru", shape = "Inlabru"), alpha = 0.2, size = 0.5) + 
  geom_point(aes(x = age.int, y = Y.mean, color = "Inlabru", fill = "Inlabru", shape = "Inlabru")) + 
  geom_point(aes(x = age.int, y = female, color = "Observed", fill = "Observed", shape = "Observed"), size = 2) + 
  facet_wrap(~ as.factor(year)) + 
  scale_color_manual(name="", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  scale_shape_manual(name = "", values = c(16,4)) + 
  theme_classic() + 
  labs(x = "Age", y = "")

ggsave("Y_by_age_in_data_drift.pdf", p.Y.age.in.data.drift, path = output.path, dpi = "retina", height = 5, width = 8)  


p.Y.year.drift <- ggplot(Y.inlabru.drift %>% filter(age.int >= 30)) + 
  geom_ribbon(aes(x = year, ymin = Y.0.025, ymax = Y.0.975, color = "Inlabru", fill = "Inlabru", shape = "Inlabru"), alpha = 0.2, size = 0.5) + 
  geom_point(aes(x = year, y = Y.mean, color = "Inlabru", fill = "Inlabru", shape = "Inlabru")) + 
  geom_point(aes(x = year, y = female, color = "Observed", fill = "Observed", shape = "Observed"), size = 2) + 
  geom_vline(aes(xintercept = 2011, color="Predicted period", fill = "Predicted period", shape = "Predicted period")) + 
  facet_wrap(~ as.factor(age)) + 
  scale_color_manual(name="", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  scale_shape_manual(name = "", values = c(16,4,3)) + 
  theme_classic() + 
  labs(x = "Year", y = "")

ggsave("Y_by_year_drift.pdf", p.Y.year.drift, path = output.path, dpi = "retina", height = 5, width = 8)  

#   ----   Plot random effects   ----
period.samples.drift <- generate(res.inlabru.drift, female.lung.cancer, ~ phi + kappa, n.samples = 10000)
period.samples.drift.df <- data.frame(period.samples.drift)

period.drift <- data.frame(t = res.inlabru.drift$summary.random$kappa$ID,
                           period.mean = apply(period.samples.drift.df, 1, mean),
                           period.0.025 = apply(period.samples.drift.df, 1, quantile, 0.025),
                           period.0.975 = apply(period.samples.drift.df, 1, quantile, 0.975))

p.alpha.drift <- ggplot(data.frame(res.inlabru.drift$summary.random$alpha), aes(x = ID)) + 
  geom_ribbon(aes(ymin = X0.025quant, ymax = X0.975quant), color = palette[1], fill = palette[1], alpha = 0.3) + 
  geom_point(aes(y = mean), color = palette[1]) + 
  theme_classic() + 
  labs(x = "x", y = "", title = "Alpha")

p.beta.drift <- ggplot(data.frame(res.inlabru.drift$summary.random$beta), aes(x = ID)) + 
  geom_errorbar(aes(ymin = X0.025quant, ymax = X0.975quant), color = palette[1], fill = palette[1], alpha = 0.7) + 
  geom_point(aes(y = mean), color = palette[1]) + 
  theme_classic() + 
  labs(x = "x", y = "", title = "Beta")

# p.kappa.drift <- ggplot(data.frame(res.inlabru.drift$summary.random$kappa), aes(x = ID)) + 
#   geom_vline(aes(xintercept = 12), color = palette[3]) + 
#   geom_ribbon(aes(ymin = X0.025quant, ymax = X0.975quant), color = palette[1], fill = palette[1], alpha = 0.3) + 
#   geom_point(aes(y = mean), color = palette[1]) + 
#   theme_classic() + 
#   #scale_color_manual(name="", values = c(palette[3])) + 
#   labs(x = "t", y = "", title = "Kappa")

p.period.drift <- ggplot(period.drift) + 
  geom_ribbon(aes(x = t, ymin = period.0.025, ymax = period.0.975), color = palette[1], fill = palette[1], alpha = 0.3) + 
  geom_point(aes(x = t, y = period.mean), color = palette[1]) + 
  geom_vline(aes(xintercept = 12), color = palette[3]) + 
  theme_classic() + 
  labs(x = "t", y= "", title = "Phi + Kappa")

p.gamma.drift <- ggplot(data.frame(res.inlabru.drift$summary.random$gamma), aes(x = ID)) + 
  geom_vline(aes(xintercept = 11), color = palette[3]) + 
  geom_vline(aes(xintercept  = 97), color = palette[4]) + 
  geom_ribbon(aes(ymin = X0.025quant, ymax = X0.975quant), color = palette[1], fill = palette[1], alpha = 0.3) + 
  geom_point(aes(y = mean), color = palette[1]) + 
  theme_classic() + 
  labs(x = "c", y = "", title = "Gamma")


p.random.drift <- (p.alpha.drift | p.beta.drift)/(p.period.drift | p.gamma.drift) + plot_layout(guides = "collect")
ggsave("random_drift.pdf", p.random.drift, path = output.path, dpi = "retina", height = 5, width = 8)

#   ----   Heatmap of epsilon   ----

p.epsilon.drift <- ggplot(data.frame(res.inlabru.drift$summary.random$epsilon) %>% left_join(observed %>% select(c("x","t","xt")), by = c("ID" = "xt"))) + 
  geom_tile(aes(x = x, y = t, fill = mean)) + 
  scale_fill_gradient(low = palette[2], high = palette[6]) + 
  theme_classic() + 
  labs(x = "x", y = "t")

ggsave("epsilon_drift.pdf", p.epsilon.drift, path = output.path, dpi = "retina", height = 3, width = 4)

#   ----   Plot precisions of random effects   ----

p.alpha.prec.drift <- ggplot(data.frame(res.inlabru.drift$marginals.hyperpar$`Precision for alpha`) %>% filter(x < 4)) + 
  geom_area(aes(x=x, y=y), color = palette[1], fill = palette[1], alpha = 0.5) + 
  theme_classic() + 
  labs(title = "Alpha", x="", y = "")

p.beta.prec.drift <- ggplot(data.frame(res.inlabru.drift$marginals.hyperpar$`Precision for beta`) %>% filter(x < 750)) + 
  geom_area(aes(x=x, y=y), color = palette[1], fill = palette[1], alpha = 0.5) + 
  theme_classic() + 
  labs(title = "Beta", x="", y = "")

p.kappa.prec.drift <- ggplot(data.frame(res.inlabru.drift$marginals.hyperpar$`Precision for kappa`) %>% filter(x < 600)) + 
  geom_area(aes(x=x, y=y), color = palette[1], fill = palette[1], alpha = 0.5) + 
  theme_classic() + 
  labs(title = "Kappa", x="", y = "")

p.gamma.prec.drift <- ggplot(data.frame(res.inlabru.drift$marginals.hyperpar$`Precision for gamma`) %>% filter(x < 100)) + 
  geom_area(aes(x=x, y=y), color = palette[1], fill = palette[1], alpha = 0.5) + 
  theme_classic() + 
  labs(title = "Gamma", x="", y = "")

p.epsilon.prec.drift <- ggplot(data.frame(res.inlabru.drift$marginals.hyperpar$`Precision for epsilon`) %>% filter(x < 150000)) + 
  geom_area(aes(x=x, y=y), color = palette[1], fill = palette[1], alpha = 0.5) + 
  theme_classic() + 
  labs(title = "Epsilon", x="", y = "")

p.prec.drift <- (p.alpha.prec.drift | p.beta.prec.drift)/(p.kappa.prec.drift | p.gamma.prec.drift |p.epsilon.prec.drift) + plot_layout(guides = "collect")
ggsave("precisions_drift.pdf", p.prec.drift, path = output.path, dpi = "retina", height = 5, width = 8)