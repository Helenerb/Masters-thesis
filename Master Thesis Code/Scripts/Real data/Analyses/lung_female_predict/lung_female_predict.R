# Running full inlabru analysis on female lung cancer data

#   ----   Load libraries and set workspace   ----

setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code")

output.path = file.path("Scripts/Real\ Data/Analyses", "lung_female_predict")

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

#   ----   Define inlabru analysis   ----

run.inlabru <- function(obs, max_iter = 100){
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

res.inlabru <- run.inlabru(female.lung.cancer)

observed <- female.lung.cancer

source("Scripts/Misc/palette.R")

# look at fitted values:
ggplot(data.frame(mean = res.inlabru$summary.fitted.values$mean[1:324],
                  x = female.lung.cancer$x,
                  year = female.lung.cancer$year,
                  true = female.lung.cancer$Y/female.lung.cancer$E)) + 
  geom_point(aes(x = x, y = mean, color = "Inlabru")) + 
  geom_point(aes(x = x, y = true, color = "True")) + 
  facet_wrap(~as.factor(year)) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette)  +
  labs(title = "Fitted values")

#   ----   Generate samples for Y   ----

Y.samples <- generate(res.inlabru, female.lung.cancer, ~ E*exp(alpha + beta*kappa + gamma + epsilon), n.samples = 1000)
Y.samples.df <- data.frame(Y.samples) 

Y.inlabru <- female.lung.cancer %>%
  mutate(Y.mean = apply(Y.samples.df, 1, mean)) %>%
  mutate(Y.0.025 = apply(Y.samples.df, 1, quantile, 0.025)) %>%
  mutate(Y.0.975 = apply(Y.samples.df, 1, quantile, 0.975)) %>%
  mutate(Y.sd = apply(Y.samples.df, 1, sd)) %>%
  mutate(DSS = ((female - Y.mean)/Y.sd)^2 + 2*log(Y.sd))

MDSS.all <- mean(Y.inlabru$DSS)
MDSS.x.above.5 <- mean({Y.inlabru %>% filter(x > 5)}$DSS)

MDSS.all.in.data <- mean({Y.inlabru %>% filter(year %in% 1999:2010)}$DSS)
MDSS.x.above.5.in.data <- mean({Y.inlabru %>% filter(x > 5) %>% filter(year %in% 1999:2010)}$DSS)

MDSS.all.out.data <- mean({Y.inlabru %>% filter(year %in% 2011:2016)}$DSS)
MDSS.x.above.5.out.data <- mean({Y.inlabru %>% filter(x > 5) %>% filter(year %in% 2011:2016)}$DSS)

write.table(list(MDSS.all = MDSS.all,
                 MDSS.x.above.5 = MDSS.x.above.5,
                 all.in.data = MDSS.all.in.data,
                 above.5.in.data = MDSS.x.above.5.in.data,
                 all.out.data = MDSS.all.out.data,
                 above.5.out.data = MDSS.x.above.5.out.data), file = file.path(output.path, "DSS.txt"))


p.Y.age <- ggplot(Y.inlabru %>% filter(year %in% 2011:2016)) + 
  geom_ribbon(aes(x = age.int, ymin = Y.0.025, ymax = Y.0.975, color = "Inlabru", fill = "Inlabru", shape = "Inlabru"), alpha = 0.2, size = 0.5) + 
  geom_point(aes(x = age.int, y = Y.mean, color = "Inlabru", fill = "Inlabru", shape = "Inlabru")) + 
  geom_point(aes(x = age.int, y = female, color = "Observed", fill = "Observed", shape = "Observed"), size = 2) + 
  facet_wrap(~ as.factor(year)) + 
  scale_color_manual(name="", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  scale_shape_manual(name = "", values = c(16,4)) + 
  theme_classic() + 
  labs(x = "Age", y = "")

ggsave("Y_by_age.pdf", p.Y.age, path = output.path, dpi = "retina", height = 5, width = 8)  

p.Y.age.in.data <- ggplot(Y.inlabru %>% filter(year %in% 1999:2010)) + 
  geom_ribbon(aes(x = age.int, ymin = Y.0.025, ymax = Y.0.975, color = "Inlabru", fill = "Inlabru", shape = "Inlabru"), alpha = 0.2, size = 0.5) + 
  geom_point(aes(x = age.int, y = Y.mean, color = "Inlabru", fill = "Inlabru", shape = "Inlabru")) + 
  geom_point(aes(x = age.int, y = female, color = "Observed", fill = "Observed", shape = "Observed"), size = 2) + 
  facet_wrap(~ as.factor(year)) + 
  scale_color_manual(name="", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  scale_shape_manual(name = "", values = c(16,4)) + 
  theme_classic() + 
  labs(x = "Age", y = "")

ggsave("Y_by_age_in_data.pdf", p.Y.age.in.data, path = output.path, dpi = "retina", height = 5, width = 8)  


p.Y.year <- ggplot(Y.inlabru %>% filter(age.int >= 30)) + 
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

ggsave("Y_by_year.pdf", p.Y.year, path = output.path, dpi = "retina", height = 5, width = 8)  

#   ----   Plot random effects   ----

p.alpha <- ggplot(data.frame(res.inlabru$summary.random$alpha), aes(x = ID)) + 
  geom_ribbon(aes(ymin = X0.025quant, ymax = X0.975quant), color = palette[1], fill = palette[1], alpha = 0.3) + 
  geom_point(aes(y = mean), color = palette[1]) + 
  theme_classic() + 
  labs(x = "x", y = "", title = "Alpha")

p.beta <- ggplot(data.frame(res.inlabru$summary.random$beta), aes(x = ID)) + 
  geom_errorbar(aes(ymin = X0.025quant, ymax = X0.975quant), color = palette[1], fill = palette[1], alpha = 0.7) + 
  geom_point(aes(y = mean), color = palette[1]) + 
  theme_classic() + 
  labs(x = "x", y = "", title = "Beta")

p.kappa <- ggplot(data.frame(res.inlabru$summary.random$kappa), aes(x = ID)) + 
  geom_vline(aes(xintercept = 12), color = palette[3]) + 
  geom_ribbon(aes(ymin = X0.025quant, ymax = X0.975quant), color = palette[1], fill = palette[1], alpha = 0.3) + 
  geom_point(aes(y = mean), color = palette[1]) + 
  theme_classic() + 
  #scale_color_manual(name="", values = c(palette[3])) + 
  labs(x = "t", y = "", title = "Kappa")

observed.cohorts <- observed %>% select(c(c, year)) %>% mutate(pred.int = if_else(year %in% 2011:2016, 0, 1)) %>%
  group_by(c) %>%
  summarise(avg = mean(pred.int))

p.gamma <- ggplot(data.frame(res.inlabru$summary.random$gamma), aes(x = ID)) + 
  geom_vline(aes(xintercept = 11), color = palette[3]) + 
  geom_vline(aes(xintercept  = 97), color = palette[4]) + 
  geom_ribbon(aes(ymin = X0.025quant, ymax = X0.975quant), color = palette[1], fill = palette[1], alpha = 0.3) + 
  geom_point(aes(y = mean), color = palette[1]) + 
  theme_classic() + 
  labs(x = "c", y = "", title = "Gamma")

p.random <- (p.alpha | p.beta)/(p.kappa | p.gamma) + plot_layout(guides = "collect")
ggsave("random.pdf", p.random, path = output.path, dpi = "retina", height = 5, width = 8)

#   ----   Heatmap of epsilon   ----

p.epsilon <- ggplot(data.frame(res.inlabru$summary.random$epsilon) %>% left_join(observed %>% select(c("x","t","xt")), by = c("ID" = "xt"))) + 
  geom_tile(aes(x = x, y = t, fill = mean)) + 
  scale_fill_gradient(low = palette[2], high = palette[6]) + 
  theme_classic() + 
  labs(x = "x", y = "t")

ggsave("epsilon.pdf", p.epsilon, path = output.path, dpi = "retina", height = 3, width = 4)

#   ----   Plot precisions of random effects   ----

p.alpha.prec <- ggplot(data.frame(res.inlabru$marginals.hyperpar$`Precision for alpha`) %>% filter(x < 4)) + 
  geom_area(aes(x=x, y=y), color = palette[1], fill = palette[1], alpha = 0.5) + 
  theme_classic() + 
  labs(title = "Alpha", x="", y = "")

p.beta.prec <- ggplot(data.frame(res.inlabru$marginals.hyperpar$`Precision for beta`) %>% filter(x < 750)) + 
  geom_area(aes(x=x, y=y), color = palette[1], fill = palette[1], alpha = 0.5) + 
  theme_classic() + 
  labs(title = "Beta", x="", y = "")

p.kappa.prec <- ggplot(data.frame(res.inlabru$marginals.hyperpar$`Precision for kappa`) %>% filter(x < 600)) + 
  geom_area(aes(x=x, y=y), color = palette[1], fill = palette[1], alpha = 0.5) + 
  theme_classic() + 
  labs(title = "Kappa", x="", y = "")

p.gamma.prec <- ggplot(data.frame(res.inlabru$marginals.hyperpar$`Precision for gamma`) %>% filter(x < 100)) + 
  geom_area(aes(x=x, y=y), color = palette[1], fill = palette[1], alpha = 0.5) + 
  theme_classic() + 
  labs(title = "Gamma", x="", y = "")

p.epsilon.prec <- ggplot(data.frame(res.inlabru$marginals.hyperpar$`Precision for epsilon`) %>% filter(x < 150000)) + 
  geom_area(aes(x=x, y=y), color = palette[1], fill = palette[1], alpha = 0.5) + 
  theme_classic() + 
  labs(title = "Epsilon", x="", y = "")

p.prec <- (p.alpha.prec | p.beta.prec)/(p.kappa.prec | p.gamma.prec |p.epsilon.prec) + plot_layout(guides = "collect")
ggsave("precisions.pdf", p.prec, path = output.path, dpi = "retina", height = 5, width = 8)
