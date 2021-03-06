# Running full inlabru analysis on male stomach cancer data

#   ----   Load libraries and set workspace   ----

setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code")

output.path = file.path("Scripts/Real\ Data/Analyses", "stomach_male_predict")

library("tidyverse")
library("inlabru")
library("ggplot2")
library("INLA")
library("patchwork")
library("rstan")

#    ----    Load data   ----
load("Data/population-germany.Rda")
load("Data/stomachCancer-germany.Rda")

#   ----   Format data   ----
male.stomach.cancer <- cancer.data %>% select(c(age, year, male, t, age.int, x, x.c, xt, cohort, c, birth.year, male.t)) %>%
  mutate(Y = replace(male, year %in% 2011:2016, NA), E = male.t) %>%
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

res.inlabru <- run.inlabru(male.stomach.cancer)

observed <- male.stomach.cancer

source("Scripts/Misc/palette.R")

#   ----   Generate samples for Y   ----

lambda.samples <- generate(res.inlabru, male.stomach.cancer, ~ E*exp(alpha + beta*kappa + gamma + epsilon), n.samples = 10000)
Y.samples <- matrix(rpois(lambda.samples, n = 324*10000), nrow = 324, ncol = 10000)
Y.samples.df <- data.frame(Y.samples) 

Y.inlabru <- male.stomach.cancer %>%
  mutate(Y.mean = apply(Y.samples.df, 1, mean)) %>%
  mutate(Y.0.025 = apply(Y.samples.df, 1, quantile, 0.025)) %>%
  mutate(Y.0.975 = apply(Y.samples.df, 1, quantile, 0.975)) %>%
  mutate(Y.sd = apply(Y.samples.df, 1, sd)) %>%
  mutate(DSS = ((male - Y.mean)/Y.sd)^2 + 2*log(Y.sd)) %>%
  mutate(contained = if_else(male >= Y.0.025 & male < Y.0.975, 1, 0))

MDSS.all <- mean(Y.inlabru$DSS)
MDSS.x.above.5 <- mean({Y.inlabru %>% filter(x > 5)}$DSS)

MDSS.all.in.data <- mean({Y.inlabru %>% filter(year %in% 1999:2010)}$DSS)
MDSS.x.above.5.in.data <- mean({Y.inlabru %>% filter(x > 5) %>% filter(year %in% 1999:2010)}$DSS)

MDSS.all.out.data <- mean({Y.inlabru %>% filter(year %in% 2011:2016)}$DSS)
MDSS.x.above.5.out.data <- mean({Y.inlabru %>% filter(x > 5) %>% filter(year %in% 2011:2016)}$DSS)

contained.all <- mean(Y.inlabru$contained)

contained.95.in.data <- mean({Y.inlabru %>% filter(year %in% 1999:2010)}$contained)
contained.95.out.data <- mean({Y.inlabru %>% filter(year %in% 2011:2016)}$contained)

contained.95.a.5.in.data <- mean({Y.inlabru %>% filter(x > 5) %>% filter(year %in% 1999:2010)}$contained)
contained.95.a.5.out.data <- mean({Y.inlabru %>% filter(x > 5) %>% filter(year %in% 2011:2016)}$contained)

write.table(list(MDSS.all = MDSS.all,
                 MDSS.x.above.5 = MDSS.x.above.5,
                 all.in.data = MDSS.all.in.data,
                 above.5.in.data = MDSS.x.above.5.in.data,
                 all.out.data = MDSS.all.out.data,
                 above.5.out.data = MDSS.x.above.5.out.data,
                 contained.95.in.data = contained.95.in.data,
                 contained.95.a.5.in.data = contained.95.a.5.in.data,
                 contained.95.out.data = contained.95.out.data,
                 contained.95.a.5.out.data = contained.95.a.5.out.data,
                 contained.all = contained.all), file = file.path(output.path, "DSS.txt"))

p.Y.age <- ggplot(Y.inlabru %>% filter(year %in% 2011:2016)) + 
  geom_ribbon(aes(x = age.int, ymin = Y.0.025, ymax = Y.0.975, color = "Estimated", fill = "Estimated", shape = "Estimated"), alpha = 0.2, size = 0.5) + 
  geom_point(aes(x = age.int, y = Y.mean, color = "Estimated", fill = "Estimated", shape = "Estimated")) + 
  geom_point(aes(x = age.int, y = male, color = "Observed", fill = "Observed", shape = "Observed"), size = 2) + 
  facet_wrap(~ as.factor(year)) + 
  scale_color_manual(name="", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  scale_shape_manual(name = "", values = c(4,16)) + 
  theme_classic() + 
  labs(x = "Age", y = "")

ggsave("Y_by_age.pdf", p.Y.age, path = output.path, dpi = "retina", height = 5, width = 8)  

p.Y.age.in.data <- ggplot(Y.inlabru %>% filter(year %in% 1999:2010)) + 
  geom_ribbon(aes(x = age.int, ymin = Y.0.025, ymax = Y.0.975, color = "Estimated", fill = "Estimated", shape = "Estimated"), alpha = 0.2, size = 0.5) + 
  geom_point(aes(x = age.int, y = Y.mean, color = "Estimated", fill = "Estimated", shape = "Estimated")) + 
  geom_point(aes(x = age.int, y = male, color = "Observed", fill = "Observed", shape = "Observed"), size = 2) + 
  facet_wrap(~ as.factor(year)) + 
  scale_color_manual(name="", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  scale_shape_manual(name = "", values = c(4,16)) + 
  theme_classic() + 
  labs(x = "Age", y = "")

ggsave("Y_by_age_in_data.pdf", p.Y.age.in.data, path = output.path, dpi = "retina", height = 5, width = 8)  


p.Y.year <- ggplot(Y.inlabru %>% filter(age.int >= 50)) + 
  geom_ribbon(aes(x = year, ymin = Y.0.025, ymax = Y.0.975, color = "Estimated", fill = "Estimated", shape = "Estimated"), alpha = 0.2, size = 0.5) + 
  geom_point(aes(x = year, y = Y.mean, color = "Estimated", fill = "Estimated", shape = "Estimated")) + 
  geom_point(aes(x = year, y = male, color = "Observed", fill = "Observed", shape = "Observed"), size = 2) + 
  geom_vline(aes(xintercept = 2011, color="Predicted period", fill = "Predicted period", shape = "Predicted period")) + 
  facet_wrap(~ as.factor(age), ncol = 4) + 
  scale_color_manual(name="", values = palette) + 
  scale_fill_manual(name = "", values = palette) + 
  scale_shape_manual(name = "", values = c(4,16,3)) + 
  theme_classic() + 
  labs(x = "Year", y = "") + 
  theme(legend.position="bottom")

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

p.gamma.prec <- ggplot(data.frame(res.inlabru$marginals.hyperpar$`Precision for gamma`) %>% filter(x < 400)) + 
  geom_area(aes(x=x, y=y), color = palette[1], fill = palette[1], alpha = 0.5) + 
  theme_classic() + 
  labs(title = "Gamma", x="", y = "")

p.epsilon.prec <- ggplot(data.frame(res.inlabru$marginals.hyperpar$`Precision for epsilon`) %>% filter(x < 150000)) + 
  geom_area(aes(x=x, y=y), color = palette[1], fill = palette[1], alpha = 0.5) + 
  theme_classic() + 
  labs(title = "Epsilon", x="", y = "")

p.prec <- (p.alpha.prec | p.beta.prec)/(p.kappa.prec | p.gamma.prec |p.epsilon.prec) + plot_layout(guides = "collect")
ggsave("precisions.pdf", p.prec, path = output.path, dpi = "retina", height = 5, width = 8)


#   -------   Testing   --------

test.lambdas <- rpois(n = 200000, lambda = 6000); q1.l = quantile(test.lambdas, 0.025); q2.l = quantile(test.lambdas, 0.975)
test.Y <- rpois(n = 200000, lambda = test.lambdas); q1.Y = quantile(test.Y, 0.025); q2.Y = quantile(test.Y, 0.975)
ggplot(data.frame(l = test.lambdas, Y = test.Y)) + 
  geom_density(aes(x = l), color = palette[1], fill = palette[1], alpha = 0.5) + 
  geom_density(aes(x = Y), color = palette[2], fill = palette[2], alpha = 0.5) + 
  geom_vline(aes(xintercept = q1.l), color = palette[1]) + 
  geom_vline(aes(xintercept = q2.l), color = palette[1]) + 
  geom_vline(aes(xintercept = q1.Y), color = palette[2]) + 
  geom_vline(aes(xintercept = q2.Y), color = palette[2]) + 
  theme_classic()