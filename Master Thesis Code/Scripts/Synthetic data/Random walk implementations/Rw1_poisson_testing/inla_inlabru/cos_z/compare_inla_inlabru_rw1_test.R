# script for comparison of inla and inlabru  - random walk models

set_workspace <- function(markov=TRUE){
  if(markov){
    .libPaths("~/Documents/R_libraries")
    setwd("~/Documents/GitHub/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Random walk implementations/Rw1_poisson_testing/inla_inlabru/cos_z")
  } else {
    setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data/Random walk implementations/Rw1_poisson_testing/inla_inlabru/cos_z")
  }
}

set_workspace(markov=F)

library("rstan")
library("inlabru")
library("tidyverse")
library("ggplot2")
library("patchwork")
library("INLA")

# generate test data:
n=100
z=seq(0,6,length.out=n)
#y=sin(z)+rnorm(n,mean=0,sd=0.5)
offset = rep(50, length(z))
# y=rpois(length(z), offset*exp(sin(z)))
# data=data.frame(y=y,z=z, x=z, eta = sin(z))
y=rpois(length(z), offset*exp(cos(z)))
data=data.frame(y=y,z=z, x=z, eta = cos(z))

run.inlabru.fh <- function(data, E, initial = log(75)){
  # quickly estimate in inlabru:
  # components = ~ -1 +
  #   eta(z, model = "rw1", constr = TRUE, hyper = list(prec = list(intitial = log(75), fixed = TRUE)), scale.model = T)
  
  # scale with 1/pi...
  hyper.eta <- list(prec = list(initial = initial, fixed = T))
  
  # components = ~ -1 +
  #   eta(x, model = "rw1", constr = TRUE, hyper = hyper.eta)
  
  components = ~ eta(x, model = "rw1", constr = TRUE, hyper = hyper.eta, scale.model = T)
  
  formula = y ~ eta
  likelihood = like(formula = formula, family = "poisson", data = data, E = E)
  
  c.compute <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE, return.marginals.predictor = TRUE)  # control.compute
  
  results <- bru(components = components,
                 likelihood = likelihood, 
                 options = list(verbose = F, 
                                bru_verbose = 1,
                                control.compute = c.compute))
  return(results)
}

res.inlabru.fh <- run.inlabru.fh(data = data, E = offset, initial = 8)

run.inla.fh <- function(data, E, initial = log(75)){
  hyper.eta.inla <- list(prec = list(initial = initial, fixed = T))
  formula  = y ~f(x, model = "rw1", constr = TRUE, hyper = hyper.eta.inla, scale.model = T)
  result <- inla(formula, family  ="poisson", data = data, E = E)
  
  return(result)
}

res.inla.fh <- run.inla.fh(data=data, E = offset, initial = 8)

palette <- c('#70A4D4', '#ECC64B', '#93AD80', '#da9124', '#696B8D',
             '#3290c1', '#5d8060', '#D7B36A', '#826133', '#A85150')

p.inla_inlabru.fh <- ggplot() + 
  geom_point(data = data.frame(res.inlabru.fh$summary.random$eta), aes(x = ID, y = mean, color = "Inlabru")) + 
  geom_point(data = data.frame(res.inla.fh$summary.random$x), aes(x = ID, y = mean, color = "Inla")) + 
  geom_point(data = data, aes(x = z, y = eta, color = "True"), shape = 1) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  labs(title = "Comparison between Inlabru and INLA, fixed hypers, initial = 8", x = "", y = "")

ggsave("inla_inlabru_fh.pdf", plot = p.inla_inlabru.fh, dpi="retina", device = "pdf", height = 5, width = 8)

###   ----   Plot marginals and intercept:   ----

p.int.fh <- ggplot() + 
  geom_area(data = data.frame(res.inlabru.fh$marginals.fixed$Intercept), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_area(data = data.frame(res.inla.fh$marginals.fixed$`(Intercept)`), aes(x = x, y = y, color = "Inla", fill = "Inla"), alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name="", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title = "Intercept", x = "", y = "")

p.eta.fh.0 <- ggplot() + 
  geom_area(data = data.frame(res.inlabru.fh$marginals.random$eta$index.1), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_area(data = data.frame(res.inla.fh$marginals.random$x$index.1), aes(x = x, y = y, color = "Inla", fill = "Inla"), alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name="", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title ="z = 0", x = "", y = "")

p.eta.fh.2 <- ggplot() + 
  geom_area(data = data.frame(res.inlabru.fh$marginals.random$eta$index.34), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_area(data = data.frame(res.inla.fh$marginals.random$x$index.34), aes(x = x, y = y, color = "Inla", fill = "Inla"), alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name="", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title="z = 2", x = "", y = "")

p.eta.fh.4 <- ggplot() + 
  geom_area(data = data.frame(res.inlabru.fh$marginals.random$eta$index.67), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_area(data = data.frame(res.inla.fh$marginals.random$x$index.67), aes(x = x, y = y, color = "Inla", fill = "Inla"), alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name="", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title="z = 4", x = "", y = "")

p.eta.fh.6 <- ggplot() + 
  geom_area(data = data.frame(res.inlabru.fh$marginals.random$eta$index.100), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_area(data = data.frame(res.inla.fh$marginals.random$x$index.100), aes(x = x, y = y, color = "Inla", fill = "Inla"), alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name="", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title = "z = 6", x = "", y = "")

p.margs.fh <- (p.int.fh | p.eta.fh.0 | p.eta.fh.2)/(p.eta.fh.4 | p.eta.fh.6) + 
  plot_layout(guides = "collect")

ggsave("margs_fh.pdf", plot = p.margs.fh, dpi="retina", device = "pdf", height = 5, width = 8)

####   ---- both with pc priors:   ----
run.inlabru.pc <- function(data, E){
  
  hyper.eta <- list(prec = list(prior = "pc.prec", param = c(5, 0.01)))
  
  # components = ~ -1 +
  #   eta(x, model = "rw1", constr = TRUE, hyper = hyper.eta)
  
  components = ~ eta(x, model = "rw1", constr = TRUE, hyper = hyper.eta, scale.model = T)
  
  formula = y ~ eta
  likelihood = like(formula = formula, family = "poisson", data = data, E = E)
  
  c.compute <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE, return.marginals.predictor = TRUE)  # control.compute
  
  results <- bru(components = components,
                 likelihood = likelihood, 
                 options = list(verbose = F, 
                                bru_verbose = 1,
                                control.compute = c.compute))
  return(results)
}

res.inlabru.pc <- run.inlabru.pc(data = data, E = offset)

run.inla.pc <- function(data, E, initial = log(75)){
  hyper.eta.inla <- list(prec = list(prior = "pc.prec", param = c(5, 0.01)))
  formula  = y ~ f(x, model = "rw1", constr = TRUE, hyper = hyper.eta.inla, scale.model = T)
  result <- inla(formula, family  ="poisson", data = data, E = E)
  
  return(result)
}

res.inla.pc <- run.inla.fh(data=data, E = offset)

p.inla_inlabru.pc <- ggplot() + 
  geom_point(data = data.frame(res.inlabru.pc$summary.random$eta), aes(x = ID, y = mean, color = "Inlabru")) + 
  geom_point(data = data.frame(res.inla.pc$summary.random$x), aes(x = ID, y = mean, color = "Inla")) + 
  geom_point(data = data, aes(x = z, y = eta, color = "True"), shape = 1) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  labs(title = "Comparison between Inlabru and INLA, pc priors", x = "", y = "")

ggsave("inla_inlabru_pc.pdf", plot = p.inla_inlabru.pc, dpi="retina", device = "pdf", height = 5, width = 8)

###   ----   Plot marginals and intercept:   ----

p.int.pc <- ggplot() + 
  geom_area(data = data.frame(res.inlabru.pc$marginals.fixed$Intercept), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_area(data = data.frame(res.inla.pc$marginals.fixed$`(Intercept)`), aes(x = x, y = y, color = "Inla", fill = "Inla"), alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name="", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title = "Intercept", x = "", y = "")

p.eta.pc.0 <- ggplot() + 
  geom_area(data = data.frame(res.inlabru.pc$marginals.random$eta$index.1), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_area(data = data.frame(res.inla.pc$marginals.random$x$index.1), aes(x = x, y = y, color = "Inla", fill = "Inla"), alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name="", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title ="z = 0", x = "", y = "")

p.eta.pc.2 <- ggplot() + 
  geom_area(data = data.frame(res.inlabru.pc$marginals.random$eta$index.34), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_area(data = data.frame(res.inla.pc$marginals.random$x$index.34), aes(x = x, y = y, color = "Inla", fill = "Inla"), alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name="", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title="z = 2", x = "", y = "")

p.eta.pc.4 <- ggplot() + 
  geom_area(data = data.frame(res.inlabru.pc$marginals.random$eta$index.67), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_area(data = data.frame(res.inla.pc$marginals.random$x$index.67), aes(x = x, y = y, color = "Inla", fill = "Inla"), alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name="", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title="z = 4", x = "", y = "")

p.eta.pc.6 <- ggplot() + 
  geom_area(data = data.frame(res.inlabru.pc$marginals.random$eta$index.100), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_area(data = data.frame(res.inla.pc$marginals.random$x$index.100), aes(x = x, y = y, color = "Inla", fill = "Inla"), alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name="", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title = "z = 6", x = "", y = "")

p.margs.pc <- (p.int.pc | p.eta.pc.0 | p.eta.pc.2)/(p.eta.pc.4 | p.eta.pc.6) + 
  plot_layout(guides = "collect")

ggsave("margs_pc.pdf", plot = p.margs.pc, dpi="retina", device = "pdf", height = 5, width = 8)



###   ----  both with default priors:   ----
run.inlabru.default <- function(data, E){
  
  #hyper.eta <- list(prec = list(prior = "pc.prec", param = c(5, 0.01)))
  
  components = ~ eta(x, model = "rw1", constr = TRUE, scale.model = T)
  
  formula = y ~ eta
  likelihood = like(formula = formula, family = "poisson", data = data, E = E)
  
  c.compute <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE, return.marginals.predictor = TRUE)  # control.compute
  
  results <- bru(components = components,
                 likelihood = likelihood, 
                 options = list(verbose = F, 
                                bru_verbose = 1,
                                control.compute = c.compute))
  return(results)
}

res.inlabru.default <- run.inlabru.default(data = data, E = offset)

run.inla.default <- function(data, E, initial = log(75)){
  #hyper.eta.inla <- list(prec = list(prior = "pc.prec", param = c(5, 0.01)))
  formula  = y ~f(x, model = "rw1", constr = TRUE, scale.model = T)
  result <- inla(formula, family  ="poisson", data = data, E = E)
  
  return(result)
}

res.inla.default <- run.inla.default(data=data, E = offset)

p.inla_inlabru.default <- ggplot() + 
  geom_point(data = data.frame(res.inlabru.default$summary.random$eta), aes(x = ID, y = mean, color = "Inlabru")) + 
  geom_point(data = data.frame(res.inla.default$summary.random$x), aes(x = ID, y = mean, color = "Inla")) + 
  geom_point(data = data, aes(x = z, y = eta, color = "True"), shape = 1) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  labs(title = "Comparison between Inlabru and INLA, default priors", x = "", y = "")

ggsave("inla_inlabru_default.pdf", plot = p.inla_inlabru.default, dpi="retina", device = "pdf", height = 5, width = 8)

###   ----   Plot marginals and intercept:   ----

res.inlabru.df <- res.inlabru.default
res.inla.df <- res.inla.default

p.int.df <- ggplot() + 
  geom_area(data = data.frame(res.inlabru.df$marginals.fixed$Intercept), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_area(data = data.frame(res.inla.df$marginals.fixed$`(Intercept)`), aes(x = x, y = y, color = "Inla", fill = "Inla"), alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name="", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title = "Intercept", x = "", y = "")

p.prec.df <- ggplot() + 
  geom_area(data = data.frame(res.inlabru.df$marginals.hyperpar$`Precision for eta`) %>% filter(x < 30), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_area(data = data.frame(res.inla.df$marginals.hyperpar$`Precision for x`) %>% filter(x < 30), aes(x = x, y = y, color = "Inla", fill = "Inla"), alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name="", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title = "Precision", x = "", y = "")

p.eta.df.0 <- ggplot() + 
  geom_area(data = data.frame(res.inlabru.df$marginals.random$eta$index.1), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_area(data = data.frame(res.inla.df$marginals.random$x$index.1), aes(x = x, y = y, color = "Inla", fill = "Inla"), alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name="", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title ="z = 0", x = "", y = "")

p.eta.df.2 <- ggplot() + 
  geom_area(data = data.frame(res.inlabru.df$marginals.random$eta$index.34), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_area(data = data.frame(res.inla.df$marginals.random$x$index.34), aes(x = x, y = y, color = "Inla", fill = "Inla"), alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name="", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title="z = 2", x = "", y = "")

p.eta.df.4 <- ggplot() + 
  geom_area(data = data.frame(res.inlabru.df$marginals.random$eta$index.67), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_area(data = data.frame(res.inla.df$marginals.random$x$index.67), aes(x = x, y = y, color = "Inla", fill = "Inla"), alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name="", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title="z = 4", x = "", y = "")

p.eta.df.6 <- ggplot() + 
  geom_area(data = data.frame(res.inlabru.df$marginals.random$eta$index.100), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_area(data = data.frame(res.inla.df$marginals.random$x$index.100), aes(x = x, y = y, color = "Inla", fill = "Inla"), alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name="", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title = "z = 6", x = "", y = "")

p.margs.df <- (p.int.df | p.eta.df.0 | p.eta.df.2)/(p.eta.df.4 | p.eta.df.6 | p.prec.df) + 
  plot_layout(guides = "collect")

ggsave("margs_df.pdf", plot = p.margs.df, dpi="retina", device = "pdf", height = 5, width = 8)






###   ----  both with "default" log-gamma priors:   ----
run.inlabru.lg <- function(data, E){
  
  hyper.eta <- list(prec = list(prior = "loggamma", param = c(1, 0.00005)))
  
  components = ~ eta(x, model = "rw1", constr = TRUE, hyper = hyper.eta, scale.model = T)
  
  formula = y ~ eta
  likelihood = like(formula = formula, family = "poisson", data = data, E = E)
  
  c.compute <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE, return.marginals.predictor = TRUE)  # control.compute
  
  results <- bru(components = components,
                 likelihood = likelihood, 
                 options = list(verbose = F, 
                                bru_verbose = 1,
                                control.compute = c.compute))
  return(results)
}

res.inlabru.lg <- run.inlabru.lg(data = data, E = offset)

run.inla.lg <- function(data, E){
  hyper.eta.inla <- list(prec = list(prior = "loggamma", param = c(1, 0.00005)))
  formula  = y ~f(x, model = "rw1", constr = TRUE, hyper = hyper.eta.inla, scale.model = T)
  result <- inla(formula, family  ="poisson", data = data, E = E)
  
  return(result)
}

res.inla.lg <- run.inla.lg(data=data, E = offset)

p.inla_inlabru.lg <- ggplot() + 
  geom_point(data = data.frame(res.inlabru.lg$summary.random$eta), aes(x = ID, y = mean, color = "Inlabru")) + 
  geom_point(data = data.frame(res.inla.lg$summary.random$x), aes(x = ID, y = mean, color = "Inla")) + 
  geom_point(data = data, aes(x = z, y = eta, color = "True"), shape = 1) + 
  theme_classic() + 
  scale_color_manual(name = "", values = palette) + 
  labs(title = "Comparison between Inlabru and INLA, default loggamma priors", x = "", y = "")

ggsave("inla_inlabru_lg.pdf", plot = p.inla_inlabru.lg, dpi="retina", device = "pdf", height = 5, width = 8)

###   ----   Plot marginals and intercept:   ----

p.int.lg <- ggplot() + 
  geom_area(data = data.frame(res.inlabru.lg$marginals.fixed$Intercept), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_area(data = data.frame(res.inla.lg$marginals.fixed$`(Intercept)`), aes(x = x, y = y, color = "Inla", fill = "Inla"), alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name="", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title = "Intercept", x = "", y = "")

p.prec.lg <- ggplot() + 
  geom_area(data = data.frame(res.inlabru.lg$marginals.hyperpar$`Precision for eta`) %>% filter(x < 30), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_area(data = data.frame(res.inla.lg$marginals.hyperpar$`Precision for x`) %>% filter(x < 30), aes(x = x, y = y, color = "Inla", fill = "Inla"), alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name="", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title = "Precision", x = "", y = "")

p.eta.lg.0 <- ggplot() + 
  geom_area(data = data.frame(res.inlabru.lg$marginals.random$eta$index.1), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_area(data = data.frame(res.inla.lg$marginals.random$x$index.1), aes(x = x, y = y, color = "Inla", fill = "Inla"), alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name="", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title ="z = 0", x = "", y = "")

p.eta.lg.2 <- ggplot() + 
  geom_area(data = data.frame(res.inlabru.lg$marginals.random$eta$index.34), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_area(data = data.frame(res.inla.lg$marginals.random$x$index.34), aes(x = x, y = y, color = "Inla", fill = "Inla"), alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name="", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title="z = 2", x = "", y = "")

p.eta.lg.4 <- ggplot() + 
  geom_area(data = data.frame(res.inlabru.lg$marginals.random$eta$index.67), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_area(data = data.frame(res.inla.lg$marginals.random$x$index.67), aes(x = x, y = y, color = "Inla", fill = "Inla"), alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name="", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title="z = 4", x = "", y = "")

p.eta.lg.6 <- ggplot() + 
  geom_area(data = data.frame(res.inlabru.lg$marginals.random$eta$index.100), aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_area(data = data.frame(res.inla.lg$marginals.random$x$index.100), aes(x = x, y = y, color = "Inla", fill = "Inla"), alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(name="", values = palette) + 
  scale_color_manual(name = "", values = palette) + 
  labs(title = "z = 6", x = "", y = "")

p.margs.lg <- (p.int.lg | p.eta.lg.0 | p.eta.lg.2)/(p.eta.lg.4 | p.eta.lg.6 | p.prec.lg) + 
  plot_layout(guides = "collect")

ggsave("margs_lg.pdf", plot = p.margs.lg, dpi="retina", device = "pdf", height = 5, width = 8)




