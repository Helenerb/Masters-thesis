# functions for plotting analysis results for real data
library(ggplot2)
library(patchwork)

source("Scripts/Misc/palette.R")
source("Scripts/Functions/plotters.R")

plot.inlabru.real <- function(res.inlabru, cancer.data, save=FALSE, path.to.storage=""){
  #' Plots inlabru results from real data
  #' 
  #' @param res.inlabru (bru object) result from inlabru analysis
  #' @param cancer.data (data.frame) data frame containing observations used in analyisis
  #' @param save (boolean) whether or not to save plots locally
  #' @param path.to.storage ("/path/to/storage") path to local folder where figures should be stored if save=TRUE
  
  #   ----   plot random effects   ----
  
  p.alpha <- ggplot(data = res.inlabru$summary.random$alpha, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.4, fill = palette[1]) + 
    geom_point(aes(y = mean), size = 0.5, color = palette[1], fill = palette[1]) + 
    labs(title="Alpha", x = "x", y='')
  
  p.beta <- ggplot(data = res.inlabru$summary.random$beta, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.4, fill = palette[1]) + 
    geom_point(aes(y = mean), size = 0.5, color = palette[1], fill = palette[1]) + 
    labs(x = "x", y = "beta", title = "Beta")
  
  p.kappa <- ggplot(data = res.inlabru$summary.random$kappa, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.4, fill = palette[1]) + 
    geom_point(aes(y = mean), size = 0.5, color = palette[1], fill = palette[1]) + 
    labs(x = "t", y = "kappa", title = "Kappa")
  
  p.intercept <- ggplot(data.frame(res.inlabru$marginals.fixed)) + 
    geom_area(aes(x = Int.x, y = Int.y), alpha = 0.4, fill = palette[1]) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1]), color = palette[1], fill=palette[1]) + 
    labs(x = "Value of intercept", y = " ", title = "Intercept")
  
  p.gamma <- ggplot(data = res.inlabru$summary.random$gamma, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.4, fill = palette[1]) + 
    geom_point(aes(y = mean), size = 0.5, color = palette[1], fill = palette[1]) + 
    labs(x = "t", y = "gamma", title = "Gamma")
  
  p.random.effects <- (p.intercept | p.alpha | p.beta) / (p.kappa | p.gamma) + 
    plot_layout(guides = "collect")
  
  if(save){
    save.figure(p.random.effects, name = "random_effects_inlabru", path = path.to.storage)
  }
  
  #   ----   plot estimated mortality rate with the observed mortality rate   ----
  
  data.mr <- data.frame(mr.sim = res.inlabru$summary.fitted.values$mean[1:length(cancer.data$`mortality rate`)]) %>%
    mutate(sim.975 = res.inlabru$summary.fitted.values$`0.975quant`[1:length(cancer.data$`mortality rate`)]) %>%
    mutate(sim.025 = res.inlabru$summary.fitted.values$`0.025quant`[1:length(cancer.data$`mortality rate`)]) %>%
    mutate(true.mr = cancer.data$`mortality rate`) %>%
    mutate(xt = cancer.data$xt, x = cancer.data$x, t = cancer.data$t)
  
  p.mr <- ggplot(data = data.mr) +
    geom_point(aes(x = mr.sim, y = true.mr), color = palette[1]) + 
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    labs(x="Estimated mortality rate", y="Observed mortality rate", title = "Mortality rate")
  
  p.mr.2 <- ggplot(data = data.mr) +
    geom_line(aes(x=xt, y = mr.sim, color="Estimated", fill = "Estimated"), size = 0.5) +
    geom_ribbon(aes(x = xt, ymin = sim.025, ymax = sim.975, fill = "Estimated"), alpha = 0.4) + 
    geom_line(aes(x=xt, y = true.mr, color="Observed", fill = "Observed"), size = 0.5) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    labs(x=" ", y="Mortality rate", title="Mortality rate")
  
  p.mr.t <- ggplot(data = data.mr) + 
    geom_point(aes(x = x, y = mr.sim, color = "Estimated", fill = "Estimated"), size = 0.5) +
    geom_ribbon(aes(x = x, ymin = sim.025, ymax = sim.975, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(x = x, y = true.mr, color = "Observed", fill = "Observed"), size = 0.5) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    labs(x = " ", y = " ", title = "Mortality rate, for each year") + 
    facet_wrap(~t)
  
  p.mr.x <- ggplot(data = data.mr) + 
    geom_point(aes(x = t, y = mr.sim, color = "Estimated", fill = "Estimated"), size = 0.5) +
    geom_ribbon(aes(x = t, ymin = sim.025, ymax = sim.975, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(x = t, y = true.mr, color = "Observed", fill = "Observed"), size = 0.5) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    labs(x = " ", y = " ", title = "Mortality rate, for each age") + 
    facet_wrap(~x)
  
  p.mr.xt <- (p.mr | p.mr.2) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.mr.xt, name="mr_xt_inlabru", path=path.to.storage)
  }
  
  p.mr.facet <- (p.mr.t | p.mr.x) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.mr.facet, name="mr_facet_inlabru", path=path.to.storage)
  }
  
  #   ----   Plot hyperparameters   ----   
  
  plots <- list(p.alpha = p.alpha, p.beta = p.beta, p.kappa = p.kappa,
                p.intercept = p.intercept, p.gamma = p.gamma,
                p.random.effects = p.random.effects)
  return(plots)
  
}

plot.hypers.inlabru.real <- function(res.inlabru, cancer.data, cutoff_alpha = 10,
                                     cutoff_beta = 1000, cutoff_kappa = 1500, cutoff_gamma = 3000,
                                     cutoff_epsilon = 50000, save=FALSE,
                                     path.to.storage=""){
  data.hyperpar <- data.frame(res.inlabru$marginals.hyperpar)
  
  p.alpha.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.alpha.x < cutoff_alpha)) + 
    geom_area(aes(x = Precision.for.alpha.x, y = Precision.for.alpha.y), color = palette[1], fill = palette[1], alpha = 0.5) + 
    geom_vline(aes(xintercept = res.inlabru$summary.hyperpar$mean[1]), color = palette[1]) + 
    labs(x = " ", y = " ", title = "Precision of alpha")
  
  p.beta.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.beta.x < cutoff_beta)) + 
    geom_area(aes(x = Precision.for.beta.x, y = Precision.for.beta.y), color = palette[1],fill = palette[1], alpha = 0.5) + 
    geom_vline(aes(xintercept = res.inlabru$summary.hyperpar$mean[2]), color = palette[1]) + 
    labs(x = " ", y = " ", title = "Precision of beta")
  
  p.kappa.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.kappa.x < cutoff_kappa)) + 
    geom_area(aes(x = Precision.for.kappa.x, y = Precision.for.kappa.y), color = palette[1], fill = palette[1], alpha = 0.5) + 
    geom_vline(aes(xintercept = res.inlabru$summary.hyperpar$mean[3]), color = palette[1]) + 
    labs(x = " ", y = " ", title = "Precision of kappa")
  
  p.gamma.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.gamma.x < cutoff_gamma)) + 
    geom_area(aes(x = Precision.for.gamma.x, y = Precision.for.gamma.y), color = palette[1], fill = palette[1], alpha = 0.5) + 
    geom_vline(aes(xintercept = res.inlabru$summary.hyperpar$mean[4]), color = palette[1]) + 
    labs(x = " ", y = " ", title = "Precision of gamma")
  
  p.epsilon.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.epsilon.x < cutoff_epsilon)) + 
    geom_area(aes(x = Precision.for.epsilon.x, y = Precision.for.epsilon.y), color = palette[1], fill = palette[1], alpha = 0.5) + 
    geom_vline(aes(xintercept = res.inlabru$summary.hyperpar$mean[5]), color = palette[1]) + 
    labs(x = " ", y = " ", title = "Precision of epsilon")
  
  p.hyperpars <- (p.alpha.prec | p.beta.prec | p.kappa.prec)/(p.gamma.prec | p.epsilon.prec)
  
  if(save){
    save.figure(p.hyperpars, name="hypers_inlabru", path=path.to.storage)
  }
}