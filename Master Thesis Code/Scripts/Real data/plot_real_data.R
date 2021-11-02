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

plot.inlabru.real.predicted <- function(res.inlabru, cancer.data, last.obs.t, save=FALSE, path.to.storage=""){
  #' Plots inlabru results from real data
  #' 
  #' @param res.inlabru (bru object) result from inlabru analysis
  #' @param cancer.data (data.frame) data frame containing observations used in analyisis
  #' @param save (boolean) whether or not to save plots locally
  #' @param path.to.storage ("/path/to/storage") path to local folder where figures should be stored if save=TRUE
  #' @param last.obs.t (int) the last observed value, given by the t index. 
  
  
  x.age.map <- cancer.data %>% select(c())
  
  #   ----   plot random effects   ----
  
  data.alpha <- left_join(res.inlabru$summary.random$alpha, {cancer.data %>% select(age, x, age.int)}, by=c("ID" = "x"))
  p.alpha <- ggplot(data = data.alpha, aes(x = age.int)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = palette[1], alpha = 0.4) + 
    geom_point(aes(y = mean), size = 0.5, color = palette[1], fill = palette[1]) + 
    labs(title="Alpha", x = "x", y='')
  
  if(save){
    save.figure(p.alpha, name = "alpha_inlabru", path = path.to.storage)
  }
  
  data.beta <- left_join(res.inlabru$summary.random$beta, {cancer.data %>% select(age, x, age.int)}, by=c("ID" = "x"))
  p.beta <- ggplot(data = data.beta, aes(x = age.int)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.4, fill = palette[1]) + 
    geom_point(aes(y = mean), size = 0.5, color = palette[1], fill = palette[1]) + 
    labs(x = "x", y = "beta", title = "Beta")
  if(save){
    save.figure(p.beta, name = "beta_inlabru", path = path.to.storage)
  }
  
  data.kappa <- left_join(res.inlabru$summary.random$kappa, {cancer.data %>% select(year, t, predict)}, by=c("ID" = "t")) %>%
    mutate(year = parse_integer(year))
  print(data.kappa)
  p.kappa <- ggplot(data = data.kappa, aes(x = year)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.4, fill = palette[1]) + 
    geom_point(aes(y = mean, color = predict, fill = predict, shape=predict), size = 0.5) + 
    #geom_vline(aes(xintercept = last.obs.t), color = palette[2]) + 
    scale_color_manual(values = palette ) +
    scale_fill_manual(values = palette ) +
    scale_shape_manual(values = c(3,2)) + 
    labs(x = "t", y = "kappa", title = "Kappa")
  if(save){
    save.figure(p.kappa, name = "kappa_inlabru", path = path.to.storage)
  }
  
  p.intercept <- ggplot(data.frame(res.inlabru$marginals.fixed)) + 
    geom_area(aes(x = Int.x, y = Int.y), alpha = 0.4, fill = palette[1]) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1]), color = palette[1], fill=palette[1]) + 
    labs(x = "Value of intercept", y = " ", title = "Intercept")
  if(save){
    save.figure(p.intercept, name = "intercept_inlabru", path = path.to.storage)
  }
  
  data.gamma <- res.inlabru$summary.random$gamma
  c.predicted <- {cancer.data %>% filter(predict == "predicted")}$c
  c.observed <- {cancer.data %>% filter(predict == "observed")}$c
  
  fully.observed <- data.gamma %>%
    filter(!ID %in% c.predicted) %>% mutate(predicted = "Fully observed")
  
  partially.observed <- data.gamma %>%
    filter(ID %in% c.predicted & ID %in% c.observed) %>%
    mutate(predicted = "Partially predicted")
  
  fully.predicted <- data.gamma %>%
    filter(!ID %in% c.observed) %>% mutate(predicted = "Fully predicted")
  
  data.gamma <- rbind(fully.observed, partially.observed, fully.predicted)
  
  p.gamma <- ggplot(data = data.gamma, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.4, fill = palette[1]) + 
    geom_point(aes(y = mean, color = predicted, fill = predicted, shape=predicted), size = 0.5) + 
    scale_color_manual(values = palette ) +
    scale_fill_manual(values = palette ) +
    scale_shape_manual(values = c(3,2,4)) + 
    labs(x = "cohort", y = "gamma", title = "Gamma")
  
  if(save){
    save.figure(p.gamma, name = "gamma_inlabru", path = path.to.storage)
  }
  
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
    mutate(xt = cancer.data$xt, x = cancer.data$x, t = cancer.data$t) %>%
    mutate(year = cancer.data$year, age=cancer.data$age, age.int = cancer.data$age.int) %>%
    mutate(year = parse_integer(year)) %>%
    mutate(predict = cancer.data$predict)
  
  p.mr <- ggplot(data = data.mr) +
    geom_point(aes(x = mr.sim, y = true.mr, color=predict, shape=predict)) + 
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    scale_shape_manual(name="", values=c(3,2)) + 
    labs(x="Estimated mortality rate", y="Observed mortality rate", title = "Mortality rate")
  
  if(save){
    save.figure(p.mr, name = "mr_inlabru", path = path.to.storage)
  }
  
  p.mr.2 <- ggplot(data = data.mr) +
    geom_line(aes(x = xt, y = true.mr, color="Observed rate", fill = "Observed rate"), size = 0.5) +
    geom_point(aes(x = xt, y = mr.sim, color="Estimated", fill = "Estimated", shape = predict), size = 0.5) +
    geom_ribbon(aes(x = xt, ymin = sim.025, ymax = sim.975, fill = "Estimated"), alpha = 0.4) + 
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    scale_shape_manual(name = "", values=c(3,2)) + 
    labs(x = " ", y = "Mortality rate", title = "Mortality rate")
  
  if(save){
    save.figure(p.mr.2, name = "mr_2_inlabru", path = path.to.storage)
  }
  
  p.mr.t <- ggplot(data = data.mr) + 
    geom_ribbon(aes(x = age.int, ymin = sim.025, ymax = sim.975, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(x = age.int, y = mr.sim, color = "Estimated", fill = "Estimated", shape=predict), size = 0.5) +
    geom_line(aes(x = age.int, y = true.mr, color = "Observed rate", fill = "Observed rate"), size = 0.5, alpha = 0.9) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    scale_shape_manual(name = "", values=c(3,2,4)) + 
    labs(x = " ", y = " ", title = "Mortality rate, for each year") + 
    facet_wrap(~year)
  
  if(save){
    save.figure(p.mr.t, name = "mr_t_inlabru", path = path.to.storage)
  }
  
  p.mr.x <- ggplot(data = data.mr) + 
    geom_ribbon(aes(x = year, ymin = sim.025, ymax = sim.975, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(x = year, y = mr.sim, color = "Estimated", fill = "Estimated", shape=predict), size = 0.5) +
    geom_line(aes(x = year, y = true.mr, color = "Observed rate", fill = "Observed rate"), size = 0.5, alpha = 0.9) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    scale_shape_manual(name="", values=c(3,2,4)) + 
    labs(x = " ", y = " ", title = "Mortality rate, for each age") + 
    facet_wrap(~age)
  
  if(save){
    save.figure(p.mr.x, name = "mr_x_inlabru", path = path.to.storage)
  }
  
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

#' plot.comparison.real(res.inlabru, res.stan, cancer.data, cohort = TRUE){
#'   #'Plots comparison of inlarbu and stan results for real data
#'   #'
#'   #'@param res.inlabru (bru.object): results of run with inlabru
#'   #'@param res.stan (data.frame): results of run with stan, saved as a summary dataframe
#'   #'@param cander.cata (data.frame): observed cancer (and population) data
#'   #'@param cohort (bool): Whether the analysis included a cohort effect
#'   
#'   
#' }