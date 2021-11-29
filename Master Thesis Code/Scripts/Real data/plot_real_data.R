# functions for plotting analysis results for real data
library(ggplot2)
library(patchwork)
library(tidyverse)

source("Scripts/Misc/palette.R")
source("Scripts/Functions/plotters.R")

plot.inlabru.real <- function(res.inlabru, cancer.data, save=FALSE, path.to.storage="", cohort=TRUE){
  #' Plots inlabru results from real data
  #' 
  #' @param res.inlabru (bru object) result from inlabru analysis
  #' @param cancer.data (data.frame) data frame containing observations used in analyisis
  #' @param save (boolean) whether or not to save plots locally
  #' @param path.to.storage ("/path/to/storage") path to local folder where figures should be stored if save=TRUE
  #' @param cohort (boolean) whether or not a cohort effect is included in the analysis
  
  #   ----   plot random effects   ----
  
  data.alpha <- left_join(res.inlabru$summary.random$alpha, {cancer.data %>% select(age, x, age.int)}, by=c("ID" = "x"))
  p.alpha <- ggplot(data = data.alpha, aes(x = age.int)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.4, fill = palette[1]) + 
    geom_point(aes(y = mean), size = 0.5, color = palette[1], fill = palette[1]) + 
    labs(title="Alpha", x = "age", y='')
  
  data.beta <- left_join(res.inlabru$summary.random$beta, {cancer.data %>% select(age, x, age.int)}, by=c("ID" = "x"))
  p.beta <- ggplot(data = data.beta, aes(x = age.int)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.4, fill = palette[1]) + 
    geom_point(aes(y = mean), size = 0.5, color = palette[1], fill = palette[1]) + 
    labs(x = "x", y = "", title = "Beta")
  
  data.kappa <- left_join(res.inlabru$summary.random$kappa, {cancer.data %>% select(year, t)}, by=c("ID" = "t")) %>%
    mutate(year = parse_integer(year))
  p.kappa <- ggplot(data = data.kappa, aes(x = year)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.4, fill = palette[1]) + 
    geom_point(aes(y = mean), size = 0.5, color = palette[1], fill = palette[1]) + 
    labs(x = "t", y = "", title = "Kappa")
  
  p.intercept <- ggplot(data.frame(res.inlabru$marginals.fixed)) + 
    geom_area(aes(x = Int.x, y = Int.y), alpha = 0.4, fill = palette[1]) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1]), color = palette[1], fill=palette[1]) + 
    labs(x = "Value of intercept", y = " ", title = "Intercept")
  
  if(cohort){
    p.gamma <- ggplot(data = res.inlabru$summary.random$gamma, aes(x = ID)) + 
      geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.4, fill = palette[1]) + 
      geom_point(aes(y = mean), size = 0.5, color = palette[1], fill = palette[1]) + 
      labs(x = "cohort", y = "", title = "Gamma")
    
    p.random.effects <- (p.intercept | p.alpha | p.beta) / (p.kappa | p.gamma) + 
      plot_layout(guides = "collect")
  } else {
    p.random.effects <- (p.intercept | p.alpha) / (p.beta | p.kappa) + 
      plot_layout(guides = "collect")
  }
  
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
    mutate(year = parse_integer(year)) 
  
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
    geom_point(aes(x = age.int, y = mr.sim, color = "Estimated", fill = "Estimated"), size = 0.5) +
    geom_ribbon(aes(x = age.int, ymin = sim.025, ymax = sim.975, fill = "Estimated"), alpha = 0.4) + 
    geom_line(aes(x = age.int, y = true.mr, color = "Observed", fill = "Observed"), size = 0.5) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    labs(x = " ", y = " ", title = "Mortality rate, for each year") + 
    facet_wrap(~year)
  
  if(save){
    save.figure(p.mr.t, name="mr_t_inlabru", path=path.to.storage)
  }
  
  p.mr.x <- ggplot(data = data.mr) + 
    geom_point(aes(x = year, y = mr.sim, color = "Estimated", fill = "Estimated"), size = 0.5) +
    geom_ribbon(aes(x = year, ymin = sim.025, ymax = sim.975, fill = "Estimated"), alpha = 0.4) + 
    geom_line(aes(x = year, y = true.mr, color = "Observed", fill = "Observed"), size = 0.5) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    labs(x = " ", y = " ", title = "Mortality rate, for each age") + 
    facet_wrap(~age)
  
  if(save){
    save.figure(p.mr.x, name="mr_x_inlabru", path=path.to.storage)
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

plot.hypers.inlabru.real <- function(res.inlabru, cancer.data, cutoff_alpha = 10,
                                     cutoff_beta = 1000, cutoff_kappa = 1500, cutoff_gamma = 3000,
                                     cutoff_epsilon = 50000, save=FALSE,
                                     path.to.storage="", cohort=TRUE){
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
  
  if(cohort){
    p.gamma.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.gamma.x < cutoff_gamma)) + 
      geom_area(aes(x = Precision.for.gamma.x, y = Precision.for.gamma.y), color = palette[1], fill = palette[1], alpha = 0.5) + 
      geom_vline(aes(xintercept = res.inlabru$summary.hyperpar$mean[4]), color = palette[1]) + 
      labs(x = " ", y = " ", title = "Precision of gamma")
  }
  
  p.epsilon.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.epsilon.x < cutoff_epsilon)) + 
    geom_area(aes(x = Precision.for.epsilon.x, y = Precision.for.epsilon.y), color = palette[1], fill = palette[1], alpha = 0.5) + 
    geom_vline(aes(xintercept = res.inlabru$summary.hyperpar$mean[5]), color = palette[1]) + 
    labs(x = " ", y = " ", title = "Precision of epsilon")
  
  if(cohort){
    p.hyperpars <- (p.alpha.prec | p.beta.prec | p.kappa.prec)/(p.gamma.prec | p.epsilon.prec)
  } else {
    p.hyperpars <- (p.alpha.prec | p.beta.prec )/(p.kappa.prec | p.epsilon.prec)
  }
  
  
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

plot.comparison.real <- function(res.inlabru, res.stan, stan.marginals, cancer.data, path.to.storage, cohort = FALSE, save=FALSE){
  #'Plots comparison of inlarbu and stan results for real data
  #'
  #'@param res.inlabru (bru.object): results of run with inlabru
  #'@param res.stan (data.frame): results of run with stan, saved as a summary dataframe
  #'@param stan.marginals (list<array>): list of draws from stan for random effects and hyperparams
  #'@param cander.cata (data.frame): observed cancer (and population) data
  #'@param path.to.storage (string): file path where output should be stored
  #'@param cohort (bool): Whether the analysis included a cohort effect
  #'@param save (bool): whether figures should be stored
  
  intercept.marginal <- data.frame(int = stan.marginals$intercept_draws)
  intercept.inlabru <- data.frame(res.inlabru$marginals.fixed)
  
  #  ----   intercept   ----
  p.intercept <- ggplot() + 
    geom_histogram(data = intercept.marginal, aes(x = int, y = after_stat(density), color = "Stan", fill = "Stan"), bins=200, alpha = 0.5) + 
    geom_area(data=intercept.inlabru, aes(x = Int.x, y = Int.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.4, size = 0.5) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of intercept", y = " ", title = "Intercept")
  
  # ----   alpha   ----
  alpha.inlabru <- left_join(res.inlabru$summary.random$alpha, {cancer.data %>% select(age, x, age.int)}, by=c("ID" = "x"))
  
  alpha.stan <- res.stan %>%
    rownames_to_column("parameter") %>% filter(grepl('alpha', parameter)) %>%
    filter(!grepl('tau_alpha', parameter)) %>%
    filter(!grepl('alpha_raw', parameter)) %>%
    mutate(index = parse_number(parameter)) %>% mutate(index = index - 1) %>%
    left_join({cancer.data %>% select(age, x, age.int)}, by = c("index" = "x"))
  
  p.alpha <- ggplot() + 
    geom_ribbon(data=alpha.inlabru, aes(x = age.int, ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
    geom_point(data=alpha.inlabru, aes(x = age.int, y = mean, color = "Inlabru", fill = "Inlabru"), size = 0.5) + 
    geom_point(data=alpha.stan, aes(x=age.int, y=mean, fill="Stan", color="Stan"), size=0.5) + 
    geom_line(data=alpha.stan, aes(x=age.int, y=`2.5%`, fill = "Stan", color="Stan"), alpha=0.5) + 
    geom_line(data=alpha.stan, aes(x=age.int, y=`97.5%`, fill="Stan", color="Stan"), alpha=0.5) +
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Alpha", x = "age", y='')
  
  # ----   beta   ----
  beta.inlabru <- left_join(res.inlabru$summary.random$beta, {cancer.data %>% select(age, x, age.int)}, by=c("ID" = "x"))
  
  beta.stan <- res.stan %>%
    rownames_to_column("parameter") %>%
    filter(grepl('beta', parameter)) %>% 
    filter(!grepl('tau_beta', parameter)) %>%
    filter(!grepl('beta_raw', parameter)) %>%
    mutate(index = parse_number(parameter)) %>% mutate(index = index - 1) %>%
    left_join({cancer.data %>% select(age, x, age.int)}, by = c("index" = "x"))
  
  p.beta <- ggplot() + 
    geom_ribbon(data=beta.inlabru, aes(x = age.int, ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
    geom_point(data=beta.inlabru, aes(x = age.int, y = mean, color = "Inlabru", fill = "Inlabru"), size = 0.5) + 
    geom_point(data=beta.stan, aes(x=age.int, y=mean, fill="Stan", color="Stan"), size=0.5) + 
    geom_line(data=beta.stan, aes(x=age.int, y=`2.5%`, fill="Stan", color="Stan"), alpha=0.5) + 
    geom_line(data=beta.stan, aes(x=age.int, y=`97.5%`, fill="Stan", color="Stan"), alpha=0.5) +
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Beta", x = "age", y='')
  
  #   ----   kappa   ----  
  kappa.inlabru <- left_join(res.inlabru$summary.random$kappa, {cancer.data %>% select(year, t)}, by=c("ID" = "t")) %>%
    mutate(year = parse_integer(year))
  
  kappa.stan <- res.stan %>%
    rownames_to_column("parameter") %>%
    filter(grepl('kappa', parameter)) %>%
    filter(!grepl('tau_kappa', parameter)) %>%
    filter(!grepl('kappa_raw', parameter)) %>%
    mutate(index = parse_number(parameter)) %>% mutate(index = index - 1) %>%
    left_join({cancer.data %>% select(year, t)}, by=c("index" = "t")) %>%
    mutate(year = parse_integer(year))
  
  p.kappa <- ggplot() + 
    geom_ribbon(data=kappa.inlabru, aes(x = year, ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
    geom_point(data=kappa.inlabru, aes(x = year, y = mean, color = "Inlabru", fill = "Inlabru"), size = 0.5) + 
    geom_point(data=kappa.stan, aes(x=year, y=mean, fill="Stan",color="Stan"), size=0.5) + 
    geom_line(data=kappa.stan, aes(x=year, y=`2.5%`, fill="Stan", color="Stan"), alpha=0.5) + 
    geom_line(data=kappa.stan, aes(x=year, y=`97.5%`, fill="Stan", color="Stan"), alpha=0.5) +
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Kappa", x = "t", y='')
  
  
  #   ----   gamma   ----  
  if (cohort){
    # Not tested
    gamma.inlabru <- res.inlabru$summary.random$gamma
    
    gamma.stan <- res.stan %>%
      rownames_to_column("parameter") %>%
      filter(grepl('gamma', parameter)) %>%
      filter(!grepl('tau_gamma', parameter)) %>%
      filter(!grepl('gamma_raw', parameter)) %>%
      mutate(index = parse_number(parameter)) %>% mutate(index = index - 1)
    
    p.gamma <- ggplot() + 
      geom_ribbon(data=gamma.inlabru, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
      geom_point(data=gamma.inlabru, aes(x = ID, y = mean, color = "Inlabru", fill = "Inlabru"), size = 0.5) + 
      geom_point(data=gamma.stan, aes(x=index, y=mean, fill="Stan", color="Stan"), size=0.5) + 
      geom_line(data=gamma.stan, aes(x=index, y=`2.5%`, fill="Stan", color="Stan"), alpha=0.5) + 
      geom_line(data=gamma.stan, aes(x=index, y=`97.5%`, fill="Stan", color="Stan"), alpha=0.5) +
      scale_color_manual(name = "",
                         values = palette ) +
      scale_fill_manual(name = "",
                        values = palette ) +
      labs(title="Gamma", x = "cohort", y='')
  }
  
  #   ----   mortality rate   ----
  
  eta.inlabru <- data.frame(mean = res.inlabru$summary.linear.predictor$mean[1:length(cancer.data$`mortality rate`)]) %>%
    mutate(`0.975quant` = res.inlabru$summary.linear.predictor$`0.975quant`[1:length(cancer.data$`mortality rate`)]) %>%
    mutate(`0.025quant` = res.inlabru$summary.linear.predictor$`0.025quant`[1:length(cancer.data$`mortality rate`)]) %>%
    mutate(xt = cancer.data$xt, x = cancer.data$x, t = cancer.data$t) %>%
    mutate(year = cancer.data$year, age=cancer.data$age, age.int = cancer.data$age.int) %>%
    mutate(year = parse_integer(year)) 
  
  eta.stan <- res.stan %>%
    rownames_to_column("parameter") %>%
    filter(grepl('eta', parameter)) %>%
    filter(!grepl('beta', parameter)) %>%
    mutate(index = parse_number(parameter)) %>% mutate(index = index - 1) %>%
    left_join({cancer.data %>% select(xt, Y, `mortality rate`, x, age.int, age, t, year, cohort, birth.year)}, by=c("index" = "xt")) %>%
    mutate(year = parse_integer(year))
  
  eta.both <- left_join(eta.inlabru, eta.stan, by = c("xt" = "index"), suffix = c(".inlabru", ".stan"))
  
  p.eta <- ggplot() +
    geom_point(data=eta.both, aes(x = mean.inlabru, y = mean.stan), color = palette[1]) + 
    labs(x="Eta - inlabru", y="Eta - stan", title = "Eta")
  
  p.eta.2 <- ggplot() +
    geom_line(data = eta.inlabru, aes(x=xt, y = mean, color="Inlabru", fill = "Inlabru")) +
    geom_ribbon(data = eta.inlabru, aes(x = xt, ymin=`0.025quant`, ymax=`0.975quant`, fill="Inlabru"), alpha = 0.5) +
    geom_line(data = eta.stan, aes(x=index, y = mean, color="Stan", fill="Stan")) +
    geom_ribbon(data = eta.stan, aes(x=index, ymin=`2.5%`, ymax=`97.5%`, fill="Stan"), alpha=0.5) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    labs(x=" ", y="", title="Eta")
  
  p.eta.t <- ggplot() + 
    geom_line(data = eta.inlabru, aes(x = age.int, y = mean, color = "Inlabru", fill = "Inlabru")) +
    geom_ribbon(data = eta.inlabru, aes(x = age.int, ymin=`0.025quant`, ymax=`0.975quant`, fill="Inlabru"), alpha = 0.5) +
    geom_line(data = eta.stan, aes(x = age.int, y = mean, color = "Stan", fill="Stan")) +
    geom_ribbon(data = eta.stan, aes(x = age.int, ymin=`2.5%`, ymax=`97.5%`, fill="Stan"), alpha=0.5) +
    scale_color_manual(name = "", values = palette) +
    scale_fill_manual(name = "", values = palette ) +
    labs(x = " ", y = " ", title = "Eta, for each year") + 
    facet_wrap(~year)
  
  p.eta.x <- ggplot() + 
    geom_line(data = eta.inlabru, aes(x = year, y = mean, color = "Inlabru")) +
    geom_ribbon(data = eta.inlabru, aes(x = year, ymin=`0.025quant`, ymax=`0.975quant`, fill="Inlabru"), alpha = 0.5) +
    geom_line(data=eta.stan, aes(x = year, y = mean, color = "Stan")) +
    geom_ribbon(data = eta.stan, aes(x = year, ymin=`2.5%`, ymax=`97.5%`, fill="Stan"), alpha=0.5) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    labs(x = " ", y = " ", title = "Eta, for each age group") + 
    facet_wrap(~age)
  
  #   ----   hyperparameters: precisions   ----   
  
  #  tau alpha 
  
  tau.alpha.stan <- data.frame(tau = stan.marginals$tau_alpha_draws)
  tau.alpha.inlabru <- data.frame(res.inlabru$marginals.hyperpar$`Precision for alpha`)
  
  p.tau.alpha <- ggplot() + 
    geom_area(data = tau.alpha.inlabru, aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = tau.alpha.stan, aes(x = tau, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins=100) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of precision of alpha", y = " ", title = "Precision of Alpha")
  
  #  tau beta
  
  tau.beta.stan <- data.frame(tau = stan.marginals$tau_beta_draws)
  tau.beta.inlabru <- data.frame(res.inlabru$marginals.hyperpar$`Precision for beta`)
  
  p.tau.beta <- ggplot() + 
    geom_area(data = tau.beta.inlabru, aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = tau.beta.stan, aes(x = tau, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins=100) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of precision of beta", y = " ", title = "Precision of Beta")
  
  # tau kappa
  tau.kappa.stan <- data.frame(tau = stan.marginals$tau_kappa_draws)
  tau.kappa.inlabru <- data.frame(res.inlabru$marginals.hyperpar$`Precision for kappa`)
  
  p.tau.kappa <- ggplot() + 
    geom_area(data = tau.kappa.inlabru, aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = tau.kappa.stan, aes(x = tau, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins = 100) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of precision of kappa", y = " ", title = "Precision of Kappa")
  
  if (cohort){
    # tau gamma
    tau.gamma.stan <- data.frame(tau = stan.marginals$tau_gamma_draws)
    tau.gamma.inlabru <- data.frame(res.inlabru$marginals.hyperpar$`Precision for gamma`)
    
    p.tau.gamma <- ggplot() + 
      geom_area(data = tau.gamma.inlabru, aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
      geom_histogram(data = tau.gamma.stan, aes(x = tau, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins = 100) + 
      scale_color_manual(name = " ", values = palette) + 
      scale_fill_manual(name = " ", values = palette) +
      labs(x = "Value of precision of gamma", y = " ", title = "Precision of Gamma")
  }
  
  # tau epsilon
  tau.epsilon.stan <- data.frame(tau = stan.marginals$tau_epsilon_draws)
  tau.epsilon.inlabru <- data.frame(res.inlabru$marginals.hyperpar$`Precision for epsilon`)
  
  p.tau.epsilon <- ggplot() + 
    geom_area(data = tau.epsilon.inlabru, aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = tau.epsilon.stan, aes(x = tau, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins = 100) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of precision of epsilon", y = " ", title = "Precision of Epsilon")
  
  #   ----   if save: store figures   ----
  
  if(cohort){
    p.random.effects <- (p.intercept | p.alpha | p.beta) / (p.kappa | p.gamma) + 
      plot_layout(guides = "collect")
    
    p.hypers <- (p.tau.alpha | p.tau.beta)/(p.tau.kappa | p.tau.gamma | p.tau.epsilon) + 
      plot_layout(guides = "collect")
  } else {
    p.random.effects <- (p.intercept | p.alpha ) / (p.beta | p.kappa ) + 
      plot_layout(guides = "collect")
    
    p.hypers <- (p.tau.alpha | p.tau.beta)/(p.tau.kappa | p.tau.epsilon) + 
      plot_layout(guides = "collect")
  }
  
  print("Arriving at the place in the code where we save")
  
  if(save){
    print("saving")
    save.figure(p.random.effects, name="random_effects_compared", path=path.to.storage)
    save.figure(p.eta, name="eta_xt_compared", path=path.to.storage)
    save.figure(p.eta.2, name="eta_xt_2_compared", path=path.to.storage)
    save.figure(p.eta.x, name="eta_x_compared", path=path.to.storage)
    save.figure(p.eta.t, name="eta_t_compared", path=path.to.storage)
    save.figure(p.hypers, name="hypers_compared", path=path.to.storage)
  }
  
  plots <- list(p.intercept = p.intercept, 
                p.alpha = p.alpha, 
                p.beta = p.beta,
                p.kappa = p.kappa,
                p.eta = p.eta,
                p.eta.2 = p.eta.2,
                p.eta.t = p.eta.t,
                p.eta.x = p.eta.x,
                p.tau.alpha = p.tau.alpha,
                p.tau.beta = p.tau.beta,
                p.tau.kappa = p.tau.kappa,
                p.tau.epsilon = p.tau.epsilon,
                p.random.effects=p.random.effects,
                p.hypers = p.hypers)
  if(cohort){
    plots <- c(plots, p.gamma=p.gamma)
    plots <- c(plots, p.tau.gamma=p.tau.gamma)
  }
  
  return(plots)
}

