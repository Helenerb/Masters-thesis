library(ggplot2)
library(patchwork)

#   ----   Source relevant functions
#source("../Functions/plotters.R")
#source("../Misc/palette.R")

# assume working directory at ../Master Thesis Code
source("Scripts/Functions/plotters.R")
source("Scripts/Misc/palette.R")

plot.inlabru.vs.underlying.cohort <- function(res.inlabru, underlying.effects,
                                              path.to.storage="", save=FALSE,
                                              phi.plus.kappa.func=phi.plus.kappa.v17){
  
  obs <- underlying.effects$obs
  nx <- underlying.effects$nx
  nt <- underlying.effects$nt
  
  results.period <- plot.posterior.period.effects(res.inlabru, underlying.effects, phi.plus.kappa.func)
  p.period <- results.period$posterior.plot
  
  data.alpha = cbind(res.inlabru$summary.random$alpha,
                     alpha.true = underlying.effects$alpha.true[res.inlabru$summary.random$alpha$ID + 1])
  p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = alpha.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Alpha - inlabru", x = "x", y='')
  
  data.beta = cbind(res.inlabru$summary.random$beta,
                    beta.true = underlying.effects$beta.true[res.inlabru$summary.random$beta$ID + 1])
  p.beta <- ggplot(data = data.beta, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = beta.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "x", y = "beta", title = "Beta - inlabru")
  
  data.kappa = cbind(res.inlabru$summary.random$kappa,
                     kappa.true = underlying.effects$kappa.true[res.inlabru$summary.random$kappa$ID + 1])
  p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = kappa.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "t", y = "kappa", title = "Kappa - inlabru")
  
  data.gamma = cbind(res.inlabru$summary.random$gamma,
                     gamma.true = underlying.effects$gamma.true[res.inlabru$summary.random$gamma$ID + nx])
  p.gamma <- ggplot(data = data.gamma, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = gamma.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "t", y = "gamma", title = "Gamma - inlabru")
  
  data.fixed <- data.frame(res.inlabru$marginals.fixed)
  
  p.phi <- ggplot(data.fixed) + 
    geom_area(aes(x = phi.x, y = phi.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[2], color = "Estimated", fill="Estimated")) + 
    geom_vline(aes(xintercept = underlying.effects$phi.true, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of phi", y = " ", title = "Phi - inlabru")
  
  p.intercept <- ggplot(data.fixed) + 
    geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Estimated", fill="Estimated")) + 
    geom_vline(aes(xintercept = underlying.effects$age.intercept.true, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of phi", y = " ", title = "Intercept - inlabru")
  
  p.phi.kappa <- (p.phi | p.kappa) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.phi.kappa, name = "phi_kappa_inlabru", path = path.to.storage)
  }
  
  p.random.effects <- (p.intercept | p.alpha | p.beta)/(p.period | p.gamma) + 
    plot_layout(guides = "collect")
  
  
  if(save){
    save.figure(p.random.effects, name = "random_effects_inlabru", path = path.to.storage)
  }
  
  
  data.eta <- data.frame(eta.sim = res.inlabru$summary.linear.predictor$mean[1:length(obs$eta)]) %>%
    mutate(true.eta = obs$eta) %>%
    mutate(xt = obs$xt, x = obs$x, t = obs$t)
  
  p.eta <- ggplot(data = data.eta) +
    geom_point(aes(x = eta.sim, y = true.eta), color = palette[1]) + 
    labs(x="Estimated eta", y="True value for eta", title = "Eta")
  
  p.eta.2 <- ggplot(data = data.eta) +
    geom_line(aes(x=xt, y = eta.sim, color="Estimated")) +
    geom_line(aes(x=xt, y = true.eta, color="True")) +
    labs(x=" ", y="Eta", title="Eta- inlabru")
  
  p.eta.t <- ggplot(data = data.eta) + 
    geom_line(aes(x = x, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = x, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.eta.x <- ggplot(data = data.eta) + 
    geom_line(aes(x = t, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = t, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each age") + 
    facet_wrap(~x)
  
  p.eta.xt <- (p.eta | p.eta.2) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.xt, name="eta_xt_inlabru", path=path.to.storage)
  }
  
  p.eta.facet <- (p.eta.t | p.eta.x) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.facet, name="eta_facet_inlabru", path=path.to.storage)
  }
  
  plots <- list(p.alpha = p.alpha, p.beta = p.beta, p.kappa = p.kappa,
                p.phi = p.phi, p.eta = p.eta, p.gamma = p.gamma,
                p.eta.2 = p.eta.2, p.eta.x = p.eta.x, 
                p.eta.t = p.eta.t, p.intercept = p.intercept,
                p.random.effects = p.random.effects,
                p.phi.kappa = p.phi.kappa,
                p.eta.all = p.eta.xt,
                p.eta.facet = p.eta.facet)
  
  summaries <- list(data.alpha = data.alpha,
                data.beta = data.beta,
                data.kappa = data.kappa,
                data.gamma = data.gamma,
                data.eta = data.eta,
                data.fixed = data.fixed,
                data.period = results.period$posterior.data,
                intercept = res.inlabru$summary.fixed$mean[1],
                phi = res.inlabru$summary.fixed$mean[2])
  return(list(plots = plots, summaries = summaries))
}

plot.inlabru.vs.underlying.lc <- function(res.inlabru, underlying.effects,
                                              path.to.storage="", save=FALSE,
                                              phi.plus.kappa.func=phi.plus.kappa.v17){
  
  obs <- underlying.effects$obs
  nx <- underlying.effects$nx
  nt <- underlying.effects$nt
  
  results.period <- plot.posterior.period.effects(res.inlabru, underlying.effects, phi.plus.kappa.func)
  p.period <- results.period$period.plot
  
  data.alpha = cbind(res.inlabru$summary.random$alpha,
                     alpha.true = underlying.effects$alpha.true[res.inlabru$summary.random$alpha$ID + 1])
  p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = alpha.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Alpha - inlabru", x = "x", y='')
  
  data.beta = cbind(res.inlabru$summary.random$beta,
                    beta.true = underlying.effects$beta.true[res.inlabru$summary.random$beta$ID + 1])
  p.beta <- ggplot(data = data.beta, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = beta.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "x", y = "beta", title = "Beta - inlabru")
  
  data.kappa = cbind(res.inlabru$summary.random$kappa,
                     kappa.true = underlying.effects$kappa.true[res.inlabru$summary.random$kappa$ID + 1])
  p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = kappa.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "t", y = "kappa", title = "Kappa - inlabru")
  
  data.fixed  <- data.frame(res.inlabru$marginals.fixed)
  
  p.phi <- ggplot(data.fixed) + 
    geom_area(aes(x = phi.x, y = phi.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[2], color = "Estimated", fill="Estimated")) + 
    geom_vline(aes(xintercept = underlying.effects$phi.true, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of phi", y = " ", title = "Phi - inlabru")
  
  p.intercept <- ggplot(data.fixed) + 
    geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Estimated", fill="Estimated")) + 
    geom_vline(aes(xintercept = underlying.effects$age.intercept.true, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of phi", y = " ", title = "Intercept - inlabru")
  
  p.random.effects <- (p.intercept | p.alpha ) / (p.beta | p.period ) + 
    plot_layout(guides = "collect")
  
  if(save){
    save.figure(p.random.effects, name = "random_effects_inlabru", path = path.to.storage)
  }
  
  p.phi.kappa <- (p.phi | p.kappa) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.phi.kappa, name = "phi_kappa_inlabru", path = path.to.storage)
  }
  
  data.eta <- data.frame(eta.sim = res.inlabru$summary.linear.predictor$mean[1:length(obs$eta)]) %>%
    mutate(true.eta = obs$eta) %>%
    mutate(xt = obs$xt, x = obs$x, t = obs$t)
  
  p.eta <- ggplot(data = data.eta) +
    geom_point(aes(x = eta.sim, y = true.eta), color = palette[1]) + 
    labs(x="Estimated eta", y="True value for eta", title = "Eta")
  
  p.eta.2 <- ggplot(data = data.eta) +
    geom_line(aes(x=xt, y = eta.sim, color="Estimated")) +
    geom_line(aes(x=xt, y = true.eta, color="True")) +
    labs(x=" ", y="Eta", title="Eta- inlabru")
  
  p.eta.t <- ggplot(data = data.eta) + 
    geom_line(aes(x = x, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = x, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.eta.x <- ggplot(data = data.eta) + 
    geom_line(aes(x = t, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = t, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each age") + 
    facet_wrap(~x)
  
  p.eta.xt <- (p.eta | p.eta.2) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.xt, name="eta_xt_inlabru", path=path.to.storage)
  }
  
  p.eta.facet <- (p.eta.t | p.eta.x) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.facet, name="eta_facet_inlabru", path=path.to.storage)
  }
  
  plots <- list(p.alpha = p.alpha, p.beta = p.beta, p.kappa = p.kappa,
                p.phi = p.phi, p.eta = p.eta,
                p.eta.2 = p.eta.2, p.eta.x = p.eta.x, 
                p.eta.t = p.eta.t, p.intercept = p.intercept,
                p.random.effects = p.random.effects,
                p.phi.kappa = p.phi.kappa,
                p.eta.all = p.eta.xt,
                p.eta.facet = p.eta.facet)
  
  summaries <- list(data.alpha = data.alpha,
                data.beta = data.beta,
                data.kappa = data.kappa,
                data.eta = data.eta,
                data.fixed = data.fixed,
                data.period = results.period$posterior.data,
                intercept = res.inlabru$summary.fixed$mean[1],
                phi = res.inlabru$summary.fixed$mean[2]
                )
  return(list(plots=plots, summaries=summaries))
}

plot.inlabru.vs.underlying.synthetic.cancer<- function(
  res.inlabru, underlying.effects, path.to.storage="",
  cutoff_alpha = 1000, cutoff_beta = 1000, cutoff_kappa = 1000, cutoff_epsilon=1000,
  save=FALSE, pdf = TRUE, png = TRUE, cohort = FALSE){
  
  obs <- underlying.effects$obs
  nx <- underlying.effects$nx
  nt <- underlying.effects$nt
  
  data.alpha = cbind(res.inlabru$summary.random$alpha,
                     alpha.true = underlying.effects$alpha.true[res.inlabru$summary.random$alpha$ID + 1])
  p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = alpha.true, color = "True value", fill = "True value"), size = 1) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 1) + 
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(title="Alpha - inlabru", x = "x", y='')
  
  data.beta = cbind(res.inlabru$summary.random$beta,
                    beta.true = underlying.effects$beta.true[res.inlabru$summary.random$beta$ID + 1])
  print(data.beta)
  p.beta <- ggplot(data = data.beta, aes(x = ID)) + 
    geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated", color = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = beta.true, color = "True value", fill = "True value"), size = 1) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 1) + 
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(x = "x", y = "beta", title = "Beta - inlabru")
  
  data.kappa = cbind(res.inlabru$summary.random$kappa,
                     kappa.true = underlying.effects$kappa.true[res.inlabru$summary.random$kappa$ID + 1])
  p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = kappa.true, color = "True value", fill = "True value"), size = 1) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 1) + 
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(x = "t", y = "kappa", title = "Kappa - inlabru")
  
  data.fixed  <- data.frame(res.inlabru$marginals.fixed)
  
  p.intercept <- ggplot(data.fixed) + 
    geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Estimated", fill="Estimated")) + 
    geom_vline(aes(xintercept = underlying.effects$intercept, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    theme_classic() + 
    labs(x = "Value of phi", y = " ", title = "Intercept - inlabru")
  
  p.random.effects <- (p.intercept | p.alpha ) / (p.beta | p.kappa ) + 
    plot_layout(guides = "collect")
  
  if(save){
    save.figure(p.random.effects, name = "random_effects_inlabru", path = path.to.storage,
                pdf=pdf, png=png)
  }
  
  data.eta <- data.frame(eta.sim = res.inlabru$summary.linear.predictor$mean[1:length(obs$eta)],
                         `0.025quant` = res.inlabru$summary.linear.predictor$`0.025quant`[1:length(obs$eta)],
                         `0.975quant` = res.inlabru$summary.linear.predictor$`0.975quant`[1:length(obs$eta)]) %>%
    mutate(true.eta = obs$eta) %>%
    mutate(xt = obs$xt, x = obs$x, t = obs$t)
  
  p.eta <- ggplot(data = data.eta) +
    geom_point(aes(x = eta.sim, y = true.eta), color = palette[1]) + 
    labs(x="Estimated eta", y="True value for eta", title = "Eta")
  
  p.eta.2 <- ggplot(data = data.eta) +
    geom_line(aes(x=xt, y = eta.sim, color="Estimated")) +
    geom_line(aes(x=xt, y = true.eta, color="True")) +
    labs(x=" ", y="Eta", title="Eta- inlabru")
  
  p.eta.t <- ggplot(data = data.eta) + 
    geom_line(aes(x = x, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = x, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.eta.x <- ggplot(data = data.eta) + 
    geom_line(aes(x = t, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = t, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each age") + 
    facet_wrap(~x)
  
  p.eta.xt <- (p.eta | p.eta.2) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.xt, name="eta_xt_inlabru", path=path.to.storage, pdf=pdf, png=png)
  }
  
  p.eta.facet <- (p.eta.t | p.eta.x) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.facet, name="eta_facet_inlabru", path=path.to.storage, pdf=pdf, png = png)
  }
  
  data.mr <- data.frame(mr.sim = res.inlabru$summary.fitted.values$mean[1:length(obs$eta)]) %>%
    mutate(true.mr = obs$mr) %>%
    mutate(xt = obs$xt, x = obs$x, t = obs$t)
  
  p.mr <- ggplot(data = data.mr) +
    geom_point(aes(x = mr.sim, y = true.mr), color = palette[1]) + 
    labs(x="Estimated mortality rate", y="Observed mortality rate", title = "Mortality rate")
  
  p.mr.2 <- ggplot(data = data.mr) +
    geom_line(aes(x=xt, y = mr.sim, color="Estimated")) +
    geom_line(aes(x=xt, y = true.mr, color="True")) +
    labs(x=" ", y="Mortality rate", title="Mortality rate - inlabru")
  
  p.mr.t <- ggplot(data = data.mr) + 
    geom_line(aes(x = x, y = mr.sim, color = "Estimated")) +
    geom_line(aes(x = x, y = true.mr, color = "True")) +
    labs(x = " ", y = " ", title = "Mortality rate - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.mr.x <- ggplot(data = data.mr) + 
    geom_line(aes(x = t, y = mr.sim, color = "Estimated")) +
    geom_line(aes(x = t, y = true.mr, color = "True")) +
    labs(x = " ", y = " ", title = "Mortality rate - inlabru, for each age") + 
    facet_wrap(~x)
  
  p.mr.xt <- (p.mr | p.mr.2) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.mr.xt, name="mr_xt_inlabru", path=path.to.storage, pdf=pdf, png=png)
  }
  
  p.mr.facet <- (p.mr.t | p.mr.x) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.mr.facet, name="mr_facet_inlabru", path=path.to.storage, pdf = pdf, png=png)
  }
  
  data.hyperpar <- data.frame(res.inlabru$marginals.hyperpar)
  
  p.alpha.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.alpha.x < cutoff_alpha)) + 
    geom_area(aes(x = Precision.for.alpha.x, y = Precision.for.alpha.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_vline(aes(xintercept = res.inlabru$summary.hyperpar$mean[1], color = "Inlabru", fill = "Inlabru")) +
    geom_vline(aes(xintercept = underlying.effects$tau.alpha.true, color = "Observed", fill = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = " ", y = " ", title = "Precision of alpha")
  
  if(save){
    save.figure(p.alpha.prec, name="alpha_prec_inlabru", path=path.to.storage, pdf = pdf, png=png)
  }
  
  p.beta.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.beta.x < cutoff_beta)) + 
    geom_area(aes(x = Precision.for.beta.x, y = Precision.for.beta.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_vline(aes(xintercept = res.inlabru$summary.hyperpar$mean[2], color = "Inlabru", fill = "Inlabru")) + 
    geom_vline(aes(xintercept = underlying.effects$tau.beta.true, color = "Observed", fill = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = " ", y = " ", title = "Precision of beta")
  
  if(save){
    save.figure(p.beta.prec, name="beta_prec_inlabru", path=path.to.storage, pdf = pdf, png=png)
  }
  
  p.kappa.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.kappa.x < cutoff_kappa)) + 
    geom_area(aes(x = Precision.for.kappa.x, y = Precision.for.kappa.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_vline(aes(xintercept = res.inlabru$summary.hyperpar$mean[3], color = "Inlabru", fill = "Inlabru")) + 
    geom_vline(aes(xintercept = underlying.effects$tau.kappa.true, color = "Observed", fill = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = " ", y = " ", title = "Precision of kappa")
  
  if(cohort){
    p.gamma.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.gamma.x < cutoff_gamma)) + 
      geom_area(aes(x = Precision.for.gamma.x, y = Precision.for.gamma.y), color = palette[1], fill = palette[1], alpha = 0.5) + 
      geom_vline(aes(xintercept = res.inlabru$summary.hyperpar$mean[4]), color = palette[1]) + 
      labs(x = " ", y = " ", title = "Precision of gamma")
  }
  
  p.epsilon.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.epsilon.x < cutoff_epsilon)) + 
    geom_area(aes(x = Precision.for.epsilon.x, y = Precision.for.epsilon.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$tau.epsilon.true, color = "Observed", fill = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = " ", y = " ", title = "Precision of epsilon")
  
  if(save){
    save.figure(p.epsilon.prec, name="beta_epsilon_inlabru", path=path.to.storage, pdf = pdf, png=png)
  }
  
  if(cohort){
    p.hyperpars <- (p.alpha.prec | p.beta.prec | p.kappa.prec)/(p.gamma.prec | p.epsilon.prec) + plot_layout(guides = "collect")
  } else {
    p.hyperpars <- (p.alpha.prec | p.beta.prec )/(p.kappa.prec | p.epsilon.prec) + plot_layout(guides = "collect")
  }
  
  if(save){
    save.figure(p.hyperpars, name="hypers_inlabru", path=path.to.storage, pdf = pdf, png=png)
  }
  
  
  plots <- list(p.alpha = p.alpha, p.beta = p.beta, p.kappa = p.kappa,
                p.eta = p.eta,
                p.eta.2 = p.eta.2, p.eta.x = p.eta.x, 
                p.eta.t = p.eta.t, p.intercept = p.intercept,
                p.random.effects = p.random.effects,
                p.eta.all = p.eta.xt,
                p.eta.facet = p.eta.facet)
  
  summaries <- list(data.alpha = data.alpha,
                    data.beta = data.beta,
                    data.kappa = data.kappa,
                    data.eta = data.eta,
                    intercept = res.inlabru$summary.fixed$mean[1]
  )
  return(list(plots=plots, summaries=summaries))
}

plot.inlabru.vs.underlying.synthetic.cancer.fixed.effects <- function(
  res.inlabru, underlying.effects, path.to.storage="",
  cutoff_alpha = 1000, cutoff_beta = 1000, cutoff_kappa = 1000, cutoff_epsilon=1000,
  save=FALSE, pdf = TRUE, png = TRUE, cohort = FALSE){
  
  obs <- underlying.effects$obs
  nx <- underlying.effects$nx
  nt <- underlying.effects$nt
  
  data.alpha = cbind(res.inlabru$summary.random$alpha,
                     alpha.true = underlying.effects$alpha.true[res.inlabru$summary.random$alpha$ID + 1])
  p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = alpha.true, color = "True value", fill = "True value"), size = 1) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 1) + 
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(title="Alpha - inlabru", x = "x", y='')
  
  data.beta = cbind(res.inlabru$summary.random$beta,
                    beta.true = underlying.effects$beta.true[res.inlabru$summary.random$beta$ID + 1])
  print(data.beta)
  p.beta <- ggplot(data = data.beta, aes(x = ID)) + 
    geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated", color = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = beta.true, color = "True value", fill = "True value"), size = 1) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 1) + 
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(x = "x", y = "beta", title = "Beta - inlabru")
  
  data.kappa = cbind(res.inlabru$summary.random$kappa,
                     kappa.true = underlying.effects$kappa.true[res.inlabru$summary.random$kappa$ID + 1])
  p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = kappa.true, color = "True value", fill = "True value"), size = 1) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 1) + 
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(x = "t", y = "kappa", title = "Kappa - inlabru")
  
  data.fixed  <- data.frame(res.inlabru$marginals.fixed)
  
  p.intercept <- ggplot(data.fixed) + 
    geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Estimated", fill="Estimated")) + 
    geom_vline(aes(xintercept = underlying.effects$intercept, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    theme_classic() + 
    labs(x = "Value of phi", y = " ", title = "Intercept - inlabru")
  
  p.random.effects <- (p.intercept | p.alpha ) / (p.beta | p.kappa ) + 
    plot_layout(guides = "collect")
  
  if(save){
    save.figure(p.random.effects, name = "random_effects_inlabru", path = path.to.storage,
                pdf=pdf, png=png)
  }
  
  data.eta <- data.frame(eta.sim = res.inlabru$summary.linear.predictor$mean[1:length(obs$eta)],
                         `0.025quant` = res.inlabru$summary.linear.predictor$`0.025quant`[1:length(obs$eta)],
                         `0.975quant` = res.inlabru$summary.linear.predictor$`0.975quant`[1:length(obs$eta)],
                          fitted.values = res.inlabru$summary.fitted.values$mean[1:length(obs$eta)]) %>%
    mutate(true.eta = obs$eta) %>%
    mutate(xt = obs$xt, x = obs$x, t = obs$t)
  
  p.eta <- ggplot(data = data.eta) +
    geom_point(aes(x = eta.sim, y = true.eta), color = palette[1]) + 
    labs(x="Estimated eta", y="True value for eta", title = "Eta")
  
  p.eta.2 <- ggplot(data = data.eta) +
    geom_line(aes(x=xt, y = eta.sim, color="Estimated")) +
    geom_line(aes(x=xt, y = true.eta, color="True")) +
    labs(x=" ", y="Eta", title="Eta- inlabru")
  
  p.eta.t <- ggplot(data = data.eta) + 
    geom_line(aes(x = x, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = x, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.eta.x <- ggplot(data = data.eta) + 
    geom_line(aes(x = t, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = t, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each age") + 
    facet_wrap(~x)
  
  p.eta.xt <- (p.eta | p.eta.2) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.xt, name="eta_xt_inlabru", path=path.to.storage, pdf=pdf, png=png)
  }
  
  p.eta.facet <- (p.eta.t | p.eta.x) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.facet, name="eta_facet_inlabru", path=path.to.storage, pdf=pdf, png = png)
  }
  
  data.mr <- data.frame(mr.sim = res.inlabru$summary.fitted.values$mean[1:length(obs$eta)]) %>%
    mutate(true.mr = obs$mr) %>%
    mutate(xt = obs$xt, x = obs$x, t = obs$t)
  
  p.mr <- ggplot(data = data.mr) +
    geom_point(aes(x = mr.sim, y = true.mr), color = palette[1]) + 
    labs(x="Estimated mortality rate", y="Observed mortality rate", title = "Mortality rate")
  
  p.mr.2 <- ggplot(data = data.mr) +
    geom_line(aes(x=xt, y = mr.sim, color="Estimated")) +
    geom_line(aes(x=xt, y = true.mr, color="True")) +
    labs(x=" ", y="Mortality rate", title="Mortality rate - inlabru")
  
  p.mr.t <- ggplot(data = data.mr) + 
    geom_line(aes(x = x, y = mr.sim, color = "Estimated")) +
    geom_line(aes(x = x, y = true.mr, color = "True")) +
    labs(x = " ", y = " ", title = "Mortality rate - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.mr.x <- ggplot(data = data.mr) + 
    geom_line(aes(x = t, y = mr.sim, color = "Estimated")) +
    geom_line(aes(x = t, y = true.mr, color = "True")) +
    labs(x = " ", y = " ", title = "Mortality rate - inlabru, for each age") + 
    facet_wrap(~x)
  
  p.mr.xt <- (p.mr | p.mr.2) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.mr.xt, name="mr_xt_inlabru", path=path.to.storage, pdf=pdf, png=png)
  }
  
  p.mr.facet <- (p.mr.t | p.mr.x) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.mr.facet, name="mr_facet_inlabru", path=path.to.storage, pdf = pdf, png=png)
  }
  
  # data.hyperpar <- data.frame(res.inlabru$marginals.hyperpar)
  # 
  # p.alpha.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.alpha.x < cutoff_alpha)) + 
  #   geom_area(aes(x = Precision.for.alpha.x, y = Precision.for.alpha.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  #   geom_vline(aes(xintercept = res.inlabru$summary.hyperpar$mean[1], color = "Inlabru", fill = "Inlabru")) +
  #   geom_vline(aes(xintercept = underlying.effects$tau.alpha.true, color = "Observed", fill = "Observed")) + 
  #   scale_color_manual(name = " ", values = palette) + 
  #   scale_fill_manual(name = " ", values = palette) +
  #   labs(x = " ", y = " ", title = "Precision of alpha")
  
  p.alpha.prec <- ggplot() + labs(title = "Fixed effects")
  
  # if(save){
  #   save.figure(p.alpha.prec, name="alpha_prec_inlabru", path=path.to.storage, pdf = pdf, png=png)
  # }
  # 
  # p.beta.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.beta.x < cutoff_beta)) + 
  #   geom_area(aes(x = Precision.for.beta.x, y = Precision.for.beta.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  #   geom_vline(aes(xintercept = res.inlabru$summary.hyperpar$mean[2], color = "Inlabru", fill = "Inlabru")) + 
  #   geom_vline(aes(xintercept = underlying.effects$tau.beta.true, color = "Observed", fill = "Observed")) + 
  #   scale_color_manual(name = " ", values = palette) + 
  #   scale_fill_manual(name = " ", values = palette) +
  #   labs(x = " ", y = " ", title = "Precision of beta")
  
  p.beta.prec <- ggplot() + labs(title = "Fixed effects")
  
  # if(save){
  #   save.figure(p.beta.prec, name="beta_prec_inlabru", path=path.to.storage, pdf = pdf, png=png)
  # }
  # 
  # p.kappa.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.kappa.x < cutoff_kappa)) + 
  #   geom_area(aes(x = Precision.for.kappa.x, y = Precision.for.kappa.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  #   geom_vline(aes(xintercept = res.inlabru$summary.hyperpar$mean[3], color = "Inlabru", fill = "Inlabru")) + 
  #   geom_vline(aes(xintercept = underlying.effects$tau.kappa.true, color = "Observed", fill = "Observed")) + 
  #   scale_color_manual(name = " ", values = palette) + 
  #   scale_fill_manual(name = " ", values = palette) +
  #   labs(x = " ", y = " ", title = "Precision of kappa")
  
  p.kappa.prec <- ggplot() + labs(title = "Fixed effects")
  
  # if(cohort){
  #   p.gamma.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.gamma.x < cutoff_gamma)) + 
  #     geom_area(aes(x = Precision.for.gamma.x, y = Precision.for.gamma.y), color = palette[1], fill = palette[1], alpha = 0.5) + 
  #     geom_vline(aes(xintercept = res.inlabru$summary.hyperpar$mean[4]), color = palette[1]) + 
  #     labs(x = " ", y = " ", title = "Precision of gamma")
  # }
  # 
  # p.epsilon.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.epsilon.x < cutoff_epsilon)) + 
  #   geom_area(aes(x = Precision.for.epsilon.x, y = Precision.for.epsilon.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  #   geom_vline(aes(xintercept = underlying.effects$tau.epsilon.true, color = "Observed", fill = "Observed")) + 
  #   scale_color_manual(name = " ", values = palette) + 
  #   scale_fill_manual(name = " ", values = palette) +
  #   labs(x = " ", y = " ", title = "Precision of epsilon")
  # 
  # if(save){
  #   save.figure(p.epsilon.prec, name="beta_epsilon_inlabru", path=path.to.storage, pdf = pdf, png=png)
  # }
  
  p.epsilon.prec <- ggplot() + labs(title = "Fixed effects")
  
  # if(cohort){
  #   p.hyperpars <- (p.alpha.prec | p.beta.prec | p.kappa.prec)/(p.gamma.prec | p.epsilon.prec) + plot_layout(guides = "collect")
  # } else {
  #   p.hyperpars <- (p.alpha.prec | p.beta.prec )/(p.kappa.prec | p.epsilon.prec) + plot_layout(guides = "collect")
  # }
  
  # if(save){
  #   save.figure(p.hyperpars, name="hypers_inlabru", path=path.to.storage, pdf = pdf, png=png)
  # }
  
  p.hyperpars <- ggplot() + labs(title = "Fixed effects")
  
  plots <- list(p.alpha = p.alpha, p.beta = p.beta, p.kappa = p.kappa,
                p.eta = p.eta,
                p.eta.2 = p.eta.2, p.eta.x = p.eta.x, 
                p.eta.t = p.eta.t, p.intercept = p.intercept,
                p.random.effects = p.random.effects,
                p.eta.all = p.eta.xt,
                p.eta.facet = p.eta.facet)
  
  summaries <- list(data.alpha = data.alpha,
                    data.beta = data.beta,
                    data.kappa = data.kappa,
                    data.eta = data.eta,
                    intercept = res.inlabru$summary.fixed$mean[1]
  )
  return(list(plots=plots, summaries=summaries))
}

plot.inlabru.vs.underlying.traditional.lc <- function(
  res.inlabru, underlying.effects, path.to.storage="",
  cutoff_alpha = 1000, cutoff_beta = 1000, cutoff_kappa = 1000, cutoff_epsilon=1000,
  save=FALSE, pdf = TRUE, png = TRUE, cohort = FALSE){
  
  obs <- underlying.effects$obs
  nx <- underlying.effects$nx
  nt <- underlying.effects$nt
  
  data.alpha = cbind(res.inlabru$summary.random$alpha,
                     alpha.true = underlying.effects$alpha.true[res.inlabru$summary.random$alpha$ID + 1])
  p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = alpha.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Alpha - inlabru", x = "x", y='')
  
  data.beta = cbind(res.inlabru$summary.random$beta,
                    beta.true = underlying.effects$beta.true[res.inlabru$summary.random$beta$ID + 1])
  
  p.beta <- ggplot(data = data.beta, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = beta.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "x", y = "beta", title = "Beta - inlabru")
  
  data.kappa = cbind(res.inlabru$summary.random$kappa,
                     kappa.true = underlying.effects$kappa.true[res.inlabru$summary.random$kappa$ID + 1])
  p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = kappa.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "t", y = "kappa", title = "Kappa - inlabru")
  
  data.fixed  <- data.frame(res.inlabru$marginals.fixed)
  
  p.intercept <- ggplot(data.fixed) + 
    geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Estimated", fill="Estimated")) + 
    geom_vline(aes(xintercept = underlying.effects$intercept, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of phi", y = " ", title = "Intercept - inlabru")
  
  p.random.effects <- (p.intercept | p.alpha ) / (p.beta | p.kappa ) + 
    plot_layout(guides = "collect")
  
  if(save){
    save.figure(p.random.effects, name = "random_effects_inlabru", path = path.to.storage,
                pdf=pdf, png=png)
  }
  
  data.eta <- data.frame(eta.sim = res.inlabru$summary.linear.predictor$mean[1:length(obs$eta)],
                         `0.025quant` = res.inlabru$summary.linear.predictor$`0.025quant`[1:length(obs$eta)],
                         `0.975quant` = res.inlabru$summary.linear.predictor$`0.975quant`[1:length(obs$eta)]) %>%
    mutate(true.eta = obs$eta.no.error) %>%
    mutate(true.eta.error = obs$eta) %>%
    mutate(xt = obs$xt, x = obs$x, t = obs$t)
  
  p.eta <- ggplot(data = data.eta) +
    geom_point(aes(x = eta.sim, y = true.eta), color = palette[1]) + 
    labs(x="Estimated eta", y="True value for eta", title = "Eta")
  
  p.eta.2 <- ggplot(data = data.eta) +
    geom_line(aes(x=xt, y = eta.sim, color="Estimated", fill = "Estimated")) +
    geom_line(aes(x=xt, y = true.eta, color="True", fill = "True")) +
    geom_point(aes(x = xt, y = true.eta.error, color = "True, with error", fill = "True, with error"), alpha = 0.7) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x=" ", y="Eta", title="Eta- inlabru")
  
  p.eta.t <- ggplot(data = data.eta) + 
    geom_line(aes(x = x, y = eta.sim, color = "Estimated", fill = "Estimated")) +
    geom_line(aes(x = x, y = true.eta, color = "True", fill = "True")) +
    geom_point(aes(x = x, y = true.eta.error, color = "True, with error", fill = "True, with error"), alpha = 0.7) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.eta.x <- ggplot(data = data.eta) + 
    geom_line(aes(x = t, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = t, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each age") + 
    facet_wrap(~x)
  
  p.eta.xt <- (p.eta | p.eta.2) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.xt, name="eta_xt_inlabru", path=path.to.storage, pdf=pdf, png=png)
  }
  
  p.eta.facet <- (p.eta.t | p.eta.x) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.facet, name="eta_facet_inlabru", path=path.to.storage, pdf=pdf, png = png)
  }
  
  data.hyperpar <- data.frame(res.inlabru$marginals.hyperpar)
  
  p.alpha.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.alpha.x < cutoff_alpha)) + 
    geom_area(aes(x = Precision.for.alpha.x, y = Precision.for.alpha.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$tau.alpha.true, color = "Observed", fill = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = " ", y = " ", title = "Precision of alpha")
  
  p.beta.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.beta.x < cutoff_beta)) + 
    geom_area(aes(x = Precision.for.beta.x, y = Precision.for.beta.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$tau.beta.true, color = "Observed", fill = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = " ", y = " ", title = "Precision of beta")
  
  p.kappa.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.kappa.x < cutoff_kappa)) + 
    geom_area(aes(x = Precision.for.kappa.x, y = Precision.for.kappa.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$tau.kappa.true, color = "Observed", fill = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = " ", y = " ", title = "Precision of kappa")
  
  p.epsilon.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.the.Gaussian.observations.x < cutoff_epsilon)) + 
    geom_area(aes(x = Precision.for.the.Gaussian.observations.x, y = Precision.for.the.Gaussian.observations.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    #geom_vline(aes(xintercept = underlying.effects$tau.kappa.true, color = "Observed", fill = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = " ", y = " ", title = "Precision of epsilon, Gaussian observations")
  
  if(cohort){
    p.gamma.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.gamma.x < cutoff_gamma)) + 
      geom_area(aes(x = Precision.for.gamma.x, y = Precision.for.gamma.y), color = palette[1], fill = palette[1], alpha = 0.5) + 
      geom_vline(aes(xintercept = res.inlabru$summary.hyperpar$mean[4]), color = palette[1]) + 
      labs(x = " ", y = " ", title = "Precision of gamma")
  }
  
  # p.epsilon.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.epsilon.x < cutoff_epsilon)) + 
  #   geom_area(aes(x = Precision.for.epsilon.x, y = Precision.for.epsilon.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  #   geom_vline(aes(xintercept = underlying.effects$tau.epsilon.true, color = "Observed", fill = "Observed")) + 
  #   scale_color_manual(name = " ", values = palette) + 
  #   scale_fill_manual(name = " ", values = palette) +
  #   labs(x = " ", y = " ", title = "Precision of epsilon")
  # 
  # if(save){
  #   save.figure(p.epsilon.prec, name="beta_epsilon_inlabru", path=path.to.storage, pdf = pdf, png=png)
  # }
  
  if(cohort){
    p.hyperpars <- (p.alpha.prec | p.beta.prec | p.kappa.prec)/(p.gamma.prec) + plot_layout(guides = "collect")
  } else {
    p.hyperpars <- (p.alpha.prec | p.epsilon.prec)/(p.beta.prec | p.kappa.prec ) + plot_layout(guides = "collect")
  }
  
  if(save){
    save.figure(p.hyperpars, name="hypers_inlabru", path=path.to.storage, pdf = pdf, png=png)
  }
  
  
  plots <- list(p.alpha = p.alpha, p.beta = p.beta, p.kappa = p.kappa,
                p.eta = p.eta,
                p.eta.2 = p.eta.2, p.eta.x = p.eta.x, 
                p.eta.t = p.eta.t, p.intercept = p.intercept,
                p.random.effects = p.random.effects,
                p.eta.all = p.eta.xt,
                p.eta.facet = p.eta.facet)
  
  summaries <- list(data.alpha = data.alpha,
                    data.beta = data.beta,
                    data.kappa = data.kappa,
                    data.eta = data.eta,
                    intercept = res.inlabru$summary.fixed$mean[1]
  )
  return(list(plots=plots, summaries=summaries))
}

plot.inlabru.vs.underlying.traditional.lc.fixed.effects <- function(
  res.inlabru, underlying.effects, path.to.storage="",
  cutoff_alpha = 1000, cutoff_beta = 1000, cutoff_kappa = 1000, cutoff_epsilon=1000,
  save=FALSE, pdf = TRUE, png = TRUE, cohort = FALSE){
  
  obs <- underlying.effects$obs
  nx <- underlying.effects$nx
  nt <- underlying.effects$nt
  
  data.alpha = cbind(res.inlabru$summary.random$alpha,
                     alpha.true = underlying.effects$alpha.true[res.inlabru$summary.random$alpha$ID + 1])
  p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = alpha.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Alpha - inlabru", x = "x", y='')
  
  data.beta = cbind(res.inlabru$summary.random$beta,
                    beta.true = underlying.effects$beta.true[res.inlabru$summary.random$beta$ID + 1])
  
  p.beta <- ggplot(data = data.beta, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = beta.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "x", y = "beta", title = "Beta - inlabru")
  
  data.kappa = cbind(res.inlabru$summary.random$kappa,
                     kappa.true = underlying.effects$kappa.true[res.inlabru$summary.random$kappa$ID + 1])
  p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = kappa.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "t", y = "kappa", title = "Kappa - inlabru")
  
  data.fixed  <- data.frame(res.inlabru$marginals.fixed)
  
  p.intercept <- ggplot(data.fixed) + 
    geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Estimated", fill="Estimated")) + 
    geom_vline(aes(xintercept = underlying.effects$intercept, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of phi", y = " ", title = "Intercept - inlabru")
  
  p.random.effects <- (p.intercept | p.alpha ) / (p.beta | p.kappa ) + 
    plot_layout(guides = "collect")
  
  if(save){
    save.figure(p.random.effects, name = "random_effects_inlabru", path = path.to.storage,
                pdf=pdf, png=png)
  }
  
  data.eta <- data.frame(eta.sim = res.inlabru$summary.linear.predictor$mean[1:length(obs$eta)],
                         `0.025quant` = res.inlabru$summary.linear.predictor$`0.025quant`[1:length(obs$eta)],
                         `0.975quant` = res.inlabru$summary.linear.predictor$`0.975quant`[1:length(obs$eta)],
                         fitted.values = res.inlabru$summary.fitted.values$mean[1:length(obs$eta)]) %>%
    mutate(true.eta = obs$eta.no.error) %>%
    mutate(true.eta.error = obs$eta) %>%
    mutate(xt = obs$xt, x = obs$x, t = obs$t)
  
  p.eta <- ggplot(data = data.eta) +
    geom_point(aes(x = eta.sim, y = true.eta), color = palette[1]) + 
    labs(x="Estimated eta", y="True value for eta", title = "Eta")
  
  p.eta.2 <- ggplot(data = data.eta) +
    geom_point(aes(x=xt, y = eta.sim, color="Estimated", fill = "Estimated")) +
    geom_ribbon(aes(x = xt, ymin = X0.025quant, ymax = X0.975quant, fill = "Estimated"), alpha = 0.5) + 
    geom_point(aes(x=xt, y = true.eta, color="True", fill = "True")) +
    geom_point(aes(x = xt, y = true.eta.error, color = "True, with error", fill = "True, with error"), shape = 4) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x=" ", y="Eta", title="Eta- inlabru")
  
  p.eta.t <- ggplot(data = data.eta) + 
    geom_line(aes(x = x, y = eta.sim, color = "Estimated", fill = "Estimated")) +
    geom_line(aes(x = x, y = true.eta, color = "True", fill = "True")) +
    geom_point(aes(x = x, y = true.eta.error, color = "True, with error", fill = "True, with error"), alpha = 0.7) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.eta.x <- ggplot(data = data.eta) + 
    geom_line(aes(x = t, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = t, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each age") + 
    facet_wrap(~x)
  
  p.eta.xt <- (p.eta | p.eta.2) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  p.diff.eta <- ggplot(data = data.eta %>% mutate(diff.eta.inlabru = eta.sim - true.eta)) +
    geom_point(aes(x = xt, y = diff.eta.inlabru)) + 
    labs(title = "Difference of true eta and eta simulated by inlabru")
  
  if(save){
    save.figure(p.diff.eta, name="eta_diff", path=path.to.storage, pdf=pdf, png=png)
  }
  
  if(save){
    save.figure(p.eta.xt, name="eta_xt_inlabru", path=path.to.storage, pdf=pdf, png=png)
  }
  
  p.eta.facet <- (p.eta.t | p.eta.x) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.facet, name="eta_facet_inlabru", path=path.to.storage, pdf=pdf, png = png)
  }
  
  data.hyperpar <- data.frame(res.inlabru$marginals.hyperpar)
  
  # p.alpha.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.alpha.x < cutoff_alpha)) + 
  #   geom_area(aes(x = Precision.for.alpha.x, y = Precision.for.alpha.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  #   geom_vline(aes(xintercept = underlying.effects$tau.alpha.true, color = "Observed", fill = "Observed")) + 
  #   scale_color_manual(name = " ", values = palette) + 
  #   scale_fill_manual(name = " ", values = palette) +
  #   labs(x = " ", y = " ", title = "Precision of alpha")
  # 
  # p.beta.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.beta.x < cutoff_beta)) + 
  #   geom_area(aes(x = Precision.for.beta.x, y = Precision.for.beta.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  #   geom_vline(aes(xintercept = underlying.effects$tau.beta.true, color = "Observed", fill = "Observed")) + 
  #   scale_color_manual(name = " ", values = palette) + 
  #   scale_fill_manual(name = " ", values = palette) +
  #   labs(x = " ", y = " ", title = "Precision of beta")
  # 
  # p.kappa.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.kappa.x < cutoff_kappa)) + 
  #   geom_area(aes(x = Precision.for.kappa.x, y = Precision.for.kappa.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  #   geom_vline(aes(xintercept = underlying.effects$tau.kappa.true, color = "Observed", fill = "Observed")) + 
  #   scale_color_manual(name = " ", values = palette) + 
  #   scale_fill_manual(name = " ", values = palette) +
  #   labs(x = " ", y = " ", title = "Precision of kappa")
  
  # if(cohort){
  #   p.gamma.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.gamma.x < cutoff_gamma)) + 
  #     geom_area(aes(x = Precision.for.gamma.x, y = Precision.for.gamma.y), color = palette[1], fill = palette[1], alpha = 0.5) + 
  #     geom_vline(aes(xintercept = res.inlabru$summary.hyperpar$mean[4]), color = palette[1]) + 
  #     labs(x = " ", y = " ", title = "Precision of gamma")
  # }
  
  # p.epsilon.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.epsilon.x < cutoff_epsilon)) + 
  #   geom_area(aes(x = Precision.for.epsilon.x, y = Precision.for.epsilon.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  #   geom_vline(aes(xintercept = underlying.effects$tau.epsilon.true, color = "Observed", fill = "Observed")) + 
  #   scale_color_manual(name = " ", values = palette) + 
  #   scale_fill_manual(name = " ", values = palette) +
  #   labs(x = " ", y = " ", title = "Precision of epsilon")
  # 
  # if(save){
  #   save.figure(p.epsilon.prec, name="beta_epsilon_inlabru", path=path.to.storage, pdf = pdf, png=png)
  # }
  
  # if(cohort){
  #   p.hyperpars <- (p.alpha.prec | p.beta.prec | p.kappa.prec)/(p.gamma.prec) + plot_layout(guides = "collect")
  # } else {
  #   p.hyperpars <- (p.alpha.prec )/(p.beta.prec | p.kappa.prec ) + plot_layout(guides = "collect")
  # }
  # 
  # if(save){
  #   save.figure(p.hyperpars, name="hypers_inlabru", path=path.to.storage, pdf = pdf, png=png)
  # }
  
  
  plots <- list(p.alpha = p.alpha, p.beta = p.beta, p.kappa = p.kappa,
                p.eta = p.eta,
                p.eta.2 = p.eta.2, p.eta.x = p.eta.x, 
                p.eta.t = p.eta.t, p.intercept = p.intercept,
                p.random.effects = p.random.effects,
                p.eta.all = p.eta.xt,
                p.eta.facet = p.eta.facet)
  
  summaries <- list(data.alpha = data.alpha,
                    data.beta = data.beta,
                    data.kappa = data.kappa,
                    data.eta = data.eta,
                    intercept = res.inlabru$summary.fixed$mean[1]
  )
  return(list(plots=plots, summaries=summaries))
}

plot.inlabru.vs.underlying.traditional.lc.fixed.effects.no.beta <- function(
  res.inlabru, underlying.effects, path.to.storage="",
  cutoff_alpha = 1000, cutoff_beta = 1000, cutoff_kappa = 1000, cutoff_epsilon=1000,
  save=FALSE, pdf = TRUE, png = TRUE, cohort = FALSE){
  
  obs <- underlying.effects$obs
  nx <- underlying.effects$nx
  nt <- underlying.effects$nt
  
  data.alpha = cbind(res.inlabru$summary.random$alpha,
                     alpha.true = underlying.effects$alpha.true[res.inlabru$summary.random$alpha$ID + 1])
  p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = alpha.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Alpha - inlabru", x = "x", y='')
  
  data.beta = cbind(res.inlabru$summary.random$beta,
                    beta.true = underlying.effects$beta.true[res.inlabru$summary.random$beta$ID + 1])
  
  # p.beta <- ggplot(data = data.beta, aes(x = ID)) + 
  #   geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
  #   geom_point(aes(y = beta.true, color = "True value", fill = "True value"), size = 0.5) + 
  #   geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
  #   scale_color_manual(name = "",
  #                      values = palette ) +
  #   scale_fill_manual(name = "",
  #                     values = palette ) +
  #   labs(x = "x", y = "beta", title = "Beta - inlabru")
  
  p.beta <- ggplot() + labs(title = "No beta")
  
  data.kappa = cbind(res.inlabru$summary.random$kappa,
                     kappa.true = underlying.effects$kappa.true[res.inlabru$summary.random$kappa$ID + 1])
  p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = kappa.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "t", y = "kappa", title = "Kappa - inlabru")
  
  data.fixed  <- data.frame(res.inlabru$marginals.fixed)
  
  p.intercept <- ggplot(data.fixed) + 
    geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Estimated", fill="Estimated")) + 
    geom_vline(aes(xintercept = underlying.effects$intercept, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of phi", y = " ", title = "Intercept - inlabru")
  
  p.random.effects <- (p.intercept | p.alpha ) / (p.beta | p.kappa ) + 
    plot_layout(guides = "collect")
  
  if(save){
    save.figure(p.random.effects, name = "random_effects_inlabru", path = path.to.storage,
                pdf=pdf, png=png)
  }
  
  data.eta <- data.frame(eta.sim = res.inlabru$summary.linear.predictor$mean[1:length(obs$eta)],
                         `0.025quant` = res.inlabru$summary.linear.predictor$`0.025quant`[1:length(obs$eta)],
                         `0.975quant` = res.inlabru$summary.linear.predictor$`0.975quant`[1:length(obs$eta)],
                         fitted.values = res.inlabru$summary.fitted.values$mean[1:length(obs$eta)]) %>%
    mutate(true.eta = obs$eta.no.error) %>%
    mutate(true.eta.error = obs$eta) %>%
    mutate(xt = obs$xt, x = obs$x, t = obs$t)
  
  p.eta <- ggplot(data = data.eta) +
    geom_point(aes(x = eta.sim, y = true.eta), color = palette[1]) + 
    labs(x="Estimated eta", y="True value for eta", title = "Eta")
  
  p.eta.2 <- ggplot(data = data.eta) +
    geom_line(aes(x=xt, y = eta.sim, color="Estimated", fill = "Estimated")) +
    geom_line(aes(x=xt, y = true.eta, color="True", fill = "True")) +
    geom_point(aes(x = xt, y = true.eta.error, color = "True, with error", fill = "True, with error"), alpha = 0.7) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x=" ", y="Eta", title="Eta- inlabru")
  
  p.eta.t <- ggplot(data = data.eta) + 
    geom_line(aes(x = x, y = eta.sim, color = "Estimated", fill = "Estimated")) +
    geom_line(aes(x = x, y = true.eta, color = "True", fill = "True")) +
    geom_point(aes(x = x, y = true.eta.error, color = "True, with error", fill = "True, with error"), alpha = 0.7) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.eta.x <- ggplot(data = data.eta) + 
    geom_line(aes(x = t, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = t, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each age") + 
    facet_wrap(~x)
  
  p.eta.xt <- (p.eta | p.eta.2) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.xt, name="eta_xt_inlabru", path=path.to.storage, pdf=pdf, png=png)
  }
  
  p.eta.facet <- (p.eta.t | p.eta.x) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.facet, name="eta_facet_inlabru", path=path.to.storage, pdf=pdf, png = png)
  }
  
  data.hyperpar <- data.frame(res.inlabru$marginals.hyperpar)
  
  # p.alpha.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.alpha.x < cutoff_alpha)) + 
  #   geom_area(aes(x = Precision.for.alpha.x, y = Precision.for.alpha.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  #   geom_vline(aes(xintercept = underlying.effects$tau.alpha.true, color = "Observed", fill = "Observed")) + 
  #   scale_color_manual(name = " ", values = palette) + 
  #   scale_fill_manual(name = " ", values = palette) +
  #   labs(x = " ", y = " ", title = "Precision of alpha")
  # 
  # p.beta.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.beta.x < cutoff_beta)) + 
  #   geom_area(aes(x = Precision.for.beta.x, y = Precision.for.beta.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  #   geom_vline(aes(xintercept = underlying.effects$tau.beta.true, color = "Observed", fill = "Observed")) + 
  #   scale_color_manual(name = " ", values = palette) + 
  #   scale_fill_manual(name = " ", values = palette) +
  #   labs(x = " ", y = " ", title = "Precision of beta")
  # 
  # p.kappa.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.kappa.x < cutoff_kappa)) + 
  #   geom_area(aes(x = Precision.for.kappa.x, y = Precision.for.kappa.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  #   geom_vline(aes(xintercept = underlying.effects$tau.kappa.true, color = "Observed", fill = "Observed")) + 
  #   scale_color_manual(name = " ", values = palette) + 
  #   scale_fill_manual(name = " ", values = palette) +
  #   labs(x = " ", y = " ", title = "Precision of kappa")
  
  # if(cohort){
  #   p.gamma.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.gamma.x < cutoff_gamma)) + 
  #     geom_area(aes(x = Precision.for.gamma.x, y = Precision.for.gamma.y), color = palette[1], fill = palette[1], alpha = 0.5) + 
  #     geom_vline(aes(xintercept = res.inlabru$summary.hyperpar$mean[4]), color = palette[1]) + 
  #     labs(x = " ", y = " ", title = "Precision of gamma")
  # }
  
  # p.epsilon.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.epsilon.x < cutoff_epsilon)) + 
  #   geom_area(aes(x = Precision.for.epsilon.x, y = Precision.for.epsilon.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  #   geom_vline(aes(xintercept = underlying.effects$tau.epsilon.true, color = "Observed", fill = "Observed")) + 
  #   scale_color_manual(name = " ", values = palette) + 
  #   scale_fill_manual(name = " ", values = palette) +
  #   labs(x = " ", y = " ", title = "Precision of epsilon")
  # 
  # if(save){
  #   save.figure(p.epsilon.prec, name="beta_epsilon_inlabru", path=path.to.storage, pdf = pdf, png=png)
  # }
  
  # if(cohort){
  #   p.hyperpars <- (p.alpha.prec | p.beta.prec | p.kappa.prec)/(p.gamma.prec) + plot_layout(guides = "collect")
  # } else {
  #   p.hyperpars <- (p.alpha.prec )/(p.beta.prec | p.kappa.prec ) + plot_layout(guides = "collect")
  # }
  # 
  # if(save){
  #   save.figure(p.hyperpars, name="hypers_inlabru", path=path.to.storage, pdf = pdf, png=png)
  # }
  
  
  plots <- list(p.alpha = p.alpha, p.beta = p.beta, p.kappa = p.kappa,
                p.eta = p.eta,
                p.eta.2 = p.eta.2, p.eta.x = p.eta.x, 
                p.eta.t = p.eta.t, p.intercept = p.intercept,
                p.random.effects = p.random.effects,
                p.eta.all = p.eta.xt,
                p.eta.facet = p.eta.facet)
  
  summaries <- list(data.alpha = data.alpha,
                    data.beta = data.beta,
                    data.kappa = data.kappa,
                    data.eta = data.eta,
                    intercept = res.inlabru$summary.fixed$mean[1]
  )
  return(list(plots=plots, summaries=summaries))
}

plot.inlabru.vs.underlying.traditional.lc.no.beta <- function(
  res.inlabru, underlying.effects, path.to.storage="",
  cutoff_alpha = 1000, cutoff_beta = 1000, cutoff_kappa = 1000, cutoff_epsilon=1000,
  save=FALSE, pdf = TRUE, png = TRUE, cohort = FALSE){
  
  obs <- underlying.effects$obs
  nx <- underlying.effects$nx
  nt <- underlying.effects$nt
  
  data.alpha = cbind(res.inlabru$summary.random$alpha,
                     alpha.true = underlying.effects$alpha.true[res.inlabru$summary.random$alpha$ID + 1])
  p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = alpha.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Alpha - inlabru", x = "x", y='')
  
  # data.beta = cbind(res.inlabru$summary.random$beta,
  #                   beta.true = underlying.effects$beta.true[res.inlabru$summary.random$beta$ID + 1])
  # 
  # p.beta <- ggplot(data = data.beta, aes(x = ID)) + 
  #   geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
  #   geom_point(aes(y = beta.true, color = "True value", fill = "True value"), size = 0.5) + 
  #   geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
  #   scale_color_manual(name = "",
  #                      values = palette ) +
  #   scale_fill_manual(name = "",
  #                     values = palette ) +
  #   labs(x = "x", y = "beta", title = "Beta - inlabru")
  # 
  data.kappa = cbind(res.inlabru$summary.random$kappa,
                     kappa.true = underlying.effects$kappa.true[res.inlabru$summary.random$kappa$ID + 1])
  p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = kappa.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "t", y = "kappa", title = "Kappa - inlabru")
  
  data.fixed  <- data.frame(res.inlabru$marginals.fixed)
  
  p.intercept <- ggplot(data.fixed) + 
    geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Estimated", fill="Estimated")) + 
    geom_vline(aes(xintercept = underlying.effects$intercept, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of phi", y = " ", title = "Intercept - inlabru")
  
  p.random.effects <- (p.intercept | p.alpha ) / (p.kappa ) + 
    plot_layout(guides = "collect")
  
  if(save){
    save.figure(p.random.effects, name = "random_effects_inlabru", path = path.to.storage,
                pdf=pdf, png=png)
  }
  
  data.eta <- data.frame(eta.sim = res.inlabru$summary.linear.predictor$mean[1:length(obs$eta)],
                         `0.025quant` = res.inlabru$summary.linear.predictor$`0.025quant`[1:length(obs$eta)],
                         `0.975quant` = res.inlabru$summary.linear.predictor$`0.975quant`[1:length(obs$eta)]) %>%
    mutate(true.eta = obs$eta.no.error) %>%
    mutate(true.eta.error = obs$eta) %>%
    mutate(xt = obs$xt, x = obs$x, t = obs$t)
  
  p.eta <- ggplot(data = data.eta) +
    geom_point(aes(x = eta.sim, y = true.eta), color = palette[1]) + 
    labs(x="Estimated eta", y="True value for eta", title = "Eta")
  
  p.eta.2 <- ggplot(data = data.eta) +
    geom_line(aes(x=xt, y = eta.sim, color="Estimated", fill = "Estimated")) +
    geom_line(aes(x=xt, y = true.eta, color="True", fill = "True")) +
    geom_point(aes(x = xt, y = true.eta.error, color = "True, with error", fill = "True, with error"), alpha = 0.7) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x=" ", y="Eta", title="Eta- inlabru")
  
  p.eta.t <- ggplot(data = data.eta) + 
    geom_line(aes(x = x, y = eta.sim, color = "Estimated", fill = "Estimated")) +
    geom_line(aes(x = x, y = true.eta, color = "True", fill = "True")) +
    geom_point(aes(x = x, y = true.eta.error, color = "True, with error", fill = "True, with error"), alpha = 0.7) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.eta.x <- ggplot(data = data.eta) + 
    geom_line(aes(x = t, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = t, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each age") + 
    facet_wrap(~x)
  
  p.eta.xt <- (p.eta | p.eta.2) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.xt, name="eta_xt_inlabru", path=path.to.storage, pdf=pdf, png=png)
  }
  
  p.eta.facet <- (p.eta.t | p.eta.x) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.facet, name="eta_facet_inlabru", path=path.to.storage, pdf=pdf, png = png)
  }
  
  data.hyperpar <- data.frame(res.inlabru$marginals.hyperpar)
  
  p.alpha.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.alpha.x < cutoff_alpha)) + 
    geom_area(aes(x = Precision.for.alpha.x, y = Precision.for.alpha.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$tau.alpha.true, color = "Observed", fill = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = " ", y = " ", title = "Precision of alpha")
  
  # p.beta.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.beta.x < cutoff_beta)) + 
  #   geom_area(aes(x = Precision.for.beta.x, y = Precision.for.beta.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  #   geom_vline(aes(xintercept = underlying.effects$tau.beta.true, color = "Observed", fill = "Observed")) + 
  #   scale_color_manual(name = " ", values = palette) + 
  #   scale_fill_manual(name = " ", values = palette) +
  #   labs(x = " ", y = " ", title = "Precision of beta")
  
  p.kappa.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.kappa.x < cutoff_kappa)) + 
    geom_area(aes(x = Precision.for.kappa.x, y = Precision.for.kappa.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$tau.kappa.true, color = "Observed", fill = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = " ", y = " ", title = "Precision of kappa")
  
  if(cohort){
    p.gamma.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.gamma.x < cutoff_gamma)) + 
      geom_area(aes(x = Precision.for.gamma.x, y = Precision.for.gamma.y), color = palette[1], fill = palette[1], alpha = 0.5) + 
      geom_vline(aes(xintercept = res.inlabru$summary.hyperpar$mean[4]), color = palette[1]) + 
      labs(x = " ", y = " ", title = "Precision of gamma")
  }
  
  # p.epsilon.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.epsilon.x < cutoff_epsilon)) + 
  #   geom_area(aes(x = Precision.for.epsilon.x, y = Precision.for.epsilon.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  #   geom_vline(aes(xintercept = underlying.effects$tau.epsilon.true, color = "Observed", fill = "Observed")) + 
  #   scale_color_manual(name = " ", values = palette) + 
  #   scale_fill_manual(name = " ", values = palette) +
  #   labs(x = " ", y = " ", title = "Precision of epsilon")
  # 
  # if(save){
  #   save.figure(p.epsilon.prec, name="beta_epsilon_inlabru", path=path.to.storage, pdf = pdf, png=png)
  # }
  
  if(cohort){
    p.hyperpars <- (p.alpha.prec | p.kappa.prec)/(p.gamma.prec) + plot_layout(guides = "collect")
  } else {
    p.hyperpars <- (p.alpha.prec | p.kappa.prec ) + plot_layout(guides = "collect")
  }
  
  if(save){
    save.figure(p.hyperpars, name="hypers_inlabru", path=path.to.storage, pdf = pdf, png=png)
  }
  
  
  plots <- list(p.alpha = p.alpha, p.kappa = p.kappa,
                p.eta = p.eta,
                p.eta.2 = p.eta.2, p.eta.x = p.eta.x, 
                p.eta.t = p.eta.t, p.intercept = p.intercept,
                p.random.effects = p.random.effects,
                p.eta.all = p.eta.xt,
                p.eta.facet = p.eta.facet)
  
  summaries <- list(data.alpha = data.alpha,
                    data.kappa = data.kappa,
                    data.eta = data.eta,
                    intercept = res.inlabru$summary.fixed$mean[1]
  )
  return(list(plots=plots, summaries=summaries))
}

plot.inlabru.vs.underlying.poisson.lc.no.beta <- function(
  res.inlabru, underlying.effects, path.to.storage="",
  cutoff_alpha = 1000, cutoff_beta = 1000, cutoff_kappa = 1000, cutoff_epsilon=1000,
  save=FALSE, pdf = TRUE, png = TRUE){
  
  obs <- underlying.effects$obs
  nx <- underlying.effects$nx
  nt <- underlying.effects$nt
  
  data.alpha = cbind(res.inlabru$summary.random$alpha,
                     alpha.true = underlying.effects$alpha.true[res.inlabru$summary.random$alpha$ID + 1])
  p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = alpha.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Alpha - inlabru", x = "x", y='')
  
  data.kappa = cbind(res.inlabru$summary.random$kappa,
                     kappa.true = underlying.effects$kappa.true[res.inlabru$summary.random$kappa$ID + 1])
  p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = kappa.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "t", y = "kappa", title = "Kappa - inlabru")
  
  data.fixed  <- data.frame(res.inlabru$marginals.fixed)
  
  p.intercept <- ggplot(data.fixed) + 
    geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Estimated", fill="Estimated")) + 
    geom_vline(aes(xintercept = underlying.effects$intercept, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of phi", y = " ", title = "Intercept - inlabru")
  
  p.random.effects <- (p.intercept | p.alpha ) / (p.kappa ) + 
    plot_layout(guides = "collect")
  
  if(save){
    save.figure(p.random.effects, name = "random_effects_inlabru", path = path.to.storage,
                pdf=pdf, png=png)
  }
  
  data.eta <- data.frame(eta.sim = res.inlabru$summary.linear.predictor$mean[1:length(obs$eta)],
                         `0.025quant` = res.inlabru$summary.linear.predictor$`0.025quant`[1:length(obs$eta)],
                         `0.975quant` = res.inlabru$summary.linear.predictor$`0.975quant`[1:length(obs$eta)]) %>%
    mutate(true.eta = obs$eta) %>%
    mutate(xt = obs$xt, x = obs$x, t = obs$t)
  
  p.eta <- ggplot(data = data.eta) +
    geom_point(aes(x = eta.sim, y = true.eta), color = palette[1]) + 
    labs(x="Estimated eta", y="True value for eta", title = "Eta")
  
  p.eta.2 <- ggplot(data = data.eta) +
    geom_line(aes(x=xt, y = eta.sim, color="Estimated", fill = "Estimated")) +
    geom_line(aes(x=xt, y = true.eta, color="True", fill = "True")) +
    geom_point(aes(x = xt, y = true.eta.error, color = "True, with error", fill = "True, with error"), alpha = 0.7) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x=" ", y="Eta", title="Eta- inlabru")
  
  p.eta.t <- ggplot(data = data.eta) + 
    geom_line(aes(x = x, y = eta.sim, color = "Estimated", fill = "Estimated")) +
    geom_line(aes(x = x, y = true.eta, color = "True", fill = "True")) +
    geom_point(aes(x = x, y = true.eta.error, color = "True, with error", fill = "True, with error"), alpha = 0.7) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.eta.x <- ggplot(data = data.eta) + 
    geom_line(aes(x = t, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = t, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each age") + 
    facet_wrap(~x)
  
  p.eta.xt <- (p.eta | p.eta.2) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.xt, name="eta_xt_inlabru", path=path.to.storage, pdf=pdf, png=png)
  }
  
  p.eta.facet <- (p.eta.t | p.eta.x) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.facet, name="eta_facet_inlabru", path=path.to.storage, pdf=pdf, png = png)
  }
  
  data.hyperpar <- data.frame(res.inlabru$marginals.hyperpar)
  
  p.alpha.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.alpha.x < cutoff_alpha)) + 
    geom_area(aes(x = Precision.for.alpha.x, y = Precision.for.alpha.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$tau.alpha.true, color = "Observed", fill = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = " ", y = " ", title = "Precision of alpha")
  
  p.kappa.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.kappa.x < cutoff_kappa)) + 
    geom_area(aes(x = Precision.for.kappa.x, y = Precision.for.kappa.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$tau.kappa.true, color = "Observed", fill = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = " ", y = " ", title = "Precision of kappa")
  
  if(save){
    save.figure(p.hyperpars, name="hypers_inlabru", path=path.to.storage, pdf = pdf, png=png)
  }
  
  
  plots <- list(p.alpha = p.alpha, p.kappa = p.kappa,
                p.eta = p.eta,
                p.eta.2 = p.eta.2, p.eta.x = p.eta.x, 
                p.eta.t = p.eta.t, p.intercept = p.intercept,
                p.random.effects = p.random.effects,
                p.eta.all = p.eta.xt,
                p.eta.facet = p.eta.facet)
  
  summaries <- list(data.alpha = data.alpha,
                    data.kappa = data.kappa,
                    data.eta = data.eta,
                    intercept = res.inlabru$summary.fixed$mean[1]
  )
  return(list(plots=plots, summaries=summaries))
}

plot.inlabru.vs.underlying.poisson.lc.fixed.hypers.no.beta <- function(
  res.inlabru, underlying.effects, path.to.storage="",
  cutoff_alpha = 1000, cutoff_beta = 1000, cutoff_kappa = 1000, cutoff_epsilon=1000,
  save=FALSE, pdf = TRUE, png = TRUE){
  #' Note: this functions needs altering if it should be useful when save=TRUE
  
  obs <- underlying.effects$obs
  nx <- underlying.effects$nx
  nt <- underlying.effects$nt
  
  data.alpha = cbind(res.inlabru$summary.random$alpha,
                     alpha.true = underlying.effects$alpha.true[res.inlabru$summary.random$alpha$ID + 1])
  p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = alpha.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Alpha - inlabru", x = "x", y='')
  
  data.kappa = cbind(res.inlabru$summary.random$kappa,
                     kappa.true = underlying.effects$kappa.true[res.inlabru$summary.random$kappa$ID + 1])
  p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = kappa.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "t", y = "kappa", title = "Kappa - inlabru")
  
  data.fixed  <- data.frame(res.inlabru$marginals.fixed)
  
  p.intercept <- ggplot(data.fixed) + 
    geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Estimated", fill="Estimated")) + 
    geom_vline(aes(xintercept = underlying.effects$intercept, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of phi", y = " ", title = "Intercept - inlabru")
  
  p.random.effects <- (p.intercept | p.alpha ) / (p.kappa ) + 
    plot_layout(guides = "collect")
  
  if(save){
    save.figure(p.random.effects, name = "random_effects_inlabru", path = path.to.storage,
                pdf=pdf, png=png)
  }
  
  data.eta <- data.frame(eta.sim = res.inlabru$summary.linear.predictor$mean[1:length(obs$eta)],
                         `0.025quant` = res.inlabru$summary.linear.predictor$`0.025quant`[1:length(obs$eta)],
                         `0.975quant` = res.inlabru$summary.linear.predictor$`0.975quant`[1:length(obs$eta)]) %>%
    mutate(true.eta = obs$eta) %>%
    mutate(xt = obs$xt, x = obs$x, t = obs$t)
  
  p.eta <- ggplot(data = data.eta) +
    geom_point(aes(x = eta.sim, y = true.eta), color = palette[1]) + 
    labs(x="Estimated eta", y="True value for eta", title = "Eta")
  
  p.eta.2 <- ggplot(data = data.eta) +
    geom_line(aes(x=xt, y = eta.sim, color="Estimated", fill = "Estimated")) +
    geom_line(aes(x=xt, y = true.eta, color="True", fill = "True")) +
    geom_point(aes(x = xt, y = true.eta.error, color = "True, with error", fill = "True, with error"), alpha = 0.7) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x=" ", y="Eta", title="Eta- inlabru")
  
  p.eta.t <- ggplot(data = data.eta) + 
    geom_line(aes(x = x, y = eta.sim, color = "Estimated", fill = "Estimated")) +
    geom_line(aes(x = x, y = true.eta, color = "True", fill = "True")) +
    geom_point(aes(x = x, y = true.eta.error, color = "True, with error", fill = "True, with error"), alpha = 0.7) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.eta.x <- ggplot(data = data.eta) + 
    geom_line(aes(x = t, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = t, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each age") + 
    facet_wrap(~x)
  
  p.eta.xt <- (p.eta | p.eta.2) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.xt, name="eta_xt_inlabru", path=path.to.storage, pdf=pdf, png=png)
  }
  
  p.eta.facet <- (p.eta.t | p.eta.x) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.facet, name="eta_facet_inlabru", path=path.to.storage, pdf=pdf, png = png)
  }
  
  data.hyperpar <- data.frame(res.inlabru$marginals.hyperpar)
  
  p.alpha.prec <- ggplot() + labs(title="Fixed hypers")
  p.kappa.prec <- ggplot() + labs(title="Fixed hypers")
  
  # p.alpha.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.alpha.x < cutoff_alpha)) + 
  #   geom_area(aes(x = Precision.for.alpha.x, y = Precision.for.alpha.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  #   geom_vline(aes(xintercept = underlying.effects$tau.alpha.true, color = "Observed", fill = "Observed")) + 
  #   scale_color_manual(name = " ", values = palette) + 
  #   scale_fill_manual(name = " ", values = palette) +
  #   labs(x = " ", y = " ", title = "Precision of alpha")
  
  # p.kappa.prec <-ggplot(data = data.hyperpar %>% filter(Precision.for.kappa.x < cutoff_kappa)) + 
  #   geom_area(aes(x = Precision.for.kappa.x, y = Precision.for.kappa.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  #   geom_vline(aes(xintercept = underlying.effects$tau.kappa.true, color = "Observed", fill = "Observed")) + 
  #   scale_color_manual(name = " ", values = palette) + 
  #   scale_fill_manual(name = " ", values = palette) +
  #   labs(x = " ", y = " ", title = "Precision of kappa")
  
  # if(save){
  #   save.figure(p.hyperpars, name="hypers_inlabru", path=path.to.storage, pdf = pdf, png=png)
  # }
  
  
  plots <- list(p.alpha = p.alpha, p.kappa = p.kappa,
                p.eta = p.eta,
                p.eta.2 = p.eta.2, p.eta.x = p.eta.x, 
                p.eta.t = p.eta.t, p.intercept = p.intercept,
                p.random.effects = p.random.effects,
                p.eta.all = p.eta.xt,
                p.eta.facet = p.eta.facet)
  
  summaries <- list(data.alpha = data.alpha,
                    data.kappa = data.kappa,
                    data.eta = data.eta,
                    intercept = res.inlabru$summary.fixed$mean[1]
  )
  return(list(plots=plots, summaries=summaries))
}

plot.inlabru.vs.underlying.lc.only.kappa.2 <- function(res.inlabru, underlying.effects,
                                                     path.to.storage="", save=FALSE,
                                                     phi.plus.kappa.func=phi.plus.kappa.v17){
  
  obs <- underlying.effects$obs
  nx <- underlying.effects$nx
  nt <- underlying.effects$nt
  
  data.alpha = cbind(res.inlabru$summary.random$alpha,
                     alpha.true = underlying.effects$alpha.true[res.inlabru$summary.random$alpha$ID + 1])
  p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = alpha.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Alpha - inlabru", x = "x", y='')
  
  data.beta = cbind(res.inlabru$summary.random$beta,
                    beta.true = underlying.effects$beta.true[res.inlabru$summary.random$beta$ID + 1])
  p.beta <- ggplot(data = data.beta, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = beta.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "x", y = "beta", title = "Beta - inlabru")
  
  data.kappa = cbind(res.inlabru$summary.random$kappa,
                     kappa.drifted = underlying.effects$kappa.drifted[res.inlabru$summary.random$kappa$ID + 1])
  p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = kappa.drifted, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "t", y = "kappa", title = "Kappa - inlabru")
  
  data.fixed  <- data.frame(res.inlabru$marginals.fixed)
  
  p.intercept <- ggplot(data.fixed) + 
    geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Estimated", fill="Estimated")) + 
    geom_vline(aes(xintercept = underlying.effects$age.intercept.true, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of phi", y = " ", title = "Intercept - inlabru")
  
  p.random.effects <- (p.intercept | p.alpha ) / (p.beta | p.kappa ) + 
    plot_layout(guides = "collect")
  
  if(save){
    save.figure(p.random.effects, name = "random_effects_inlabru", path = path.to.storage)
  }
  
  data.eta <- data.frame(eta.sim = res.inlabru$summary.linear.predictor$mean[1:length(obs$eta)]) %>%
    mutate(mean  = eta.sim) %>%
    mutate(`0.975quant` = res.inlabru$summary.linear.predictor$`0.975quant`[1:length(obs$eta)]) %>%
    mutate(`0.025quant` = res.inlabru$summary.linear.predictor$`0.025quant`[1:length(obs$eta)]) %>%
    mutate(true.eta = obs$eta) %>%
    mutate(xt = obs$xt, x = obs$x, t = obs$t)
  
  p.eta <- ggplot(data = data.eta) +
    geom_point(aes(x = eta.sim, y = true.eta), color = palette[1]) + 
    labs(x="Estimated eta", y="True value for eta", title = "Eta")
  
  p.eta.2 <- ggplot(data = data.eta) +
    geom_line(aes(x=xt, y = eta.sim, color="Estimated")) +
    geom_line(aes(x=xt, y = true.eta, color="True")) +
    labs(x=" ", y="Eta", title="Eta- inlabru")
  
  p.eta.t <- ggplot(data = data.eta) + 
    geom_line(aes(x = x, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = x, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.eta.x <- ggplot(data = data.eta) + 
    geom_line(aes(x = t, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = t, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each age") + 
    facet_wrap(~x)
  
  p.eta.xt <- (p.eta | p.eta.2) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.xt, name="eta_xt_inlabru", path=path.to.storage)
  }
  
  p.eta.facet <- (p.eta.t | p.eta.x) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.facet, name="eta_facet_inlabru", path=path.to.storage)
  }
  
  plots <- list(p.alpha = p.alpha, p.beta = p.beta, p.kappa = p.kappa,
                p.eta = p.eta,
                p.eta.2 = p.eta.2, p.eta.x = p.eta.x, 
                p.eta.t = p.eta.t, p.intercept = p.intercept,
                p.random.effects = p.random.effects,
                p.eta.all = p.eta.xt,
                p.eta.facet = p.eta.facet)
  
  summaries <- list(data.alpha = data.alpha,
                    data.beta = data.beta,
                    data.kappa = data.kappa,
                    data.eta = data.eta,
                    data.fixed = data.fixed,
                    intercept = res.inlabru$summary.fixed$mean[1],
                    phi = res.inlabru$summary.fixed$mean[2]
  )
  return(list(plots=plots, summaries=summaries))
}

plot.inlabru.vs.underlying.lc.only.kappa.no.intercept <- function(res.inlabru, underlying.effects,
                                                       path.to.storage="", save=FALSE){
  
  obs <- underlying.effects$obs
  nx <- underlying.effects$nx
  nt <- underlying.effects$nt
  
  data.alpha = cbind(res.inlabru$summary.random$alpha,
                     alpha.true = underlying.effects$alpha.true[res.inlabru$summary.random$alpha$ID + 1], 
                     alpha.shifted = underlying.effects$alpha.true[res.inlabru$summary.random$alpha$ID + 1] + underlying.effects$age.intercept.true)
  p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = alpha.shifted, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Alpha - inlabru", x = "x", y='')
  
  data.beta = cbind(res.inlabru$summary.random$beta,
                    beta.true = underlying.effects$beta.true[res.inlabru$summary.random$beta$ID + 1])
  p.beta <- ggplot(data = data.beta, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = beta.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "x", y = "beta", title = "Beta - inlabru")
  
  data.kappa = cbind(res.inlabru$summary.random$kappa,
                     kappa.drifted = underlying.effects$kappa.drifted[res.inlabru$summary.random$kappa$ID + 1])
  p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = kappa.drifted, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "t", y = "kappa", title = "Kappa - inlabru")
  
  p.random.effects <- (p.alpha ) / (p.beta | p.kappa ) + 
    plot_layout(guides = "collect")
  
  if(save){
    save.figure(p.random.effects, name = "random_effects_inlabru", path = path.to.storage)
  }
  
  data.eta <- data.frame(eta.sim = res.inlabru$summary.linear.predictor$mean[1:length(obs$eta)]) %>%
    mutate(mean  = eta.sim) %>%
    mutate(`0.975quant` = res.inlabru$summary.linear.predictor$`0.975quant`[1:length(obs$eta)]) %>%
    mutate(`0.025quant` = res.inlabru$summary.linear.predictor$`0.025quant`[1:length(obs$eta)]) %>%
    mutate(true.eta = obs$eta) %>%
    mutate(xt = obs$xt, x = obs$x, t = obs$t)
  
  p.eta <- ggplot(data = data.eta) +
    geom_point(aes(x = eta.sim, y = true.eta), color = palette[1]) + 
    labs(x="Estimated eta", y="True value for eta", title = "Eta")
  
  p.eta.2 <- ggplot(data = data.eta) +
    geom_line(aes(x=xt, y = eta.sim, color="Estimated")) +
    geom_line(aes(x=xt, y = true.eta, color="True")) +
    labs(x=" ", y="Eta", title="Eta- inlabru")
  
  p.eta.t <- ggplot(data = data.eta) + 
    geom_line(aes(x = x, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = x, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.eta.x <- ggplot(data = data.eta) + 
    geom_line(aes(x = t, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = t, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each age") + 
    facet_wrap(~x)
  
  p.eta.xt <- (p.eta | p.eta.2) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.xt, name="eta_xt_inlabru", path=path.to.storage)
  }
  
  p.eta.facet <- (p.eta.t | p.eta.x) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.facet, name="eta_facet_inlabru", path=path.to.storage)
  }
  
  plots <- list(p.alpha = p.alpha, p.beta = p.beta, p.kappa = p.kappa,
                p.eta = p.eta,
                p.eta.2 = p.eta.2, p.eta.x = p.eta.x, 
                p.eta.t = p.eta.t,
                p.random.effects = p.random.effects,
                p.eta.all = p.eta.xt,
                p.eta.facet = p.eta.facet)
  
  summaries <- list(data.alpha = data.alpha,
                    data.beta = data.beta,
                    data.kappa = data.kappa,
                    data.eta = data.eta,
                    data.fixed = data.fixed,
                    phi = res.inlabru$summary.fixed$mean[2]
  )
  return(list(plots=plots, summaries=summaries))
}


plot.inlabru.vs.underlying.cohort.only.kappa <- function(res.inlabru, underlying.effects,
                                                           path.to.storage="", save=FALSE,
                                                           phi.plus.kappa.func=phi.plus.kappa.v17){
  
  obs <- underlying.effects$obs
  nx <- underlying.effects$nx
  nt <- underlying.effects$nt
  nc <- underlying.effects$nc
  
  data.alpha = cbind(res.inlabru$summary.random$alpha,
                     alpha.true = underlying.effects$alpha.true[res.inlabru$summary.random$alpha$ID + 1])
  p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = alpha.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Alpha - inlabru", x = "x", y='')
  
  data.beta = cbind(res.inlabru$summary.random$beta,
                    beta.true = underlying.effects$beta.true[res.inlabru$summary.random$beta$ID + 1])
  p.beta <- ggplot(data = data.beta, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = beta.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "x", y = "beta", title = "Beta - inlabru")
  
  data.kappa = cbind(res.inlabru$summary.random$kappa,
                     kappa.drifted = underlying.effects$kappa.true[res.inlabru$summary.random$kappa$ID + 1] + underlying.effects$phi.true*res.inlabru$summary.random$kappa$ID)
  p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = kappa.drifted, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "t", y = "kappa", title = "Kappa - inlabru")
  
  data.fixed  <- data.frame(res.inlabru$marginals.fixed)
  
  p.intercept <- ggplot(data.fixed) + 
    geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Estimated", fill="Estimated")) + 
    geom_vline(aes(xintercept = underlying.effects$age.intercept.true, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of phi", y = " ", title = "Intercept - inlabru")
  
  data.gamma = cbind(res.inlabru$summary.random$gamma,
                     gamma.true = underlying.effects$gamma.true[res.inlabru$summary.random$gamma$ID + nx])
  p.gamma <- ggplot(data = data.gamma, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = gamma.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "t", y = "gamma", title = "Gamma - inlabru")
  
  p.random.effects <- (p.intercept | p.alpha ) / (p.beta | p.kappa | p.gamma) + 
    plot_layout(guides = "collect")
  
  if(save){
    save.figure(p.random.effects, name = "random_effects_inlabru", path = path.to.storage)
  }
  
  data.eta <- data.frame(eta.sim = res.inlabru$summary.linear.predictor$mean[1:length(obs$eta)]) %>%
    mutate(true.eta = obs$eta) %>%
    mutate(xt = obs$xt, x = obs$x, t = obs$t)
  
  p.eta <- ggplot(data = data.eta) +
    geom_point(aes(x = eta.sim, y = true.eta), color = palette[1]) + 
    labs(x="Estimated eta", y="True value for eta", title = "Eta")
  
  p.eta.2 <- ggplot(data = data.eta) +
    geom_line(aes(x=xt, y = eta.sim, color="Estimated")) +
    geom_line(aes(x=xt, y = true.eta, color="True")) +
    labs(x=" ", y="Eta", title="Eta- inlabru")
  
  p.eta.t <- ggplot(data = data.eta) + 
    geom_line(aes(x = x, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = x, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.eta.x <- ggplot(data = data.eta) + 
    geom_line(aes(x = t, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = t, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each age") + 
    facet_wrap(~x)
  
  p.eta.xt <- (p.eta | p.eta.2) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.xt, name="eta_xt_inlabru", path=path.to.storage)
  }
  
  p.eta.facet <- (p.eta.t | p.eta.x) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.facet, name="eta_facet_inlabru", path=path.to.storage)
  }
  
  plots <- list(p.alpha = p.alpha, p.beta = p.beta, p.kappa = p.kappa,
                p.eta = p.eta, p.gamma = p.gamma,
                p.eta.2 = p.eta.2, p.eta.x = p.eta.x, 
                p.eta.t = p.eta.t, p.intercept = p.intercept,
                p.random.effects = p.random.effects,
                p.phi.kappa = p.phi.kappa,
                p.eta.all = p.eta.xt,
                p.eta.facet = p.eta.facet)
  
  summaries <- list(data.alpha = data.alpha,
                    data.beta = data.beta,
                    data.kappa = data.kappa,
                    data.eta = data.eta,
                    data.fixed = data.fixed,
                    data.gamma = data.gamma,
                    data.period = results.period$posterior.data,
                    intercept = res.inlabru$summary.fixed$mean[1],
                    phi = res.inlabru$summary.fixed$mean[2]
  )
  return(list(plots=plots, summaries=summaries))
}


plot.inlabru.vs.underlying.cohort.only.kappa.2 <- function(res.inlabru, underlying.effects,
                                                       path.to.storage="", save=FALSE,
                                                       phi.plus.kappa.func=phi.plus.kappa.v17){
  
  obs <- underlying.effects$obs
  nx <- underlying.effects$nx
  nt <- underlying.effects$nt
  nc <- underlying.effects$nc
  
  data.alpha = cbind(res.inlabru$summary.random$alpha,
                     alpha.true = underlying.effects$alpha.true[res.inlabru$summary.random$alpha$ID + 1])
  p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = alpha.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Alpha - inlabru", x = "x", y='')
  
  data.beta = cbind(res.inlabru$summary.random$beta,
                    beta.true = underlying.effects$beta.true[res.inlabru$summary.random$beta$ID + 1])
  p.beta <- ggplot(data = data.beta, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = beta.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "x", y = "beta", title = "Beta - inlabru")
  
  data.kappa = cbind(res.inlabru$summary.random$kappa,
                     kappa.drifted = underlying.effects$kappa.drifted[res.inlabru$summary.random$kappa$ID + 1])
  p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = kappa.drifted, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "t", y = "kappa", title = "Kappa - inlabru")
  
  data.fixed  <- data.frame(res.inlabru$marginals.fixed)
  
  p.intercept <- ggplot(data.fixed) + 
    geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Estimated", fill="Estimated")) + 
    geom_vline(aes(xintercept = underlying.effects$age.intercept.true, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of phi", y = " ", title = "Intercept - inlabru")
  
  data.gamma = cbind(res.inlabru$summary.random$gamma,
                     gamma.true = underlying.effects$gamma.true[res.inlabru$summary.random$gamma$ID + nx])
  p.gamma <- ggplot(data = data.gamma, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = gamma.true, color = "True value", fill = "True value"), size = 0.5) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "t", y = "gamma", title = "Gamma - inlabru")
  
  p.random.effects <- (p.intercept | p.alpha ) / (p.beta | p.kappa | p.gamma) + 
    plot_layout(guides = "collect")
  
  if(save){
    save.figure(p.random.effects, name = "random_effects_inlabru", path = path.to.storage)
  }
  
  data.eta <- data.frame(eta.sim = res.inlabru$summary.linear.predictor$mean[1:length(obs$eta)]) %>%
    mutate(true.eta = obs$eta) %>%
    mutate(xt = obs$xt, x = obs$x, t = obs$t)
  
  p.eta <- ggplot(data = data.eta) +
    geom_point(aes(x = eta.sim, y = true.eta), color = palette[1]) + 
    labs(x="Estimated eta", y="True value for eta", title = "Eta")
  
  p.eta.2 <- ggplot(data = data.eta) +
    geom_line(aes(x=xt, y = eta.sim, color="Estimated")) +
    geom_line(aes(x=xt, y = true.eta, color="True")) +
    labs(x=" ", y="Eta", title="Eta- inlabru")
  
  p.eta.t <- ggplot(data = data.eta) + 
    geom_line(aes(x = x, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = x, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.eta.x <- ggplot(data = data.eta) + 
    geom_line(aes(x = t, y = eta.sim, color = "Estimated")) +
    geom_line(aes(x = t, y = true.eta, color = "True")) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each age") + 
    facet_wrap(~x)
  
  p.eta.xt <- (p.eta | p.eta.2) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.xt, name="eta_xt_inlabru", path=path.to.storage)
  }
  
  p.eta.facet <- (p.eta.t | p.eta.x) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  if(save){
    save.figure(p.eta.facet, name="eta_facet_inlabru", path=path.to.storage)
  }
  
  plots <- list(p.alpha = p.alpha, p.beta = p.beta, p.kappa = p.kappa,
                p.eta = p.eta, p.gamma = p.gamma,
                p.eta.2 = p.eta.2, p.eta.x = p.eta.x, 
                p.eta.t = p.eta.t, p.intercept = p.intercept,
                p.random.effects = p.random.effects,
                p.eta.all = p.eta.xt,
                p.eta.facet = p.eta.facet)
  
  summaries <- list(data.alpha = data.alpha,
                    data.beta = data.beta,
                    data.kappa = data.kappa,
                    data.eta = data.eta,
                    data.fixed = data.fixed,
                    data.gamma = data.gamma,
                    intercept = res.inlabru$summary.fixed$mean[1]
  )
  return(list(plots=plots, summaries=summaries))
}

