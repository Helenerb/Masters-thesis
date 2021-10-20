library(ggplot2)
library(patchwork)

#   ----   Source relevant functions
source("../Functions/plotters.R")
source("../Misc/palette.R")

phi.plus.kappa.v7 <- function(){
  t = 0:(100-1)
  res = kappa + phi*t
  return(res)
}

phi.plus.kappa.v17 <- function(){
  t = 0:(20-1)
  res = kappa + phi*t
  return(res)
}

plot.posterior.period.effects <- function(res.inlabru, underlying.effects, phi.plus.kappa.func){
  
  nt = underlying.effects$nt
  
  samps = inla.posterior.sample(res.inlabru, n = 1000)
  
  posterior.phi.kappa <- inla.posterior.sample.eval(fun = phi.plus.kappa.func, samples=samps)
  
  posterior.phi.kappa.df <- data.frame(t = 1:nt,
                                       mean = apply(posterior.phi.kappa, 1, mean),
                                       q1 = apply(posterior.phi.kappa, 1, quantile, 0.025),
                                       q2 = apply(posterior.phi.kappa, 1, quantile, 0.975)) %>%
    mutate(kappa = underlying.effects$kappa.true[t]) %>%
    mutate(phi.t = underlying.effects$phi.true*(t-1)) %>%
    mutate(kappa.phi = kappa + phi.t)
  
  gg.posterior <- ggplot(data = posterior.phi.kappa.df) +
    geom_ribbon(aes(x = t, ymin = q1, ymax = q2, fill = "Estimated"), alpha = 0.5) + 
    geom_point(aes(x = t, y = mean, color = "Estimated", fill = "Estimated"), size=0.5) +
    geom_point(aes(x = t, y = kappa.phi, color = "True", fill = "True"), size=0.5) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) + 
    labs(title = "Phi*t + kappa", x = "t", y = "")
  
  return(list(posterior.plot = gg.posterior, posterior.data = posterior.phi.kappa.df))
}

plot.inlabru.vs.underlying.v1 <- function(res.inlabru, underlying.effects){
   
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
   
   data.gamma = cbind(
     res.inlabru$summary.random$gamma,
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
   
   p.phi <- ggplot(data.frame(res.inlabru$marginals.fixed)) + 
     geom_area(aes(x = phi.x, y = phi.y, fill = "Estimated"), alpha = 0.4) + 
     geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[2], color = "Inlabru", fill="Inlabru")) + 
     geom_vline(aes(xintercept = underlying.effects$phi.true, color="True", fill="True")) +
     scale_color_manual(name = " ", values = palette) + 
     scale_fill_manual(name = " ", values = palette) +
     labs(x = "Value of phi", y = " ", title = "Phi - inlabru")
   
   p.intercept <- ggplot(data.frame(res.inlabru$marginals.fixed)) + 
     geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
     geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Inlabru", fill="Inlabru")) + 
     geom_vline(data = underlying.effects, aes(xintercepet = age.intercept.true, color = "True", fill="True" ))
     scale_color_manual(name = " ", values = palette) + 
     scale_fill_manual(name = " ", values = palette) +
     labs(x = "Value of intercept", y = " ", title = "Intercept - inlabru")
   
   data.eta <- data.frame(eta.sim = res.inlabru$summary.linear.predictor$mean[1:length(obs$eta)]) %>%
     mutate(true.eta = obs$eta) %>%
     mutate(xt = obs$xt, x = obs$x, t = obs$t)
   
   p.eta <- ggplot(data = data.eta) +
     geom_point(aes(x = eta.sim, y = true.eta), color = palette[1]) + 
     labs(x="Estimated eta", y="True value for eta", title = "Eta")
   
   p.eta.2 <- ggplot(data = data.eta) +
     geom_point(aes(x=xt, y = eta.sim, color="Inlabru"), size = 0.5) +
     geom_point(aes(x=xt, y = true.eta, color="True"), size = 0.5) +
     labs(x=" ", y="Eta", title="Eta- inlabru")
   
   p.eta.t <- ggplot(data = data.eta) + 
     geom_line(aes(x = x, y = eta.sim, color = "Inlabru")) +
     geom_line(aes(x = x, y = true.eta, color = "True")) +
     labs(x = " ", y = " ", title = "Eta - inlabru, for each year") + 
     facet_wrap(~t)
   
   p.eta.x <- ggplot(data = data.eta) + 
     geom_line(aes(x = t, y = eta.sim, color = "Inlabru")) +
     geom_line(aes(x = t, y = true.eta, color = "True")) +
     labs(x = " ", y = " ", title = "Eta - inlabru, for each age") + 
     facet_wrap(~x)
   
   
   plots <- list(p.alpha = p.alpha, p.beta = p.beta, p.kappa = p.kappa,
                 p.phi = p.phi, p.gamma = p.gamma, p.eta = p.eta,
                 p.eta.2 = p.eta.2, p.eta.x = p.eta.x, p.eta.t = p.eta.t,
                 p.intercept = p.intercept)
   return(plots)
}

plot.inlabru.vs.underlying.v5 <- function(res.inlabru, underlying.effects){
  
  obs <- underlying.effects$obs
  
  data.alpha = cbind(res.inlabru$summary.random$alpha,
                     alpha.true = underlying.effects$alpha.true[res.inlabru$summary.random$alpha$ID + 1])
  p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = alpha.true, color = "True value", fill = "True value")) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Alpha - inlabru", x = "x", y='')
  
  data.beta = cbind(res.inlabru$summary.random$beta,
                    beta.true = underlying.effects$beta.true[res.inlabru$summary.random$beta$ID + 1])
  p.beta <- ggplot(data = data.beta, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = beta.true, color = "True value", fill = "True value")) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "x", y = "beta", title = "Beta - inlabru")
  
  data.kappa = cbind(res.inlabru$summary.random$kappa,
                     kappa.true = underlying.effects$kappa.true[res.inlabru$summary.random$kappa$ID + 1])
  p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = kappa.true, color = "True value", fill = "True value")) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "t", y = "kappa", title = "Kappa - inlabru")
  
  p.phi <- ggplot(data.frame(res.inlabru$marginals.fixed)) + 
    geom_area(aes(x = phi.x, y = phi.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[2], color = "Estimated", fill="Estimated")) + 
    geom_vline(aes(xintercept = underlying.effects$phi.true, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of phi", y = " ", title = "Phi - inlabru")
  
  p.intercept <- ggplot(data.frame(res.inlabru$marginals.fixed)) + 
    geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Estimated", fill="Estimated")) + 
    geom_vline(aes(xintercept = underlying.effects$age.intercept.true, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of phi", y = " ", title = "Intercept - inlabru")
  
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
  
  
  plots <- list(p.alpha = p.alpha, p.beta = p.beta, p.kappa = p.kappa,
                p.phi = p.phi, p.eta = p.eta,
                p.eta.2 = p.eta.2, p.eta.x = p.eta.x, 
                p.eta.t = p.eta.t, p.intercept = p.intercept)
  return(plots)
}

plot.inlabru.vs.underlying.v3 <- function(res.inlabru, underlying.effects){
  
  obs <- underlying.effects$obs
  
  data.alpha = cbind(res.inlabru$summary.random$alpha,
                     alpha.true = underlying.effects$alpha.true[res.inlabru$summary.random$alpha$ID])
  p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = alpha.true, color = "True value", fill = "True value")) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
    geom_point(aes(y = mean + res.inlabru$summary.fixed$mean[1], fill = "Estimated + intercept", color = "Estimated + intercept")) +
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Alpha - inlabru", x = "x", y='')
  
  data.beta = cbind(res.inlabru$summary.random$beta,
                    beta.true = underlying.effects$beta.true[res.inlabru$summary.random$beta$ID])
  p.beta <- ggplot(data = data.beta, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = beta.true, color = "True value", fill = "True value")) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "x", y = "beta", title = "Beta - inlabru")
  
  data.kappa = cbind(res.inlabru$summary.random$kappa,
                     kappa.true = underlying.effects$kappa.true[res.inlabru$summary.random$kappa$ID])
  p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = kappa.true, color = "True value", fill = "True value")) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(x = "t", y = "kappa", title = "Kappa - inlabru")
  
  p.phi <- ggplot(data.frame(res.inlabru$marginals.fixed)) + 
    geom_area(aes(x = phi.x, y = phi.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[2], color = "Estimated", fill="Estimated")) + 
    geom_vline(aes(xintercept = underlying.effects$phi.true, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of phi", y = " ", title = "Phi - inlabru")
  
  p.intercept <- ggplot(data.frame(res.inlabru$marginals.fixed)) + 
    geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Estimated", fill="Estimated")) + 
    #geom_vline(aes(xintercept = underlying.effects$.true, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of phi", y = " ", title = "Intercept - inlabru")
  
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
    
  
  plots <- list(p.alpha = p.alpha, p.beta = p.beta, p.kappa = p.kappa,
                p.phi = p.phi, p.eta = p.eta,
                p.eta.2 = p.eta.2, p.eta.x = p.eta.x, 
                p.eta.t = p.eta.t, p.intercept = p.intercept)
  return(plots)
}

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

plot.inlabru.vs.underlying.lc.only.kappa <- function(res.inlabru, underlying.effects,
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
                     kappa.true = underlying.effects$kappa.true[res.inlabru$summary.random$kappa$ID + 1],
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
  
  p.random.effects <- (p.intercept | p.alpha ) / (p.beta | p.kappa ) + 
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



# Deprecate?
plot_inlabru_results <- function(res.inlabru, underlying.effects, produce_plot_func){
  obs <- underlying.effects$obs
  
  #   ----    Plot remaining random effects    ----
  
  source("plot_inlabru_vs_underlying.R")  # ikke source, men flytt den hit"
  
  # plotting results from run with cohort effects:
  plots <- plot.inlabru.vs.underlying.v5(res.inlabru, underlying.effects)
  p.alpha <- plots$p.alpha
  pl.beta <- plots$p.beta
  p.phi <- plots$p.phi
  p.intercept <- plots$p.intercept
  p.kappa <- plots$p.kappa
  p.eta <- plots$p.eta
  p.eta.2 <- plots$p.eta.2
  p.eta.t <- plots$p.eta.t
  p.eta.t <- plots$p.eta.x
  
  p.random.effects <- (p.intercept | p.alpha)/(p.beta | p.posterior) + 
    plot_layout(guides = "collect")
  
  p.phi.and.kappa <- (p.phi | p.kappa) + 
    plot_layout(guides = "collect")
  
  p.eta <- (p.eta | p.eta.2)/(p.eta.t | p.eta.x) + 
    plot_layout(guides = "collect")
  
}