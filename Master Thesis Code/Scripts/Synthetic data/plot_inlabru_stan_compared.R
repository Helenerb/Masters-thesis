# script contatining functions for comparing inlabru and stan results

library(ggplot2)
library(patchwork)

source("../Functions/plotters.R")
source("../Misc/palette.R")

plot.inlabru.stan.compared.lc <- function(stan.summaries, inlabru.summaries, underlying.effects){
  obs <- underlying.effects$obs
  
  #  ----   intercept   ----
  print(stan.summaries$summary_fixed)
  p.intercept <- ggplot() + 
    geom_area(data=inlabru.summaries$data.fixed, aes(x = Int.x, y = Int.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.4, size = 0.5) + 
    geom_vline(aes(xintercept = inlabru.summaries$intercept, color = "Inlabru", fill="Inlabru")) + 
    geom_vline(data=stan.summaries$summary_fixed, aes(xintercept = mean[2], fill = "Stan", color = "Stan")) + 
    geom_vline(data=stan.summaries$summary_fixed, aes(xintercept = `2.5%`[2], fill = "Stan", color = "Stan"), alpha = 0.5) +
    geom_vline(data=stan.summaries$summary_fixed, aes(xintercept = `97.5%`[2], fill = "Stan", color = "Stan"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$age.intercept.true, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of intercept", y = " ", title = "Intercept")
  
  # ---   alpha   ----
  p.alpha <- ggplot() + 
    geom_ribbon(data=inlabru.summaries$data.alpha, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
    geom_point(data=inlabru.summaries$data.alpha, aes(x = ID, y = mean, color = "Inlabru", fill = "Inlabru"), size = 0.5) + 
    geom_point(data=stan.summaries$summary_alpha, aes(x=index - 1, y=mean, fill="Stan", color="Stan"), size=0.5) + 
    geom_line(data=stan.summaries$summary_alpha, aes(x=index - 1, y=`2.5%`, fill = "Stan", color="Stan"), alpha=0.5) + 
    geom_line(data=stan.summaries$summary_alpha, aes(x=index - 1, y=`97.5%`, fill="Stan", color="Stan"), alpha=0.5) +
  
    geom_point(data=inlabru.summaries$data.alpha, 
               aes(x = ID, y = underlying.effects$alpha.true, color = "True", fill = "True"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Alpha", x = "x", y='')
  
  # ---   beta   ----
  p.beta <- ggplot() + 
    geom_ribbon(data=inlabru.summaries$data.beta, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
    geom_point(data=inlabru.summaries$data.beta, aes(x = ID, y = mean, color = "Inlabru", fill = "Inlabru"), size = 0.5) + 
    
    geom_point(data=stan.summaries$summary_beta, aes(x=index - 1, y=mean, fill="Stan", color="Stan"), size=0.5) + 
    geom_line(data=stan.summaries$summary_beta, aes(x=index - 1, y=`2.5%`, fill="Stan", color="Stan"), alpha=0.5) + 
    geom_line(data=stan.summaries$summary_beta, aes(x=index - 1, y=`97.5%`, fill="Stan", color="Stan"), alpha=0.5) +
    
    geom_point(data=inlabru.summaries$data.beta, 
               aes(x = ID, y = underlying.effects$beta.true, color = "True", fill = "True"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Beta", x = "x", y='')
  
  #   ----   kappa   ----  
  p.kappa <- ggplot() + 
    geom_ribbon(data=inlabru.summaries$data.kappa, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
    geom_point(data=inlabru.summaries$data.kappa, aes(x = ID, y = mean, color = "Inlabru", fill = "Inlabru"), size = 0.5) + 
    
    geom_point(data=stan.summaries$summary_kappa, aes(x=index - 1, y=mean, fill="Stan",color="Stan"), size=0.5) + 
    geom_line(data=stan.summaries$summary_kappa, aes(x=index - 1, y=`2.5%`, fill="Stan", color="Stan"), alpha=0.5) + 
    geom_line(data=stan.summaries$summary_kappa, aes(x=index - 1, y=`97.5%`, fill="Stan", color="Stan"), alpha=0.5) +
    
    geom_point(data=inlabru.summaries$data.kappa, 
               aes(x = ID, y = underlying.effects$kappa.true, color = "True", fill = "True"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Kappa", x = "t", y='')
  
  #   ----   phi   ----
  p.phi <- ggplot() + 
    geom_area(data=inlabru.summaries$data.fixed, aes(x = phi.x, y = phi.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.4, size=0.5) + 
    geom_vline(aes(xintercept = inlabru.summaries$phi, color = "Inlabru", fill="Inlabru")) + 
    geom_vline(data=stan.summaries$summary_fixed, aes(xintercept = mean[1], fill="Stan", color = "Stan")) + 
    geom_vline(data=stan.summaries$summary_fixed, aes(xintercept = `2.5%`[1], fill="Stan", color = "Stan"), alpha = 0.5) +
    geom_vline(data=stan.summaries$summary_fixed, aes(xintercept = `97.5%`[1], fill="Stan", color = "Stan"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$phi.true, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of phi", y = " ", title = "Phi")
  
  #   ----   period   ----
  p.period <- ggplot() +
    geom_ribbon(data=inlabru.summaries$data.period, aes(x = t, ymin = q1, ymax = q2, fill = "Inlabru"), alpha = 0.5) + 
    geom_point(data=inlabru.summaries$data.period, aes(x = t, y = mean, color = "Inlabru", fill = "Inlabru"), size=0.5) +
    
    geom_point(data=stan.summaries$summary_kappa, aes(x=index, y=mean, fill="Stan", color="Stan"), size=0.5) + 
    geom_line(data=stan.summaries$summary_kappa, aes(x=index, y=`2.5%`, fill="Stan", color="Stan"), alpha=0.5) + 
    geom_line(data=stan.summaries$summary_kappa, aes(x=index, y=`97.5%`, fill="Stan", color="Stan"), alpha=0.5) +
    
    geom_point(data=inlabru.summaries$data.period, aes(x = t, y = kappa.phi, color = "True", fill = "True"), size=0.5) +
    
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) + 
    labs(title = "Period effect", x = "t", y = "")
  
  #   ----   eta   ----
  p.eta <- ggplot() +
    geom_point(data=inlabru.summaries$data.eta, aes(x = eta.sim, y = true.eta, color = "Inlabru")) + 
    geom_point(data=stan.summaries$summary_eta, aes(x = mean, y = true_eta, color = "Stan")) + 
    scale_color_manual(name = " ", values = palette) + 
    labs(x="Estimated eta", y="True value for eta", title = "Eta")
  
  p.eta.2 <- ggplot() +
    geom_line(data = inlabru.summaries$data.eta, aes(x=xt, y = eta.sim, color="Inlabru")) +
    geom_line(data = inlabru.summaries$data.eta, aes(x=xt, y = true.eta, color="True")) +
    geom_line(data = stan.summaries$summary_eta, aes(x=xt, y = mean, color="True")) +
    scale_color_manual(name = "", values = palette ) +
    labs(x=" ", y="Eta", title="Eta")
  
  p.eta.t <- ggplot() + 
    geom_line(data = inlabru.summaries$data.eta, aes(x = x, y = eta.sim, color = "Inlabru")) +
    geom_line(data = inlabru.summaries$data.eta, aes(x = x, y = true.eta, color = "True")) +
    geom_line(data=stan.summaries$summary_eta, aes(x = x, y = mean, color = "Stan")) +
    scale_color_manual(name = "", values = palette ) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.eta.x <- ggplot() + 
    geom_line(data = inlabru.summaries$data.eta, aes(x = t, y = eta.sim, color = "Inlabru")) +
    geom_line(data = inlabru.summaries$data.eta, aes(x = t, y = true.eta, color = "True")) +
    geom_line(data=stan.summaries$summary_eta, aes(x = t, y = mean, color = "Stan")) +
    scale_color_manual(name = "", values = palette ) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each age") + 
    facet_wrap(~x)
  
  plots <- list(p.intercept = p.intercept, 
                p.alpha = p.alpha, 
                p.beta = p.beta,
                p.kappa = p.kappa,
                p.phi = p.phi,
                p.period = p.period,
                p.eta = p.eta,
                p.eta.2 = p.eta.2,
                p.eta.t = p.eta.t,
                p.eta.x = p.eta.x)
  
  return(plots)
}

plot.inlabru.stan.compared.cohort <- function(stan.summaries, inlabru.summaries, underlying.effects){
  obs <- underlying.effects$obs
  
  #  ----   intercept   ----
  p.intercept <- ggplot() + 
    geom_area(data=inlabru.summaries$data.fixed, aes(x = Int.x, y = Int.y, fill = "Inlabru", color = "Inlabru"), alpha = 0.4, size = 0.5) + 
    geom_vline(aes(xintercept = inlabru.summaries$intercept, color = "Inlabru", fill="Inlabru")) + 
    geom_vline(data=stan.summaries$summary_fixed, aes(xintercept = mean[2], fill="Stan", color = "Stan")) + 
    geom_vline(data=stan.summaries$summary_fixed, aes(xintercept = `2.5%`[2], fill="Stan", color = "Stan"), alpha = 0.5) +
    geom_vline(data=stan.summaries$summary_fixed, aes(xintercept = `97.5%`[2], fill="Stan", color = "Stan"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$age.intercept.true, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of intercept", y = " ", title = "Intercept")
  
  # ---   alpha   ----
  p.alpha <- ggplot() + 
    geom_ribbon(data=inlabru.summaries$data.alpha, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
    geom_point(data=inlabru.summaries$data.alpha, aes(x = ID, y = mean, color = "Inlabru", fill = "Inlabru"), size = 0.5) + 
    
    geom_point(data=stan.summaries$summary_alpha, aes(x=index - 1, y=mean, fill="Stan", color="Stan"), size=0.5) + 
    geom_line(data=stan.summaries$summary_alpha, aes(x=index - 1, y=`2.5%`, fill="Stan", color="Stan"), alpha=0.5) + 
    geom_line(data=stan.summaries$summary_alpha, aes(x=index - 1, y=`97.5%`, fill="Stan", color="Stan"), alpha=0.5) +
    
    geom_point(data=inlabru.summaries$data.alpha, 
               aes(x = ID, y = underlying.effects$alpha.true, color = "True", fill = "True"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Alpha", x = "x", y='')
  
  # ---   beta   ----
  p.beta <- ggplot() + 
    geom_ribbon(data=inlabru.summaries$data.beta, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
    geom_point(data=inlabru.summaries$data.beta, aes(x = ID, y = mean, color = "Inlabru", fill = "Inlabru"), size = 0.5) + 
    
    geom_point(data=stan.summaries$summary_beta, aes(x=index - 1, y=mean, fill="Stan", color="Stan"), size=0.5) + 
    geom_line(data=stan.summaries$summary_beta, aes(x=index - 1, y=`2.5%`, fill="Stan", color="Stan"), alpha=0.5) + 
    geom_line(data=stan.summaries$summary_beta, aes(x=index - 1, y=`97.5%`, fill="Stan", color="Stan"), alpha=0.5) +
    
    geom_point(data=inlabru.summaries$data.beta, 
               aes(x = ID, y = underlying.effects$beta.true, color = "True", fill = "True"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Beta", x = "x", y='')
  
  #   ----   kappa   ----  
  p.kappa <- ggplot() + 
    geom_ribbon(data=inlabru.summaries$data.kappa, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
    geom_point(data=inlabru.summaries$data.kappa, aes(x = ID, y = mean, color = "Inlabru", fill = "Inlabru"), size = 0.5) + 
    
    geom_point(data=stan.summaries$summary_kappa, aes(x=index - 1, y=mean, fill="Stan", color="Stan"), size=0.5) + 
    geom_line(data=stan.summaries$summary_kappa, aes(x=index - 1, y=`2.5%`, fill="Stan", color="Stan"), alpha=0.5) + 
    geom_line(data=stan.summaries$summary_kappa, aes(x=index - 1, y=`97.5%`, fill="Stan", color="Stan"), alpha=0.5) +
    
    geom_point(data=inlabru.summaries$data.kappa, 
               aes(x = ID, y = underlying.effects$kappa.true, color = "True", fill = "True"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Kappa", x = "t", y='')
  
  #   ----   phi   ----
  p.phi <- ggplot() + 
    geom_area(data=inlabru.summaries$data.fixed, aes(x = phi.x, y = phi.y, fill = "Inlabru"), alpha = 0.4, size = 0.5) + 
    geom_vline(aes(xintercept = inlabru.summaries$phi, color = "Inlabru", fill="Inlabru")) + 
    geom_vline(data=stan.summaries$summary_fixed, aes(xintercept = mean[1], fill = "Stan", color = "Stan")) + 
    geom_vline(data=stan.summaries$summary_fixed, aes(xintercept = `2.5%`[1], fill = "Stan", color = "Stan"), alpha = 0.5) +
    geom_vline(data=stan.summaries$summary_fixed, aes(xintercept = `97.5%`[1], fill = "Stan", color = "Stan"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$phi.true, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of phi", y = " ", title = "Phi")
  
  #   ----   period   ----
  p.period <- ggplot() +
    geom_ribbon(data=inlabru.summaries$data.period, aes(x = t, ymin = q1, ymax = q2, fill = "Estimated"), alpha = 0.5) + 
    geom_point(data=inlabru.summaries$data.period, aes(x = t, y = mean, color = "Estimated", fill = "Estimated"), size=0.5) +
    
    geom_point(data=stan.summaries$summary_kappa, aes(x=index, y=mean, fill="Stan", color="Stan")) + 
    geom_line(data=stan.summaries$summary_kappa, aes(x=index, y=`2.5%`, fill="Stan", color="Stan"), alpha=0.5) + 
    geom_line(data=stan.summaries$summary_kappa, aes(x=index, y=`97.5%`, fill="Stan", color="Stan"), alpha=0.5) +
    
    geom_point(data=inlabru.summaries$data.period, aes(x = t, y = kappa.phi, color = "True", fill = "True"), size=0.5) +
    
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) + 
    labs(title = "Period effect", x = "t", y = "")
  
  #   ----   gamma   ----  
  p.gamma <- ggplot() + 
    geom_ribbon(data=inlabru.summaries$data.gamma, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
    geom_point(data=inlabru.summaries$data.gamma, aes(x = ID, y = mean, color = "Inlabru", fill = "Inlabru"), size = 0.5) + 
    
    geom_point(data=stan.summaries$summary_gamma,
               aes(x=index - underlying.effects$nx, y=mean, fill="Stan", color="Stan"), size=0.5) + 
    geom_line(data=stan.summaries$summary_gamma,
              aes(x=index - underlying.effects$nx, y=`2.5%`, fill="Stan", color="Stan"), alpha=0.5) + 
    geom_line(data=stan.summaries$summary_gamma,
              aes(x=index - underlying.effects$nx, y=`97.5%`, fill="Stan", color="Stan"), alpha=0.5) +
    
    geom_point(data=inlabru.summaries$data.gamma, 
               aes(x = ID, y = underlying.effects$gamma.true, color = "True", fill = "True"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Gamma", x = "t", y='')
  
  #   ----   eta   ----
  p.eta <- ggplot() +
    geom_point(data=inlabru.summaries$data.eta, aes(x = eta.sim, y = true.eta, color = "Inlabru")) + 
    geom_point(data=stan.summaries$summary_eta, aes(x = mean, y = true_eta, color = "Stan")) + 
    scale_color_manual(name = " ", values = palette) + 
    labs(x="Estimated eta", y="True value for eta", title = "Eta")
  
  p.eta.2 <- ggplot() +
    geom_line(data = inlabru.summaries$data.eta, aes(x=xt, y = eta.sim, color="Inlabru")) +
    geom_line(data = inlabru.summaries$data.eta, aes(x=xt, y = true.eta, color="True")) +
    geom_line(data = stan.summaries$summary_eta, aes(x=xt, y = mean, color="Stan")) +
    scale_color_manual(name = "", values = palette ) +
    labs(x=" ", y="Eta", title="Eta")
  
  p.eta.t <- ggplot() + 
    geom_line(data = inlabru.summaries$data.eta, aes(x = x, y = eta.sim, color = "Inlabru")) +
    geom_line(data = inlabru.summaries$data.eta, aes(x = x, y = true.eta, color = "True")) +
    geom_line(data=stan.summaries$summary_eta, aes(x = x, y = mean, color = "Stan")) +
    scale_color_manual(name = "", values = palette ) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.eta.x <- ggplot() + 
    geom_line(data = inlabru.summaries$data.eta, aes(x = t, y = eta.sim, color = "Inlabru")) +
    geom_line(data = inlabru.summaries$data.eta, aes(x = t, y = true.eta, color = "True")) +
    geom_line(data=stan.summaries$summary_eta, aes(x = t, y = mean, color = "Stan")) +
    scale_color_manual(name = "", values = palette ) +
    labs(x = " ", y = " ", title = "Eta - inlabru, for each age") + 
    facet_wrap(~x)
  
  plots <- list(p.intercept = p.intercept, 
                p.alpha = p.alpha, 
                p.beta = p.beta,
                p.kappa = p.kappa,
                p.phi = p.phi,
                p.period = p.period,
                p.gamma = p.gamma,
                p.eta = p.eta,
                p.eta.2 = p.eta.2,
                p.eta.t = p.eta.t,
                p.eta.x = p.eta.x)
}

save.compared.drifted.lc <- function(plots, path.to.storage){
  p.random.effects <- (plots$p.intercept | plots$p.alpha)/(plots$p.beta | plots$p.period) + 
    plot_layout(guides="collect")
  save.figure(p.random.effects, name = "random_effects_comparison", path = path.to.storage)
  
  p.phi.kappa <- (plots$p.phi | plots$p.kappa) + plot_layout(guides="collect")
  save.figure(p.phi.kappa, name = "phi_kappa_comparison", path = path.to.storage)
  
  p.eta.xt <- (plots$p.eta | plots$p.eta.2) + plot_layout(guides = "collect")
  save.figure(p.eta.xt, name = "eta_xt_comparison", path = path.to.storage)
  
  p.eta.facet <- (plots$p.eta.x | plots$p.eta.t) + plot_layout(guides = "collect")
  save.figure(p.eta.facet, name = "eta_facet_comparison", path = path.to.storage)
}

save.compared.undrifted.lc <- function(plots, path.to.storage){
  p.random.effects <- (plots$p.intercept | plots$p.alpha | plots$p.beta)/(plots$p.phi | plots$p.kappa) + 
    plot_layout(guides="collect")
  save.figure(p.random.effects, name = "random_effects_comparison", path = path.to.storage)
  
  p.eta.xt <- (plots$p.eta | plots$p.eta.2) + plot_layout(guides = "collect")
  save.figure(p.eta.xt, name = "eta_xt_comparison", path = path.to.storage)
  
  p.eta.facet <- (plots$p.eta.x | plots$p.eta.t) + plot_layout(guides = "collect")
  save.figure(p.eta.facet, name = "eta_facet_comparison", path = path.to.storage)
}

save.compared.drifted.cohort <- function(plots, path.to.storage){
  p.random.effects <- (plots$p.intercept | plots$p.alpha | plots$p.beta)/(plots$p.period | plots$p.gamma) + 
    plot_layout(guides="collect")
  save.figure(p.random.effects, name = "random_effects_comparison", path = path.to.storage)
  
  p.phi.kappa <- (plots$p.phi | plots$p.kappa) + plot_layout(guides="collect")
  save.figure(p.phi.kappa, name = "phi_kappa_comparison", path = path.to.storage)
  
  p.eta.xt <- (plots$p.eta | plots$p.eta.2) + plot_layout(guides = "collect")
  save.figure(p.eta.xt, name = "eta_xt_comparison", path = path.to.storage)
  
  p.eta.facet <- (plots$p.eta.x | plots$p.eta.t) + plot_layout(guides = "collect")
  save.figure(p.eta.facet, name = "eta_facet_comparison", path = path.to.storage)
}

save.compared.undrifted.cohort <- function(plots, path.to.storage){
  p.random.effects <- (plots$p.intercept | plots$p.alpha | plots$p.beta)/( plots$p.phi | plots$p.kappa | plots$p.gamma) + 
    plot_layout(guides="collect")
  save.figure(p.random.effects, name = "random_effects_comparison", path = path.to.storage)
  
  p.eta.xt <- (plots$p.eta | plots$p.eta.2) + plot_layout(guides = "collect")
  save.figure(p.eta.xt, name = "eta_xt_comparison", path = path.to.storage)
  
  p.eta.facet <- (plots$p.eta.x | plots$p.eta.t) + plot_layout(guides = "collect")
  save.figure(p.eta.facet, name = "eta_facet_comparison", path = path.to.storage)
}

produce.compared.plots <- function(
  stan.summaries, inlabru.summaries, underlying.effects, plot.func, save.func,
  path.to.storage){
  
  plots <- plot.func(stan.summaries, inlabru.summaries, underlying.effects)
  save.func(plots, path.to.storage)
  return(plots)
}