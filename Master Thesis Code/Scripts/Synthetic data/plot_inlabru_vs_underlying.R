library(ggplot2)

plot.inlabru.vs.underlying.v1 <- function(res.inlabru, underlying.effects){
   
   palette.basis <- c('#70A4D4', '#ECC64B', '#93AD80', '#da9124', '#696B8D',
                      '#3290c1',
                      '#5d8060', '#D7B36A', '#826133', '#A85150')
   
   obs <- underlying.effects$obs
   
   data.alpha = cbind(res.inlabru$summary.random$alpha,
                      alpha.true = underlying.effects$alpha.true[res.inlabru$summary.random$alpha$ID + 1])
   p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
     geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
     geom_point(aes(y = alpha.true, color = "True value", fill = "True value")) + 
     geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
     scale_color_manual(name = "",
                        values = palette.basis ) +
     scale_fill_manual(name = "",
                       values = palette.basis ) +
     labs(title="Alpha - inlabru", x = "x", y='')
   
   data.beta = cbind(res.inlabru$summary.random$beta,
                     beta.true = underlying.effects$beta.true[res.inlabru$summary.random$beta$ID + 1])
   p.beta <- ggplot(data = data.beta, aes(x = ID)) + 
     geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
     geom_point(aes(y = beta.true, color = "True value", fill = "True value")) + 
     geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
     scale_color_manual(name = "",
                        values = palette.basis ) +
     scale_fill_manual(name = "",
                       values = palette.basis ) +
     labs(x = "x", y = "beta", title = "Beta - inlabru")
   
   data.kappa = cbind(res.inlabru$summary.random$kappa,
                      kappa.true = underlying.effects$kappa.true[res.inlabru$summary.random$kappa$ID + 1])
   p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
     geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
     geom_point(aes(y = kappa.true, color = "True value", fill = "True value")) + 
     geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
     scale_color_manual(name = "",
                        values = palette.basis ) +
     scale_fill_manual(name = "",
                       values = palette.basis ) +
     labs(x = "t", y = "kappa", title = "Kappa - inlabru")
   
   data.gamma = cbind(
     res.inlabru$summary.random$gamma,
     gamma.true = underlying.effects$gamma.true[res.inlabru$summary.random$gamma$ID + 79 + 1])
   p.gamma <- ggplot(data = data.gamma, aes(x = ID)) + 
     geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
     geom_point(aes(y = gamma.true, color = "True value", fill = "True value")) + 
     geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
     scale_color_manual(name = "",
                        values = palette.basis ) +
     scale_fill_manual(name = "",
                       values = palette.basis ) +
     labs(x = "t", y = "gamma", title = "Gamma - inlabru")
   
   p.phi <- ggplot(data.frame(res.inlabru$marginals.fixed)) + 
     geom_area(aes(x = phi.x, y = phi.y, fill = "Estimated"), alpha = 0.4) + 
     geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[2], color = "Estimated", fill="Estimated")) + 
     geom_vline(aes(xintercept = underlying.effects$phi.true, color="True", fill="True")) +
     scale_color_manual(name = " ", values = palette.basis) + 
     scale_fill_manual(name = " ", values = palette.basis) +
     labs(x = "Value of phi", y = " ", title = "Phi - inlabru")
   
   p.intercept <- ggplot(data.frame(res.inlabru$marginals.fixed)) + 
     geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
     geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Estimated", fill="Estimated")) + 
     scale_color_manual(name = " ", values = palette.basis) + 
     scale_fill_manual(name = " ", values = palette.basis) +
     labs(x = "Value of phi", y = " ", title = "Intercept - inlabru")
   
   data.eta <- data.frame(eta.sim = res.inlabru$summary.linear.predictor$mean[1:length(obs$eta)]) %>%
     mutate(true.eta = obs$eta) %>%
     mutate(xt = obs$xt, x = obs$x, t = obs$t)
   
   p.eta <- ggplot(data = data.eta) +
     geom_point(aes(x = eta.sim, y = true.eta), color = palette.basis[1]) + 
     labs(x="Estimated eta", y="True value for eta", title = "Eta")
   
   p.eta.2 <- ggplot(data = data.eta) +
     geom_point(aes(x=xt, y = eta.sim, color="Estimated")) +
     geom_point(aes(x=xt, y = true.eta, color="True")) +
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
                 p.phi = p.phi, p.gamma = p.gamma, p.eta = p.eta,
                 p.eta.2 = p.eta.2, p.eta.x = p.eta.x, p.eta.t = p.eta.t,
                 p.intercept = p.intercept)
   return(plots)
}

plot.inlabru.vs.underlying.v5 <- function(res.inlabru, underlying.effects){
  
  palette.basis <- c('#70A4D4', '#ECC64B', '#93AD80', '#da9124', '#696B8D',
                     '#3290c1',
                     '#5d8060', '#D7B36A', '#826133', '#A85150')
  
  obs <- underlying.effects$obs
  
  data.alpha = cbind(res.inlabru$summary.random$alpha,
                     alpha.true = underlying.effects$alpha.true[res.inlabru$summary.random$alpha$ID + 1])
  p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = alpha.true, color = "True value", fill = "True value")) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
    geom_point(aes(y = mean + res.inlabru$summary.fixed$mean[1], fill = "Estimated + intercept", color = "Estimated + intercept")) +
    scale_color_manual(name = "",
                       values = palette.basis ) +
    scale_fill_manual(name = "",
                      values = palette.basis ) +
    labs(title="Alpha - inlabru", x = "x", y='')
  
  data.beta = cbind(res.inlabru$summary.random$beta,
                    beta.true = underlying.effects$beta.true[res.inlabru$summary.random$beta$ID + 1])
  p.beta <- ggplot(data = data.beta, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = beta.true, color = "True value", fill = "True value")) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
    scale_color_manual(name = "",
                       values = palette.basis ) +
    scale_fill_manual(name = "",
                      values = palette.basis ) +
    labs(x = "x", y = "beta", title = "Beta - inlabru")
  
  data.kappa = cbind(res.inlabru$summary.random$kappa,
                     kappa.true = underlying.effects$kappa.true[res.inlabru$summary.random$kappa$ID + 1])
  p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = kappa.true, color = "True value", fill = "True value")) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
    scale_color_manual(name = "",
                       values = palette.basis ) +
    scale_fill_manual(name = "",
                      values = palette.basis ) +
    labs(x = "t", y = "kappa", title = "Kappa - inlabru")
  
  p.phi <- ggplot(data.frame(res.inlabru$marginals.fixed)) + 
    geom_area(aes(x = phi.x, y = phi.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[2], color = "Estimated", fill="Estimated")) + 
    geom_vline(aes(xintercept = underlying.effects$phi.true, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette.basis) + 
    scale_fill_manual(name = " ", values = palette.basis) +
    labs(x = "Value of phi", y = " ", title = "Phi - inlabru")
  
  p.intercept <- ggplot(data.frame(res.inlabru$marginals.fixed)) + 
    geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Estimated", fill="Estimated")) + 
    #geom_vline(aes(xintercept = underlying.effects$.true, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette.basis) + 
    scale_fill_manual(name = " ", values = palette.basis) +
    labs(x = "Value of phi", y = " ", title = "Intercept - inlabru")
  
  data.eta <- data.frame(eta.sim = res.inlabru$summary.linear.predictor$mean[1:length(obs$eta)]) %>%
    mutate(true.eta = obs$eta) %>%
    mutate(xt = obs$xt, x = obs$x, t = obs$t)
  
  p.eta <- ggplot(data = data.eta) +
    geom_point(aes(x = eta.sim, y = true.eta), color = palette.basis[1]) + 
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
  
  palette.basis <- c('#70A4D4', '#ECC64B', '#93AD80', '#da9124', '#696B8D',
                     '#3290c1',
                     '#5d8060', '#D7B36A', '#826133', '#A85150')
  
  obs <- underlying.effects$obs
  
  data.alpha = cbind(res.inlabru$summary.random$alpha,
                     alpha.true = underlying.effects$alpha.true[res.inlabru$summary.random$alpha$ID])
  p.alpha <- ggplot(data = data.alpha, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = alpha.true, color = "True value", fill = "True value")) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
    geom_point(aes(y = mean + res.inlabru$summary.fixed$mean[1], fill = "Estimated + intercept", color = "Estimated + intercept")) +
    scale_color_manual(name = "",
                       values = palette.basis ) +
    scale_fill_manual(name = "",
                      values = palette.basis ) +
    labs(title="Alpha - inlabru", x = "x", y='')
  
  data.beta = cbind(res.inlabru$summary.random$beta,
                    beta.true = underlying.effects$beta.true[res.inlabru$summary.random$beta$ID])
  p.beta <- ggplot(data = data.beta, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = beta.true, color = "True value", fill = "True value")) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
    scale_color_manual(name = "",
                       values = palette.basis ) +
    scale_fill_manual(name = "",
                      values = palette.basis ) +
    labs(x = "x", y = "beta", title = "Beta - inlabru")
  
  data.kappa = cbind(res.inlabru$summary.random$kappa,
                     kappa.true = underlying.effects$kappa.true[res.inlabru$summary.random$kappa$ID])
  p.kappa <- ggplot(data = data.kappa, aes(x = ID)) + 
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`, fill = "Estimated"), alpha = 0.4) + 
    geom_point(aes(y = kappa.true, color = "True value", fill = "True value")) + 
    geom_point(aes(y = mean, color = "Estimated", fill = "Estimated")) + 
    scale_color_manual(name = "",
                       values = palette.basis ) +
    scale_fill_manual(name = "",
                      values = palette.basis ) +
    labs(x = "t", y = "kappa", title = "Kappa - inlabru")
  
  p.phi <- ggplot(data.frame(res.inlabru$marginals.fixed)) + 
    geom_area(aes(x = phi.x, y = phi.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[2], color = "Estimated", fill="Estimated")) + 
    geom_vline(aes(xintercept = underlying.effects$phi.true, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette.basis) + 
    scale_fill_manual(name = " ", values = palette.basis) +
    labs(x = "Value of phi", y = " ", title = "Phi - inlabru")
  
  p.intercept <- ggplot(data.frame(res.inlabru$marginals.fixed)) + 
    geom_area(aes(x = Int.x, y = Int.y, fill = "Estimated"), alpha = 0.4) + 
    geom_vline(data = res.inlabru$summary.fixed, aes(xintercept = mean[1], color = "Estimated", fill="Estimated")) + 
    #geom_vline(aes(xintercept = underlying.effects$.true, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette.basis) + 
    scale_fill_manual(name = " ", values = palette.basis) +
    labs(x = "Value of phi", y = " ", title = "Intercept - inlabru")
  
  data.eta <- data.frame(eta.sim = res.inlabru$summary.linear.predictor$mean[1:length(obs$eta)]) %>%
    mutate(true.eta = obs$eta) %>%
    mutate(xt = obs$xt, x = obs$x, t = obs$t)
  
  p.eta <- ggplot(data = data.eta) +
    geom_point(aes(x = eta.sim, y = true.eta), color = palette.basis[1]) + 
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