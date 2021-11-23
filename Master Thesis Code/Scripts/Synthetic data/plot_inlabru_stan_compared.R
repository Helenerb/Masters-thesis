# script contatining functions for comparing inlabru and stan results

library(ggplot2)
library(patchwork)

# source relevant functions
# source("../Functions/plotters.R")
# source("../Misc/palette.R")

# assume wokring directory at ---Master Thesis Code
source("Scripts/Functions/plotters.R")
source("Scripts/Misc/palette.R")



plot.inlabru.stan.compared.rw2 <- function(stan.summaries,
                                                  stan.marginals,
                                                  inlabru.summaries,
                                                  res.inlabru,
                                                  underlying.effects,
                                                  cohort=TRUE,
                                           tau.alpha.cutoff = 20,
                                           tau.beta.cutoff = 50000,
                                           tau.kappa.cutoff = 50000,
                                           tau.epsilon.cutoff = 1500
                                           ){
  #' Produces plots with comparison of estimation results from inlabru and STAN
  #' 
  #'@param stan.summaries (list<data.frame>) summaries of STAN results
  #'@param stan.marginals (list<array>) Hamiltonian MC samples from STAN
  #'@param inlabru.summaries (list<data.frame>) summaries of inlabru results
  #'@param res.inlabru (bru object) raw inlabru results
  #'@param underlying.effects (list) underlying data for which analysis is run
  #'@param cohort (boolean) whether analysis includes cohort effect 
  #'
  obs <- underlying.effects$obs
  
  intercept.marginal <- data.frame(int = stan.marginals$intercept_draws)
  
  inlabru.data.fixed = data.frame(res.inlabru$marginals.fixed)
  
  #  ----   intercept   ----
  if(length(stan.marginals$intercept_draws) > 0){
    p.intercept <- ggplot() + 
      geom_histogram(data = intercept.marginal, aes(x = int, y = after_stat(density), color = "Stan", fill = "Stan"), bins=200, alpha = 0.5) + 
      geom_area(data=inlabru.data.fixed, aes(x = Int.x, y = Int.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.4, size = 0.5) + 
      geom_vline(aes(xintercept = underlying.effects$age.intercept.true, color="True", fill="True")) +
      scale_color_manual(name = " ", values = palette) + 
      scale_fill_manual(name = " ", values = palette) +
      theme_classic() + 
      labs(x = "Value of intercept", y = " ", title = "Intercept")
  } else {
    p.intercept <- ggplot(data = data.frame(a = 1, b = 2)) + geom_point(aes(x = a, y = b)) + labs(title = "no intercept data for stan")
  }
  
  
  # ---   alpha   ----
  p.alpha <- ggplot() + 
    geom_ribbon(data=inlabru.summaries$data.alpha, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
    geom_point(data=inlabru.summaries$data.alpha, aes(x = ID, y = mean, color = "Inlabru", fill = "Inlabru"), size = 0.5) + 
    geom_point(data=stan.summaries$summary_alpha, aes(x=index - 1, y=mean, fill="Stan", color="Stan"), size=0.5) + 
    geom_line(data=stan.summaries$summary_alpha, aes(x=index - 1, y=`2.5%`, fill = "Stan", color="Stan"), alpha=0.5) + 
    geom_line(data=stan.summaries$summary_alpha, aes(x=index - 1, y=`97.5%`, fill="Stan", color="Stan"), alpha=0.5) +
    
    geom_point(data=inlabru.summaries$data.alpha, 
               aes(x = ID, y = underlying.effects$alpha.true, color = "True", fill = "True"), size = 0.5) + 
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(title="Alpha", x = "x", y='')
  
  # ---   beta   ----
  
  p.beta <- ggplot() + 
    geom_errorbar(data = inlabru.summaries$data.beta, aes(x = ID + 0.1, ymin = `0.025quant`, ymax = `0.975quant`, color = "Inlabru", fill = "Inlabru")) +
    #geom_point(data=inlabru.summaries$data.beta, aes(x = ID, y = mean, color = "Inlabru", fill = "Inlabru"), size = 0.5) +
    #geom_point(data=stan.summaries$summary_beta, aes(x=index - 1, y=mean, fill="Stan", color="Stan"), size=0.5) + 
    geom_errorbar(data = stan.summaries$summary_beta, aes(x = index - 1 - 0.1, ymin = `2.5%`, ymax = `97.5%`, color = "Stan", fill = "Stan")) +
    geom_point(data = inlabru.summaries$data.beta, aes(x = ID, y = underlying.effects$beta.true, color = "True", fill = "True"), shape = 4) +
    
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    scale_shape_manual(name = "") + 
    theme_classic() + 
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
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(title="Kappa", x = "t", y='')
  
  
  #   ----   gamma   ----  
  if (cohort){
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
      theme_classic() + 
      labs(title="Gamma", x = "c", y='')
  }
  
  #   ----   eta   ----
  p.eta <- ggplot() +
    geom_point(data=inlabru.summaries$data.eta, aes(x = eta.sim, y = true.eta, color = "Inlabru")) + 
    geom_point(data=stan.summaries$summary_eta, aes(x = mean, y = true_eta, color = "Stan")) + 
    scale_color_manual(name = " ", values = palette) + 
    theme_classic() + 
    labs(x="Estimated eta", y="True value for eta", title = "Eta")
  
  p.eta.2 <- ggplot() +
    geom_point(data = inlabru.summaries$data.eta, aes(x=xt, y = eta.sim, color="Inlabru", fill="Inlabru"), size=0.5) +
    #geom_ribbon(data = inlabru.summaries$data.eta, aes(x = xt, ymin = `0.025quant`, ymax=`0.975quant`, fill = "Inlabru"), alpha = 0.5)  +
    geom_ribbon(data = inlabru.summaries$data.eta, aes(x = xt, ymin = X0.025quant, ymax=X0.975quant, fill = "Inlabru"), alpha = 0.5)  +
    geom_point(data = stan.summaries$summary_eta, aes(x=xt, y = mean, color="Stan", fill="Stan"), size=0.5) +
    geom_ribbon(data = stan.summaries$summary_eta, aes(x = xt, ymin = `2.5%`, ymax=`97.5%`, fill="Stan"), alpha=0.5)  +
    geom_line(data = inlabru.summaries$data.eta, aes(x=xt, y = true.eta, color="True", fill="True")) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(x=" ", y="Eta", title="Eta")
  
  p.eta.t <- ggplot() + 
    geom_point(data = inlabru.summaries$data.eta, aes(x=x, y = eta.sim, color = "Inlabru", fill="Inlabru"), size=0.5) +
    #geom_ribbon(data = inlabru.summaries$data.eta, aes(x = x, ymin = `0.025quant`, ymax=`0.975quant`, fill = "Inlabru"), alpha = 0.5)  +
    geom_ribbon(data = inlabru.summaries$data.eta, aes(x = x, ymin = X0.025quant, ymax=X0.975quant, fill = "Inlabru"), alpha = 0.5)  +
    geom_point(data = stan.summaries$summary_eta, aes(x = x, y = mean, color = "Stan", fill="Stan"), size=0.5) +
    geom_ribbon(data = stan.summaries$summary_eta, aes(x = x, ymin = `2.5%`, ymax=`97.5%`, fill="Stan"), alpha=0.5)  +
    geom_line(data = inlabru.summaries$data.eta, aes(x = x, y = true.eta, color = "True", fill="True")) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(x = " ", y = " ", title = "Eta - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.eta.x <- ggplot() + 
    geom_point(data = inlabru.summaries$data.eta, aes(x = t, y = eta.sim, color = "Inlabru", fill="Inlabru"), size=0.5) +
    #geom_ribbon(data = inlabru.summaries$data.eta, aes(x = t, ymin = `0.025quant`, ymax=`0.975quant`, fill = "Inlabru"), alpha = 0.5)  +
    geom_ribbon(data = inlabru.summaries$data.eta, aes(x = t, ymin = X0.025quant, ymax=X0.975quant, fill = "Inlabru"), alpha = 0.5)  +
    geom_point(data=stan.summaries$summary_eta, aes(x = t, y = mean, color = "Stan", fill="Stan"), size=0.5) +
    geom_ribbon(data = stan.summaries$summary_eta, aes(x = t, ymin = `2.5%`, ymax=`97.5%`, fill="Stan"), alpha=0.5)  +
    geom_line(data = inlabru.summaries$data.eta, aes(x = t, y = true.eta, color = "True", fill="True")) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(x = " ", y = " ", title = "Eta - inlabru, for each age") + 
    facet_wrap(~x)
  
  #   ----   hyperparameters: precisions   ----   
  
  
  #  tau alpha 
  
  tau.alpha.stan <- data.frame(tau = stan.marginals$tau_alpha_draws) %>%
    filter(tau < tau.alpha.cutoff)
  tau.alpha.inlabru <- data.frame(res.inlabru$marginals.hyperpar$`Precision for alpha`) %>% 
    filter(x < tau.alpha.cutoff)
  
  p.tau.alpha <- ggplot() + 
    geom_area(data = tau.alpha.inlabru, aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = tau.alpha.stan, aes(x = tau, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins=100) + 
    geom_vline(data = tau.alpha.inlabru, aes(xintercept = underlying.effects$tau.alpha.true, color = "Observed", fill  = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    theme_classic() + 
    labs(x = "Value of precision of alpha", y = " ", title = "Precision of Alpha")
  
  #  tau beta
  
  tau.beta.stan <- data.frame(tau = stan.marginals$tau_beta_draws) %>%
    filter(tau < tau.beta.cutoff)
  tau.beta.inlabru <- data.frame(res.inlabru$marginals.hyperpar$`Precision for beta`) %>%
    filter(x < tau.beta.cutoff)
  
  p.tau.beta <- ggplot() + 
    geom_area(data = tau.beta.inlabru, aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = tau.beta.stan, aes(x = tau, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins=100) + 
    geom_vline(data = tau.beta.inlabru, aes(xintercept = underlying.effects$tau.beta.true, color = "Observed", fill = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    theme_classic() + 
    labs(x = "Value of precision of beta", y = " ", title = "Precision of Beta")
  
  # tau kappa
  tau.kappa.stan <- data.frame(tau = stan.marginals$tau_kappa_draws) %>%
    filter(tau < tau.kappa.cutoff)
  tau.kappa.inlabru <- data.frame(res.inlabru$marginals.hyperpar$`Precision for kappa`) %>%
    filter(x < tau.kappa.cutoff)
  
  p.tau.kappa <- ggplot() + 
    geom_area(data = tau.kappa.inlabru, aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = tau.kappa.stan, aes(x = tau, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins = 100) + 
    geom_vline(data = tau.kappa.inlabru, aes(xintercept = underlying.effects$tau.kappa.true, color = "Observed", fill = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    theme_classic() + 
    labs(x = "Value of precision of kappa", y = " ", title = "Precision of Kappa")
  
  if (cohort){
    # tau gamma
    tau.gamma.stan <- data.frame(tau = stan.marginals$tau_gamma_draws) %>%
      filter(tau < tau.gamma.cutoff)
    tau.gamma.inlabru <- data.frame(res.inlabru$marginals.hyperpar$`Precision for gamma`) %>%
      filter(x < tau.gamma.cutoff)
    
    p.tau.gamma <- ggplot() + 
      geom_area(data = tau.gamma.inlabru, aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
      geom_histogram(data = tau.gamma.stan, aes(x = tau, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins = 100) + 
      geom_vline(data = tau.gamma.inlabru, aes(xintercept = underlying.effects$tau.gamma.true, color = "Observed", fill = "Observed")) + 
      scale_color_manual(name = " ", values = palette) + 
      scale_fill_manual(name = " ", values = palette) +
      theme_classic() + 
      labs(x = "Value of precision of gamma", y = " ", title = "Precision of Gamma")
  }
  
  # tau epsilon
  if(length(stan.marginals$tau_epsilon_draws) > 0){
    tau.epsilon.stan <- data.frame(tau = stan.marginals$tau_epsilon_draws) %>%
      filter(tau < tau.epsilon.cutoff)
    tau.epsilon.inlabru <- data.frame(res.inlabru$marginals.hyperpar$`Precision for epsilon`) %>%
      filter(x < tau.epsilon.cutoff)
    
    p.tau.epsilon <- ggplot() + 
      geom_area(data = tau.epsilon.inlabru, aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
      geom_histogram(data = tau.epsilon.stan, aes(x = tau, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins = 100) + 
      geom_vline(data = tau.epsilon.inlabru, aes(xintercept = underlying.effects$tau.epsilon.true, color = "Observed", fill = "Observed")) + 
      scale_color_manual(name = " ", values = palette) + 
      scale_fill_manual(name = " ", values = palette) +
      theme_classic() + 
      labs(x = "Value of precision of epsilon", y = " ", title = "Precision of Epsilon")
  } else {
    print("No tau epsilon!")
    dummy.data  <- data.frame(a = c(1.0, 2.0), b = c(1.0, 2.0))
    print(dummy.data)
    p.tau.epsilon <- ggplot(data = dummy.data) + geom_point(aes(x = a, y = b)) + labs(title = "No available tau epsilon")
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
                p.tau.epsilon = p.tau.epsilon)
  if(cohort){
    plots <- c(plots, p.gamma=p.gamma)
    plots <- c(plots, p.tau.gamma=p.tau.gamma)
  }
  
  return(plots)
}

plot.inlabru.stan.traditional.lc <- function(stan.summaries,
                                           stan.marginals,
                                           inlabru.summaries,
                                           res.inlabru,
                                           underlying.effects,
                                           cohort=TRUE,
                                           tau.alpha.cutoff = 20,
                                           tau.beta.cutoff = 50000,
                                           tau.kappa.cutoff = 50000,
                                           tau.epsilon.cutoff = 1500,
                                           a45=F
){
  #' Produces plots with comparison of estimation results from inlabru and STAN
  #' 
  #'@param stan.summaries (list<data.frame>) summaries of STAN results
  #'@param stan.marginals (list<array>) Hamiltonian MC samples from STAN
  #'@param inlabru.summaries (list<data.frame>) summaries of inlabru results
  #'@param res.inlabru (bru object) raw inlabru results
  #'@param underlying.effects (list) underlying data for which analysis is run
  #'@param cohort (boolean) whether analysis includes cohort effect 
  #'
  obs <- underlying.effects$obs
  
  intercept.marginal <- data.frame(int = stan.marginals$intercept_draws)
  
  inlabru.data.fixed = data.frame(res.inlabru$marginals.fixed)
  
  #  ----   intercept   ----
  p.intercept <- ggplot() + 
    geom_histogram(data = intercept.marginal, aes(x = int, y = after_stat(density), color = "Stan", fill = "Stan"), bins=200, alpha = 0.5) + 
    geom_area(data=inlabru.data.fixed, aes(x = Int.x, y = Int.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.4, size = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$age.intercept.true, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    theme_classic() + 
    labs(x = "Value of intercept", y = " ", title = "Intercept")
  
  # ---   alpha   ----
  if(!a45){
    p.alpha <- ggplot() + 
      geom_ribbon(data=inlabru.summaries$data.alpha, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
      geom_point(data=inlabru.summaries$data.alpha, aes(x = ID, y = mean, color = "Inlabru", fill = "Inlabru"), size = 0.5) + 
      geom_point(data=stan.summaries$summary_alpha, aes(x= index - 1, y=mean, fill="Stan", color="Stan"), size=0.5) + 
      geom_line(data=stan.summaries$summary_alpha, aes(x= index - 1, y=`2.5%`, fill = "Stan", color="Stan"), alpha=0.5) + 
      geom_line(data=stan.summaries$summary_alpha, aes(x= index - 1, y=`97.5%`, fill="Stan", color="Stan"), alpha=0.5) +
      
      geom_point(data=inlabru.summaries$data.alpha, 
                 aes(x = ID, y = underlying.effects$alpha.true, color = "True", fill = "True"), size = 0.5) + 
      scale_color_manual(name = "", values = palette ) +
      scale_fill_manual(name = "", values = palette ) +
      theme_classic() + 
      labs(title="Alpha", x = "x", y='')
  } else {
    p.alpha <- ggplot() + 
      geom_ribbon(data=inlabru.summaries$data.alpha, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
      geom_point(data=inlabru.summaries$data.alpha, aes(x = ID, y = mean, color = "Inlabru", fill = "Inlabru"), size = 0.5) + 
      geom_point(data=stan.summaries$summary_alpha, aes(x= index - 1 + 9, y=mean, fill="Stan", color="Stan"), size=0.5) + 
      geom_line(data=stan.summaries$summary_alpha, aes(x= index - 1 + 9, y=`2.5%`, fill = "Stan", color="Stan"), alpha=0.5) + 
      geom_line(data=stan.summaries$summary_alpha, aes(x= index - 1 + 9, y=`97.5%`, fill="Stan", color="Stan"), alpha=0.5) +
      
      geom_point(data=inlabru.summaries$data.alpha, 
                 aes(x = ID, y = underlying.effects$alpha.true, color = "True", fill = "True"), size = 0.5) + 
      scale_color_manual(name = "", values = palette ) +
      scale_fill_manual(name = "", values = palette ) +
      theme_classic() + 
      labs(title="Alpha", x = "x", y='')
  }
  
  
  # ---   beta   ----
  
  if(!a45){
    p.beta <- ggplot() + 
      geom_errorbar(data = inlabru.summaries$data.beta, aes(x = ID + 0.1, ymin = `0.025quant`, ymax = `0.975quant`, color = "Inlabru", fill = "Inlabru")) +
      geom_errorbar(data = stan.summaries$summary_beta, aes(x = index - 1 - 0.1, ymin = `2.5%`, ymax = `97.5%`, color = "Stan", fill = "Stan")) +
      geom_point(data = inlabru.summaries$data.beta, aes(x = ID, y = underlying.effects$beta.true, color = "True", fill = "True"), shape = 4) +
      
      scale_color_manual(name = "", values = palette ) +
      scale_fill_manual(name = "", values = palette ) +
      scale_shape_manual(name = "") + 
      theme_classic() + 
      labs(title="Beta", x = "x", y='')
  } else {
    p.beta <- ggplot() + 
      geom_errorbar(data = inlabru.summaries$data.beta, aes(x = ID + 0.1, ymin = `0.025quant`, ymax = `0.975quant`, color = "Inlabru", fill = "Inlabru")) +
      geom_errorbar(data = stan.summaries$summary_beta, aes(x = index - 1 - 0.1 + 9, ymin = `2.5%`, ymax = `97.5%`, color = "Stan", fill = "Stan")) +
      geom_point(data = inlabru.summaries$data.beta, aes(x = ID, y = underlying.effects$beta.true, color = "True", fill = "True"), shape = 4) +
      
      scale_color_manual(name = "", values = palette ) +
      scale_fill_manual(name = "", values = palette ) +
      scale_shape_manual(name = "") + 
      theme_classic() + 
      labs(title="Beta", x = "x", y='')
  }
  
  #   ----   kappa   ----  
  p.kappa <- ggplot() + 
    geom_ribbon(data=inlabru.summaries$data.kappa, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
    geom_point(data=inlabru.summaries$data.kappa, aes(x = ID, y = mean, color = "Inlabru", fill = "Inlabru"), size = 0.5) + 
    
    geom_point(data=stan.summaries$summary_kappa, aes(x=index - 1, y=mean, fill="Stan",color="Stan"), size=0.5) + 
    geom_line(data=stan.summaries$summary_kappa, aes(x=index - 1, y=`2.5%`, fill="Stan", color="Stan"), alpha=0.5) + 
    geom_line(data=stan.summaries$summary_kappa, aes(x=index - 1, y=`97.5%`, fill="Stan", color="Stan"), alpha=0.5) +
    
    geom_point(data=inlabru.summaries$data.kappa,
               aes(x = ID, y = underlying.effects$kappa.true, color = "True", fill = "True"), size = 0.5) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(title="Kappa", x = "t", y='')
  
  
  #   ----   gamma   ----  
  if (cohort){
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
      theme_classic() + 
      labs(title="Gamma", x = "c", y='')
  }
  
  #   ----   eta   ----
  p.eta <- ggplot() +
    geom_point(data=inlabru.summaries$data.eta, aes(x = eta.sim, y = true.eta, color = "Inlabru")) + 
    geom_point(data=stan.summaries$summary_eta, aes(x = mean, y = true_eta, color = "Stan")) + 
    scale_color_manual(name = " ", values = palette) + 
    theme_classic() + 
    labs(x="Estimated eta", y="True value for eta", title = "Eta")
  
  p.eta.2 <- ggplot() +
    geom_point(data = inlabru.summaries$data.eta, aes(x=xt, y = eta.sim, color="Inlabru", fill="Inlabru"), size=0.5) +
    #geom_ribbon(data = inlabru.summaries$data.eta, aes(x = xt, ymin = `0.025quant`, ymax=`0.975quant`, fill = "Inlabru"), alpha = 0.5)  +
    geom_ribbon(data = inlabru.summaries$data.eta, aes(x = xt, ymin = X0.025quant, ymax=X0.975quant, fill = "Inlabru"), alpha = 0.5)  +
    geom_point(data = stan.summaries$summary_eta, aes(x=xt, y = mean, color="Stan", fill="Stan"), size=0.5) +
    geom_ribbon(data = stan.summaries$summary_eta, aes(x = xt, ymin = `2.5%`, ymax=`97.5%`, fill="Stan"), alpha=0.5)  +
    geom_line(data = inlabru.summaries$data.eta, aes(x=xt, y = true.eta, color="True", fill="True")) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(x=" ", y="Eta", title="Eta")
  
  p.eta.t <- ggplot() + 
    geom_point(data = inlabru.summaries$data.eta, aes(x=x, y = eta.sim, color = "Inlabru", fill="Inlabru"), size=0.5) +
    #geom_ribbon(data = inlabru.summaries$data.eta, aes(x = x, ymin = `0.025quant`, ymax=`0.975quant`, fill = "Inlabru"), alpha = 0.5)  +
    geom_ribbon(data = inlabru.summaries$data.eta, aes(x = x, ymin = X0.025quant, ymax=X0.975quant, fill = "Inlabru"), alpha = 0.5)  +
    geom_point(data = stan.summaries$summary_eta, aes(x = x, y = mean, color = "Stan", fill="Stan"), size=0.5) +
    geom_ribbon(data = stan.summaries$summary_eta, aes(x = x, ymin = `2.5%`, ymax=`97.5%`, fill="Stan"), alpha=0.5)  +
    geom_line(data = inlabru.summaries$data.eta, aes(x = x, y = true.eta, color = "True", fill="True")) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(x = " ", y = " ", title = "Eta - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.eta.x <- ggplot() + 
    geom_point(data = inlabru.summaries$data.eta, aes(x = t, y = eta.sim, color = "Inlabru", fill="Inlabru"), size=0.5) +
    #geom_ribbon(data = inlabru.summaries$data.eta, aes(x = t, ymin = `0.025quant`, ymax=`0.975quant`, fill = "Inlabru"), alpha = 0.5)  +
    geom_ribbon(data = inlabru.summaries$data.eta, aes(x = t, ymin = X0.025quant, ymax=X0.975quant, fill = "Inlabru"), alpha = 0.5)  +
    geom_point(data=stan.summaries$summary_eta, aes(x = t, y = mean, color = "Stan", fill="Stan"), size=0.5) +
    geom_ribbon(data = stan.summaries$summary_eta, aes(x = t, ymin = `2.5%`, ymax=`97.5%`, fill="Stan"), alpha=0.5)  +
    geom_line(data = inlabru.summaries$data.eta, aes(x = t, y = true.eta, color = "True", fill="True")) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(x = " ", y = " ", title = "Eta - inlabru, for each age") + 
    facet_wrap(~x)
  
  #   ----   predictor - marginals   ----
  
  p.eta.marg.1 <- 
  
  #   ----   hyperparameters: precisions   ----   
  
  
  #  tau alpha 
  
  tau.alpha.stan <- data.frame(tau = stan.marginals$tau_alpha_draws) %>%
    filter(tau < tau.alpha.cutoff)
  tau.alpha.inlabru <- data.frame(res.inlabru$marginals.hyperpar$`Precision for alpha`) %>% 
    filter(x < tau.alpha.cutoff)
  
  p.tau.alpha <- ggplot() + 
    geom_area(data = tau.alpha.inlabru, aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = tau.alpha.stan, aes(x = tau, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins=100) + 
    geom_vline(data = tau.alpha.inlabru, aes(xintercept = underlying.effects$tau.alpha.true, color = "Observed", fill  = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    theme_classic() + 
    labs(x = "Value of precision of alpha", y = " ", title = "Precision of Alpha")
  
  #  tau beta
  
  tau.beta.stan <- data.frame(tau = stan.marginals$tau_beta_draws) %>%
    filter(tau < tau.beta.cutoff)
  tau.beta.inlabru <- data.frame(res.inlabru$marginals.hyperpar$`Precision for beta`) %>%
    filter(x < tau.beta.cutoff)
  
  p.tau.beta <- ggplot() + 
    geom_area(data = tau.beta.inlabru, aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = tau.beta.stan, aes(x = tau, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins=100) + 
    geom_vline(data = tau.beta.inlabru, aes(xintercept = underlying.effects$tau.beta.true, color = "Observed", fill = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    theme_classic() + 
    labs(x = "Value of precision of beta", y = " ", title = "Precision of Beta")
  
  # tau kappa
  tau.kappa.stan <- data.frame(tau = stan.marginals$tau_kappa_draws) %>%
    filter(tau < tau.kappa.cutoff)
  tau.kappa.inlabru <- data.frame(res.inlabru$marginals.hyperpar$`Precision for kappa`) %>%
    filter(x < tau.kappa.cutoff)
  
  p.tau.kappa <- ggplot() + 
    geom_area(data = tau.kappa.inlabru, aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = tau.kappa.stan, aes(x = tau, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins = 100) + 
    geom_vline(data = tau.kappa.inlabru, aes(xintercept = underlying.effects$tau.kappa.true, color = "Observed", fill = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    theme_classic() + 
    labs(x = "Value of precision of kappa", y = " ", title = "Precision of Kappa")
  
  if (cohort){
    # tau gamma
    tau.gamma.stan <- data.frame(tau = stan.marginals$tau_gamma_draws) %>%
      filter(tau < tau.gamma.cutoff)
    tau.gamma.inlabru <- data.frame(res.inlabru$marginals.hyperpar$`Precision for gamma`) %>%
      filter(x < tau.gamma.cutoff)
    
    p.tau.gamma <- ggplot() + 
      geom_area(data = tau.gamma.inlabru, aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
      geom_histogram(data = tau.gamma.stan, aes(x = tau, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins = 100) + 
      geom_vline(data = tau.gamma.inlabru, aes(xintercept = underlying.effects$tau.gamma.true, color = "Observed", fill = "Observed")) + 
      scale_color_manual(name = " ", values = palette) + 
      scale_fill_manual(name = " ", values = palette) +
      theme_classic() + 
      labs(x = "Value of precision of gamma", y = " ", title = "Precision of Gamma")
  }
  
  # tau epsilon
  tau.epsilon.stan <- data.frame(tau = stan.marginals$tau_epsilon_draws) %>%
    filter(tau < tau.epsilon.cutoff)
  tau.epsilon.inlabru <- data.frame(res.inlabru$marginals.hyperpar$`Precision for the Gaussian observations`) %>%
    filter(x < tau.epsilon.cutoff)
  
  p.tau.epsilon <- ggplot() + 
    geom_area(data = tau.epsilon.inlabru, aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = tau.epsilon.stan, aes(x = tau, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins = 100) + 
    geom_vline(data = tau.epsilon.inlabru, aes(xintercept = underlying.effects$tau.epsilon.true, color = "Observed", fill = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    theme_classic() + 
    labs(x = "Value of precision of epsilon", y = " ", title = "Precision of Epsilon")
  
  
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
                p.tau.epsilon = p.tau.epsilon)
  if(cohort){
    plots <- c(plots, p.gamma=p.gamma)
    plots <- c(plots, p.tau.gamma=p.tau.gamma)
  }
  
  return(plots)
}

plot.inlabru.stan.traditional.lc.no.beta <- function(stan.summaries,
                                             stan.marginals,
                                             inlabru.summaries,
                                             res.inlabru,
                                             underlying.effects,
                                             cohort=TRUE,
                                             tau.alpha.cutoff = 20,
                                             tau.beta.cutoff = 50000,
                                             tau.kappa.cutoff = 50000,
                                             tau.epsilon.cutoff = 1500,
                                             a45=F
){
  #' Produces plots with comparison of estimation results from inlabru and STAN
  #' 
  #'@param stan.summaries (list<data.frame>) summaries of STAN results
  #'@param stan.marginals (list<array>) Hamiltonian MC samples from STAN
  #'@param inlabru.summaries (list<data.frame>) summaries of inlabru results
  #'@param res.inlabru (bru object) raw inlabru results
  #'@param underlying.effects (list) underlying data for which analysis is run
  #'@param cohort (boolean) whether analysis includes cohort effect 
  #'
  obs <- underlying.effects$obs
  
  intercept.marginal <- data.frame(int = stan.marginals$intercept_draws)
  
  inlabru.data.fixed = data.frame(res.inlabru$marginals.fixed)
  
  #  ----   intercept   ----
  p.intercept <- ggplot() + 
    geom_histogram(data = intercept.marginal, aes(x = int, y = after_stat(density), color = "Stan", fill = "Stan"), bins=200, alpha = 0.5) + 
    geom_area(data=inlabru.data.fixed, aes(x = Int.x, y = Int.y, color = "Inlabru", fill = "Inlabru"), alpha = 0.4, size = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$age.intercept.true, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    theme_classic() + 
    labs(x = "Value of intercept", y = " ", title = "Intercept")
  
  # ---   alpha   ----
  if(!a45){
    p.alpha <- ggplot() + 
      geom_ribbon(data=inlabru.summaries$data.alpha, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
      geom_point(data=inlabru.summaries$data.alpha, aes(x = ID, y = mean, color = "Inlabru", fill = "Inlabru"), size = 0.5) + 
      geom_point(data=stan.summaries$summary_alpha, aes(x= index - 1, y=mean, fill="Stan", color="Stan"), size=0.5) + 
      geom_line(data=stan.summaries$summary_alpha, aes(x= index - 1, y=`2.5%`, fill = "Stan", color="Stan"), alpha=0.5) + 
      geom_line(data=stan.summaries$summary_alpha, aes(x= index - 1, y=`97.5%`, fill="Stan", color="Stan"), alpha=0.5) +
      
      geom_point(data=inlabru.summaries$data.alpha, 
                 aes(x = ID, y = underlying.effects$alpha.true, color = "True", fill = "True"), size = 0.5) + 
      scale_color_manual(name = "", values = palette ) +
      scale_fill_manual(name = "", values = palette ) +
      theme_classic() + 
      labs(title="Alpha", x = "x", y='')
  } else {
    p.alpha <- ggplot() + 
      geom_ribbon(data=inlabru.summaries$data.alpha, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
      geom_point(data=inlabru.summaries$data.alpha, aes(x = ID, y = mean, color = "Inlabru", fill = "Inlabru"), size = 0.5) + 
      geom_point(data=stan.summaries$summary_alpha, aes(x= index - 1 + 9, y=mean, fill="Stan", color="Stan"), size=0.5) + 
      geom_line(data=stan.summaries$summary_alpha, aes(x= index - 1 + 9, y=`2.5%`, fill = "Stan", color="Stan"), alpha=0.5) + 
      geom_line(data=stan.summaries$summary_alpha, aes(x= index - 1 + 9, y=`97.5%`, fill="Stan", color="Stan"), alpha=0.5) +
      
      geom_point(data=inlabru.summaries$data.alpha, 
                 aes(x = ID, y = underlying.effects$alpha.true, color = "True", fill = "True"), size = 0.5) + 
      scale_color_manual(name = "", values = palette ) +
      scale_fill_manual(name = "", values = palette ) +
      theme_classic() + 
      labs(title="Alpha", x = "x", y='')
  }
  
  
  # ---   beta   ----
  
  # if(!a45){
  #   p.beta <- ggplot() + 
  #     geom_errorbar(data = inlabru.summaries$data.beta, aes(x = ID + 0.1, ymin = `0.025quant`, ymax = `0.975quant`, color = "Inlabru", fill = "Inlabru")) +
  #     geom_errorbar(data = stan.summaries$summary_beta, aes(x = index - 1 - 0.1, ymin = `2.5%`, ymax = `97.5%`, color = "Stan", fill = "Stan")) +
  #     geom_point(data = inlabru.summaries$data.beta, aes(x = ID, y = underlying.effects$beta.true, color = "True", fill = "True"), shape = 4) +
  #     
  #     scale_color_manual(name = "", values = palette ) +
  #     scale_fill_manual(name = "", values = palette ) +
  #     scale_shape_manual(name = "") + 
  #     theme_classic() + 
  #     labs(title="Beta", x = "x", y='')
  # } else {
  #   p.beta <- ggplot() + 
  #     geom_errorbar(data = inlabru.summaries$data.beta, aes(x = ID + 0.1, ymin = `0.025quant`, ymax = `0.975quant`, color = "Inlabru", fill = "Inlabru")) +
  #     geom_errorbar(data = stan.summaries$summary_beta, aes(x = index - 1 - 0.1 + 9, ymin = `2.5%`, ymax = `97.5%`, color = "Stan", fill = "Stan")) +
  #     geom_point(data = inlabru.summaries$data.beta, aes(x = ID, y = underlying.effects$beta.true, color = "True", fill = "True"), shape = 4) +
  #     
  #     scale_color_manual(name = "", values = palette ) +
  #     scale_fill_manual(name = "", values = palette ) +
  #     scale_shape_manual(name = "") + 
  #     theme_classic() + 
  #     labs(title="Beta", x = "x", y='')
  # }
  
  p.beta <- ggplot(data = data.frame(x = 1, y = 2)) + geom_point(aes(x = x, y = y)) + labs(title = "No beta available")
  
  #   ----   kappa   ----  
  p.kappa <- ggplot() + 
    geom_ribbon(data=inlabru.summaries$data.kappa, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
    geom_point(data=inlabru.summaries$data.kappa, aes(x = ID, y = mean, color = "Inlabru", fill = "Inlabru"), size = 0.5) + 
    
    geom_point(data=stan.summaries$summary_kappa, aes(x=index - 1, y=mean, fill="Stan",color="Stan"), size=0.5) + 
    geom_line(data=stan.summaries$summary_kappa, aes(x=index - 1, y=`2.5%`, fill="Stan", color="Stan"), alpha=0.5) + 
    geom_line(data=stan.summaries$summary_kappa, aes(x=index - 1, y=`97.5%`, fill="Stan", color="Stan"), alpha=0.5) +
    
    # geom_point(data=inlabru.summaries$data.kappa,
    #            aes(x = ID, y = underlying.effects$kappa.true, color = "True", fill = "True"), size = 0.5) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(title="Kappa", x = "t", y='')
  
  
  #   ----   gamma   ----  
  if (cohort){
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
      theme_classic() + 
      labs(title="Gamma", x = "c", y='')
  }
  
  #   ----   eta   ----
  p.eta <- ggplot() +
    geom_point(data=inlabru.summaries$data.eta, aes(x = eta.sim, y = true.eta, color = "Inlabru")) + 
    geom_point(data=stan.summaries$summary_eta, aes(x = mean, y = true_eta, color = "Stan")) + 
    scale_color_manual(name = " ", values = palette) + 
    theme_classic() + 
    labs(x="Estimated eta", y="True value for eta", title = "Eta")
  
  p.eta.2 <- ggplot() +
    geom_point(data = inlabru.summaries$data.eta, aes(x=xt, y = eta.sim, color="Inlabru", fill="Inlabru"), size=0.5) +
    #geom_ribbon(data = inlabru.summaries$data.eta, aes(x = xt, ymin = `0.025quant`, ymax=`0.975quant`, fill = "Inlabru"), alpha = 0.5)  +
    geom_ribbon(data = inlabru.summaries$data.eta, aes(x = xt, ymin = X0.025quant, ymax=X0.975quant, fill = "Inlabru"), alpha = 0.5)  +
    geom_point(data = stan.summaries$summary_eta, aes(x=xt, y = mean, color="Stan", fill="Stan"), size=0.5) +
    geom_ribbon(data = stan.summaries$summary_eta, aes(x = xt, ymin = `2.5%`, ymax=`97.5%`, fill="Stan"), alpha=0.5)  +
    geom_line(data = inlabru.summaries$data.eta, aes(x=xt, y = true.eta, color="True", fill="True")) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(x=" ", y="Eta", title="Eta")
  
  p.eta.t <- ggplot() + 
    geom_point(data = inlabru.summaries$data.eta, aes(x=x, y = eta.sim, color = "Inlabru", fill="Inlabru"), size=0.5) +
    #geom_ribbon(data = inlabru.summaries$data.eta, aes(x = x, ymin = `0.025quant`, ymax=`0.975quant`, fill = "Inlabru"), alpha = 0.5)  +
    geom_ribbon(data = inlabru.summaries$data.eta, aes(x = x, ymin = X0.025quant, ymax=X0.975quant, fill = "Inlabru"), alpha = 0.5)  +
    geom_point(data = stan.summaries$summary_eta, aes(x = x, y = mean, color = "Stan", fill="Stan"), size=0.5) +
    geom_ribbon(data = stan.summaries$summary_eta, aes(x = x, ymin = `2.5%`, ymax=`97.5%`, fill="Stan"), alpha=0.5)  +
    geom_line(data = inlabru.summaries$data.eta, aes(x = x, y = true.eta, color = "True", fill="True")) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(x = " ", y = " ", title = "Eta - inlabru, for each year") + 
    facet_wrap(~t)
  
  p.eta.x <- ggplot() + 
    geom_point(data = inlabru.summaries$data.eta, aes(x = t, y = eta.sim, color = "Inlabru", fill="Inlabru"), size=0.5) +
    #geom_ribbon(data = inlabru.summaries$data.eta, aes(x = t, ymin = `0.025quant`, ymax=`0.975quant`, fill = "Inlabru"), alpha = 0.5)  +
    geom_ribbon(data = inlabru.summaries$data.eta, aes(x = t, ymin = X0.025quant, ymax=X0.975quant, fill = "Inlabru"), alpha = 0.5)  +
    geom_point(data=stan.summaries$summary_eta, aes(x = t, y = mean, color = "Stan", fill="Stan"), size=0.5) +
    geom_ribbon(data = stan.summaries$summary_eta, aes(x = t, ymin = `2.5%`, ymax=`97.5%`, fill="Stan"), alpha=0.5)  +
    geom_line(data = inlabru.summaries$data.eta, aes(x = t, y = true.eta, color = "True", fill="True")) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    theme_classic() + 
    labs(x = " ", y = " ", title = "Eta - inlabru, for each age") + 
    facet_wrap(~x)
  
    
  #   ----   hyperparameters: precisions   ----   
  
  
  #  tau alpha 
  
  tau.alpha.stan <- data.frame(tau = stan.marginals$tau_alpha_draws) %>%
    filter(tau < tau.alpha.cutoff)
  tau.alpha.inlabru <- data.frame(res.inlabru$marginals.hyperpar$`Precision for alpha`) %>% 
    filter(x < tau.alpha.cutoff)
  
  p.tau.alpha <- ggplot() + 
    geom_area(data = tau.alpha.inlabru, aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = tau.alpha.stan, aes(x = tau, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins=100) + 
    geom_vline(data = tau.alpha.inlabru, aes(xintercept = underlying.effects$tau.alpha.true, color = "Observed", fill  = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    theme_classic() + 
    labs(x = "Value of precision of alpha", y = " ", title = "Precision of Alpha")
  
  #  tau beta
  
  # tau.beta.stan <- data.frame(tau = stan.marginals$tau_beta_draws) %>%
  #   filter(tau < tau.beta.cutoff)
  # tau.beta.inlabru <- data.frame(res.inlabru$marginals.hyperpar$`Precision for beta`) %>%
  #   filter(x < tau.beta.cutoff)
  # 
  # p.tau.beta <- ggplot() + 
  #   geom_area(data = tau.beta.inlabru, aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  #   geom_histogram(data = tau.beta.stan, aes(x = tau, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins=100) + 
  #   geom_vline(data = tau.beta.inlabru, aes(xintercept = underlying.effects$tau.beta.true, color = "Observed", fill = "Observed")) + 
  #   scale_color_manual(name = " ", values = palette) + 
  #   scale_fill_manual(name = " ", values = palette) +
  #   theme_classic() + 
  #   labs(x = "Value of precision of beta", y = " ", title = "Precision of Beta")
  
  p.tau.beta <- ggplot(data = data.frame(x = 1, y = 2)) + geom_point(aes(x = x, y = y)) + labs(title = "No beta available")
  
  # tau kappa
  tau.kappa.stan <- data.frame(tau = stan.marginals$tau_kappa_draws) %>%
    filter(tau < tau.kappa.cutoff)
  tau.kappa.inlabru <- data.frame(res.inlabru$marginals.hyperpar$`Precision for kappa`) %>%
    filter(x < tau.kappa.cutoff)
  
  p.tau.kappa <- ggplot() + 
    geom_area(data = tau.kappa.inlabru, aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = tau.kappa.stan, aes(x = tau, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins = 100) + 
    geom_vline(data = tau.kappa.inlabru, aes(xintercept = underlying.effects$tau.kappa.true, color = "Observed", fill = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    theme_classic() + 
    labs(x = "Value of precision of kappa", y = " ", title = "Precision of Kappa")
  
  if (cohort){
    # tau gamma
    tau.gamma.stan <- data.frame(tau = stan.marginals$tau_gamma_draws) %>%
      filter(tau < tau.gamma.cutoff)
    tau.gamma.inlabru <- data.frame(res.inlabru$marginals.hyperpar$`Precision for gamma`) %>%
      filter(x < tau.gamma.cutoff)
    
    p.tau.gamma <- ggplot() + 
      geom_area(data = tau.gamma.inlabru, aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
      geom_histogram(data = tau.gamma.stan, aes(x = tau, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins = 100) + 
      geom_vline(data = tau.gamma.inlabru, aes(xintercept = underlying.effects$tau.gamma.true, color = "Observed", fill = "Observed")) + 
      scale_color_manual(name = " ", values = palette) + 
      scale_fill_manual(name = " ", values = palette) +
      theme_classic() + 
      labs(x = "Value of precision of gamma", y = " ", title = "Precision of Gamma")
  }
  
  # tau epsilon
  tau.epsilon.stan <- data.frame(tau = stan.marginals$tau_epsilon_draws) %>%
    filter(tau < tau.epsilon.cutoff)
  tau.epsilon.inlabru <- data.frame(res.inlabru$marginals.hyperpar$`Precision for the Gaussian observations`) %>%
    filter(x < tau.epsilon.cutoff)
  
  p.tau.epsilon <- ggplot() + 
    geom_area(data = tau.epsilon.inlabru, aes(x = x, y = y, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
    geom_histogram(data = tau.epsilon.stan, aes(x = tau, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins = 100) + 
    geom_vline(data = tau.epsilon.inlabru, aes(xintercept = underlying.effects$tau.epsilon.true, color = "Observed", fill = "Observed")) + 
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    theme_classic() + 
    labs(x = "Value of precision of epsilon", y = " ", title = "Precision of Epsilon")
  
  
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
                p.tau.epsilon = p.tau.epsilon)
  if(cohort){
    plots <- c(plots, p.gamma=p.gamma)
    plots <- c(plots, p.tau.gamma=p.tau.gamma)
  }
  
  return(plots)
}



save.compared.rw2.lc <- function(plots, path.to.storage){
  p.random.effects <- (plots$p.intercept | plots$p.alpha)/(plots$p.beta | plots$p.kappa) + 
    plot_layout(guides="collect")
  save.figure(p.random.effects, name = "random_effects_comparison", path = path.to.storage)
  
  p.eta.xt <- (plots$p.eta | plots$p.eta.2) + plot_layout(guides = "collect")
  save.figure(p.eta.xt, name = "eta_xt_comparison", path = path.to.storage)
  
  p.eta.facet <- (plots$p.eta.x | plots$p.eta.t) + plot_layout(guides = "collect")
  save.figure(p.eta.facet, name = "eta_facet_comparison", path = path.to.storage)
}

save.compared.rw2.cohort <- function(plots, path.to.storage){
  p.random.effects <- (plots$p.intercept | plots$p.alpha | plots$p.beta )/(plots$p.kappa | plots$p.gamma) + 
    plot_layout(guides="collect")
  save.figure(p.random.effects, name = "random_effects_comparison", path = path.to.storage)
  
  p.eta.xt <- (plots$p.eta | plots$p.eta.2) + plot_layout(guides = "collect")
  save.figure(p.eta.xt, name = "eta_xt_comparison", path = path.to.storage)
  
  p.eta.facet <- (plots$p.eta.x | plots$p.eta.t) + plot_layout(guides = "collect")
  save.figure(p.eta.facet, name = "eta_facet_comparison", path = path.to.storage)
}



save.compared.rw2 <- function(plots, path.to.storage, cohort=TRUE){
  if(cohort){
    p.random.effects <- (plots$p.intercept | plots$p.alpha | plots$p.beta )/(plots$p.kappa | plots$p.gamma) + 
      plot_layout(guides="collect")
  } else {
    p.random.effects <- (plots$p.intercept | plots$p.alpha )/(plots$p.beta | plots$p.kappa) + 
      plot_layout(guides="collect")
  }
  
  save.figure(p.random.effects, name = "random_effects_comparison", path = path.to.storage)
  
  save.figure(plots$p.intercept, name = "intercept_comparison", path = path.to.storage)
  
  if(cohort){
    p.hypers <- (plots$p.tau.alpha | plots$p.tau.beta) / (plots$p.tau.kappa | plots$p.tau.gamma | plots$p.tau.epsilon) + plot_layout(guides = "collect")
  }
  else{
    p.hypers <- (plots$p.tau.alpha | plots$p.tau.beta) / (plots$p.tau.kappa | plots$p.tau.epsilon) + plot_layout(guides = "collect")
  }
  save.figure(p.hypers, name = "hypers_comparison", path = path.to.storage)
  
  p.eta.xt <- (plots$p.eta | plots$p.eta.2) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  save.figure(p.eta.xt, name = "eta_xt_comparison", path = path.to.storage)
  
  p.eta.facet <- (plots$p.eta.x | plots$p.eta.t) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  save.figure(p.eta.facet, name = "eta_facet_comparison", path = path.to.storage)
  
  save.figure(plots$p.eta.x, "eta_x_comparison", path = path.to.storage)
  save.figure(plots$p.eta.t, "eta_t_comparison", path = path.to.storage)
}

produce.compared.plots <- function(
  stan.summaries, stan.marginals, inlabru.summaries, res.inlabru, underlying.effects, plot.func, save.func,
  path.to.storage){
  
  plots <- plot.func(stan.summaries, stan.marginals, inlabru.summaries, res.inlabru, underlying.effects)
  save.func(plots, path.to.storage)
  return(plots)
}