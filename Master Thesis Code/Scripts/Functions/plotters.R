# script containing functions for plotting. 
library(ggplot2)
library(patchwork)

# assuming working directory at .../Master Thesis Code
source("Scripts/Misc/palette.R")

save.figure <- function(plot, name, path, pdf=TRUE, png=TRUE){
  #'@param plot gg object
  #'@param name string:  on the format '<name>.png'
  #'@param path string: path to location where figure should be stored
  #'@param pdf (boolean): whether to save as pdf
  #'@param png (boolean) whether to save as png
  
  if (png){
    ggsave(paste(name, '.png', sep=""),
           plot = plot,
           device = "png",
           path = path,
           height = 5, width = 8,
           dpi = "retina"
    )
  }
  if (pdf) {
    ggsave(paste(name, '.pdf', sep=""),
           plot = plot,
           device = "pdf",
           path = path,
           height = 5, width = 8,
           dpi = "retina"
    )
  }
}

trace_plot <- function(draws, iterations, warmup, chains, title){
  #' produces trace plots from array of draws from STAN run
  #' 
  #' @param draws (array<float>): HMC draws of parameter by STAN
  #' @param iterations (int): number of iterations for each chain in STAN analysis
  #' @param warmup (int): number of warm-up runs in STAN analysis
  #' @param chains (int): number of chains in STAN analysis
  #' @param title (string): title of plot
  
  draws.df <- data.frame(draw = draws) %>%
    mutate(no = rep(1:(iterations - warmup), chains)) %>%
    mutate(chain = as.character(rep(1:chains, each = iterations - warmup)))
  
  trace_plot <- ggplot(data = draws.df) + 
    geom_line(aes(x=no, y = draw, color=chain), alpha = 0.7) + 
    scale_color_manual(values=palette) + 
    labs(title=title, x = "iterations", y = " ")
    
  return(trace_plot)
}

trace_plot_matrix <- function(draws_mat, iterations, warmup, chains, title){
  draws.df <- data.frame(draws_mat) %>%
    mutate(no = rep(1:(iterations-warmup), chains)) %>%
    mutate(chain = as.character(rep(1:chains, each = iterations - warmup))) %>% 
    pivot_longer(cols = !c(no, chain), names_to = "parameter",
               values_to = "draw", names_prefix = "X",
               names_transform = list("parameter" = as.integer))
  
  trace_plot <- ggplot(data = draws.df) +
    geom_line(aes(x=no, y = draw, color = chain), alpha = 0.7) + 
    scale_color_manual(values = palette) +
    facet_wrap(~parameter) + 
    labs(title=title, x = "iterations", y = " ")
  
  return(trace_plot)
  
}

plot.counts.inlabru.stan.compared <- function(comparison.Y, path.to.storage, pdf=T, png=F){
  #' Produces and saves plots of counts
  #' 
  #' @param comparison.Y (data.frame): containts counts estimated by stan and inlabru, as well as observed counts. 
  #' @param path.to.storage (string): path to loaction where plots should be stored
  #' @param pdf (boolean): whether to save figures in pdf format
  #' @param png (boolean): whether to save figures in png format
  #' 
  p.counts.xt <- ggplot(data = comparison.Y, aes(x = xt)) + 
    geom_point(aes(y = mean.inlabru, color = "Inlabru", fill = "Inlabru")) + 
    geom_ribbon(aes(ymin = X0.025.inlabru, ymax = X0.975.inlabru, fill = "Inlabru"), alpha = 0.5) +
    geom_point(aes(y = mean.stan, color = "Stan", fill = "Stan")) +
    geom_ribbon(aes(ymin = X0.025.stan, ymax = X0.975.stan, fill = "Stan"), alpha = 0.5) +
    geom_point(aes(y = Y.observed, color = "Observed", fill = "Observed"), shape = 4) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    scale_shape_manual(name = "") +
    theme_classic() + 
    labs(title = "Estimated cancer death counts", x = "t, x", y = " ")
  
  save.figure(p.counts.xt, "counts_xt_comparison", path = path.to.storage, pdf = pdf, png = png)
  
  p.counts.t <- ggplot(data = comparison.Y, aes(x = x)) + 
    geom_point(aes(y = mean.inlabru, color = "Inlabru", fill = "Inlabru")) + 
    geom_ribbon(aes(ymin = X0.025.inlabru, ymax = X0.975.inlabru, fill = "Inlabru"), alpha = 0.5) +
    geom_point(aes(y = mean.stan, color = "Stan", fill = "Stan")) +
    geom_ribbon(aes(ymin = X0.025.stan, ymax = X0.975.stan, fill = "Stan"), alpha = 0.5) +
    geom_point(aes(y = Y.observed, color = "Observed", fill = "Observed"), shape = 4) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    scale_shape_manual(name = "") +
    theme_classic() + 
    labs(title = "Estimated cancer death counts", x = "t, x", y = " ") + 
    facet_wrap(~t)
  
  save.figure(p.counts.t, "counts_t_comparison", path = path.to.storage, pdf = pdf, png = png)
  
  p.counts.x <- ggplot(data = comparison.Y, aes(x = t)) + 
    geom_point(aes(y = mean.inlabru, color = "Inlabru", fill = "Inlabru")) + 
    geom_ribbon(aes(ymin = X0.025.inlabru, ymax = X0.975.inlabru, fill = "Inlabru"), alpha = 0.5) +
    geom_point(aes(y = mean.stan, color = "Stan", fill = "Stan")) +
    geom_ribbon(aes(ymin = X0.025.stan, ymax = X0.975.stan, fill = "Stan"), alpha = 0.5) +
    geom_point(aes(y = Y.observed, color = "Observed", fill = "Observed"), shape = 4) +
    scale_color_manual(name = "", values = palette ) +
    scale_fill_manual(name = "", values = palette ) +
    scale_shape_manual(name = "") +
    theme_classic() + 
    labs(title = "Estimated cancer death counts", x = "t, x", y = " ") + 
    facet_wrap(~x)
  
  save.figure(p.counts.x, "counts_x_comparison", path = path.to.storage, pdf = pdf, png = png)
}

# plot.predictor.inlabru.stan.compared.old <- function(
#   inlabru.predictor.df, stan.predictor.df,
#   path.to.storage, columns, pdf=T, png=F, a45=F) {
#   
#   if(!a45){
#     p.predictor.1 <- ggplot() + 
#       geom_density(data = inlabru.predictor.df, aes(x = X1,   color = "Inlabru", fill = "Inlabru"), alpha = 0.5, bins = 50) + 
#       geom_density(data = stan.predictor.df, aes(x = X1,   color = "Stan", fill = "Stan"), alpha = 0.5, bins=50) + 
#       scale_color_manual(name = " ", values = palette) + 
#       scale_fill_manual(name = " ", values = palette) +
#       theme_classic() + 
#       labs(x = "Value of predictor", y = " ", title = "Marginal of eta[1]")
#     
#     p.predictor.2 <- ggplot(data = data.frame(), aes(x = X64, y = after_stat(density))) + 
#       geom_density(data = inlabru.predictor.df, aes(color = "Inlabru", fill = "Inlabru"), alpha = 0.5, bins = 50) + 
#       geom_density(data = stan.predictor.df, aes(color = "Stan", fill = "Stan"), alpha = 0.5, bins = 50) + 
#       scale_color_manual(name = " ", values = palette) + 
#       scale_fill_manual(name = " ", values = palette) +
#       theme_classic() + 
#       labs(x = "Value of predictor", y = " ", title = "Marginal of eta[64]")
#     
#     p.predictor.3 <- ggplot(data = data.frame(), aes(x = X128, y = after_stat(density))) + 
#       geom_density(data = inlabru.predictor.df, aes(color = "Inlabru", fill = "Inlabru"), alpha = 0.5, bins = 50) + 
#       geom_density(data = stan.predictor.df, aes(color = "Stan", fill = "Stan"), alpha = 0.5, bins = 50) + 
#       scale_color_manual(name = " ", values = palette) + 
#       scale_fill_manual(name = " ", values = palette) +
#       theme_classic() + 
#       labs(x = "Value of predictor", y = " ", title = "Marginal of eta[128]")
#     
#     p.predictor.4 <- ggplot(data = data.frame(), aes(x = X192, y = after_stat(density))) + 
#       geom_density(data = inlabru.predictor.df, aes(color = "Inlabru", fill = "Inlabru"), alpha = 0.5, bins = 50) + 
#       geom_density(data = stan.predictor.df, aes(color = "Stan", fill = "Stan"), alpha = 0.5, bins = 50) + 
#       scale_color_manual(name = " ", values = palette) + 
#       scale_fill_manual(name = " ", values = palette) +
#       theme_classic() + 
#       labs(x = "Value of predictor", y = " ", title = "Marginal of eta[192]")
#     
#     p.predictor.5 <- ggplot(data = data.frame(), aes(x = X256, y = after_stat(density))) + 
#       geom_density(data = inlabru.predictor.df, aes(color = "Inlabru", fill = "Inlabru"), alpha = 0.5, bins = 50) + 
#       geom_density(data = stan.predictor.df, aes(color = "Stan", fill = "Stan"), alpha = 0.5, bins = 50) + 
#       scale_color_manual(name = " ", values = palette) + 
#       scale_fill_manual(name = " ", values = palette) +
#       theme_classic() + 
#       labs(x = "Value of predictor", y = " ", title = "Marginal of eta[256]")
#     
#     p.predictor.6 <- ggplot(data = data.frame(), aes(x = X324, y = after_stat(density))) + 
#       geom_density(data = inlabru.predictor.df, aes(color = "Inlabru", fill = "Inlabru"), alpha = 0.5, bins = 50) + 
#       geom_density(data = stan.predictor.df, aes(color = "Stan", fill = "Stan"), alpha = 0.5, bins = 50) + 
#       scale_color_manual(name = " ", values = palette) + 
#       scale_fill_manual(name = " ", values = palette) +
#       theme_classic() + 
#       labs(x = "Value of predictor", y = " ", title = "Marginal of eta[324]")
#     
#     p.predictor <- (p.predictor.1 | p.predictor.2 | p.predictor.3)/(p.predictor.4 | p.predictor.5 | p.predictor.6) + 
#       plot_layout(guides = "collect")
#     
#     save.figure(p.predictor, "predictor_marginals_comparison", path = path.to.storage, png = png, pdf = pdf)
#   }
#   if(a45){
#     
#     p.predictor.1 <- ggplot() + 
#       geom_density(data = inlabru.predictor.df, aes(x = X1,   color = "Inlabru", fill = "Inlabru"), alpha = 0.5, bins = 50) + 
#       geom_density(data = stan.predictor.df, aes(x = X1,   color = "Stan", fill = "Stan"), alpha = 0.5, bins=50) + 
#       scale_color_manual(name = " ", values = palette) + 
#       scale_fill_manual(name = " ", values = palette) +
#       theme_classic() + 
#       labs(x = "Value of predictor", y = " ", title = "Marginal of eta[1]")
#     
#     p.predictor.2 <- ggplot(data = data.frame(), aes(x = X32, y = after_stat(density))) + 
#       geom_density(data = inlabru.predictor.df, aes(color = "Inlabru", fill = "Inlabru"), alpha = 0.5, bins = 50) + 
#       geom_density(data = stan.predictor.df, aes(color = "Stan", fill = "Stan"), alpha = 0.5, bins = 50) + 
#       scale_color_manual(name = " ", values = palette) + 
#       scale_fill_manual(name = " ", values = palette) +
#       theme_classic() + 
#       labs(x = "Value of predictor", y = " ", title = "Marginal of eta[32]")
#     
#     p.predictor.3 <- ggplot(data = data.frame(), aes(x = X64, y = after_stat(density))) + 
#       geom_density(data = inlabru.predictor.df, aes(color = "Inlabru", fill = "Inlabru"), alpha = 0.5, bins = 50) + 
#       geom_density(data = stan.predictor.df, aes(color = "Stan", fill = "Stan"), alpha = 0.5, bins = 50) + 
#       scale_color_manual(name = " ", values = palette) + 
#       scale_fill_manual(name = " ", values = palette) +
#       theme_classic() + 
#       labs(x = "Value of predictor", y = " ", title = "Marginal of eta[64]")
#     
#     p.predictor.4 <- ggplot(data = data.frame(), aes(x = X96, y = after_stat(density))) + 
#       geom_density(data = inlabru.predictor.df, aes(color = "Inlabru", fill = "Inlabru"), alpha = 0.5, bins = 50) + 
#       geom_density(data = stan.predictor.df, aes(color = "Stan", fill = "Stan"), alpha = 0.5, bins = 50) + 
#       scale_color_manual(name = " ", values = palette) + 
#       scale_fill_manual(name = " ", values = palette) +
#       theme_classic() + 
#       labs(x = "Value of predictor", y = " ", title = "Marginal of eta[96]")
#     
#     p.predictor.5 <- ggplot(data = data.frame(), aes(x = X128, y = after_stat(density))) + 
#       geom_density(data = inlabru.predictor.df, aes(color = "Inlabru", fill = "Inlabru"), alpha = 0.5, bins = 50) + 
#       geom_density(data = stan.predictor.df, aes(color = "Stan", fill = "Stan"), alpha = 0.5, bins = 50) + 
#       scale_color_manual(name = " ", values = palette) + 
#       scale_fill_manual(name = " ", values = palette) +
#       theme_classic() + 
#       labs(x = "Value of predictor", y = " ", title = "Marginal of eta[128]")
#     
#     p.predictor.6 <- ggplot(data = data.frame(), aes(x = X162, y = after_stat(density))) + 
#       geom_density(data = inlabru.predictor.df, aes(color = "Inlabru", fill = "Inlabru"), alpha = 0.5, bins = 50) + 
#       geom_density(data = stan.predictor.df, aes(color = "Stan", fill = "Stan"), alpha = 0.5, bins = 50) + 
#       scale_color_manual(name = " ", values = palette) + 
#       scale_fill_manual(name = " ", values = palette) +
#       theme_classic() + 
#       labs(x = "Value of predictor", y = " ", title = "Marginal of eta[162]")
#     
#     p.predictor <- (p.predictor.1 | p.predictor.2 | p.predictor.3)/(p.predictor.4 | p.predictor.5 | p.predictor.6) + 
#       plot_layout(guides = "collect")
#     
#     save.figure(p.predictor, "predictor_marginals_comparison", path = path.to.storage, png = png, pdf = pdf)
#   }
# }

plot.predictor.inlabru.stan.compared <- function(
  res.inlabru, stan.predictor.df,
  path.to.storage, columns, pdf=T, png=F, a45=F, filename = "predictor_marginals_comparison") {
  
  if(!a45){
    pred.1.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.001)
  
    p.predictor.1 <- ggplot() + 
      geom_area(data = pred.1.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X1, color = "Stan", fill = "Stan"), alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Predictor[1]", x = " ", y = " ")
    
    pred.2.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.064)
    
    p.predictor.2 <- ggplot() + 
      geom_area(data = pred.2.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X64,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Predictor[64]", x = " ", y = " ")
    
    pred.3.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.128)
    
    p.predictor.3 <- ggplot() + 
      geom_area(data = pred.3.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X128,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Predictor[128]", x = " ", y = " ")
    
    pred.4.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.192)
    
    p.predictor.4 <- ggplot() + 
      geom_area(data = pred.4.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X192,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Predictor[192]", x = " ", y = " ")
    
    pred.5.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.256)
    
    p.predictor.5 <- ggplot() + 
      geom_area(data = pred.5.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X256,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Predictor[256]", x = " ", y = " ")
    
    pred.6.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.324)
    
    p.predictor.6 <- ggplot() + 
      geom_area(data = pred.6.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X324,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Predictor[324]", x = " ", y = " ")
    
    p.predictor <- (p.predictor.1 | p.predictor.2 | p.predictor.3)/(p.predictor.4 | p.predictor.5 | p.predictor.6) + 
      plot_layout(guides = "collect")
    
    save.figure(p.predictor, filename, path = path.to.storage, png = png, pdf = pdf)
  }
  if(a45){
    
    pred.1.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.001)
    
    p.predictor.1 <- ggplot() + 
      geom_area(data = pred.1.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X1,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Predictor[1]", x = " ", y = " ")
    
    pred.2.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.032)
    
    p.predictor.2 <- ggplot() + 
      geom_area(data = pred.2.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X32,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Predictor[32]", x = " ", y = " ")
    
    pred.3.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.064)
    
    p.predictor.3 <- ggplot() + 
      geom_area(data = pred.3.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X64,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Predictor[64]", x = " ", y = " ")
    
    pred.4.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.096)
    
    p.predictor.4 <- ggplot() + 
      geom_area(data = pred.4.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X96,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Predictor[96]", x = " ", y = " ")
    
    pred.5.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.128)
    
    p.predictor.5 <- ggplot() + 
      geom_area(data = pred.5.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X128,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Predictor[128]", x = " ", y = " ")
    
    pred.6.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.162)
    
    p.predictor.6 <- ggplot() + 
      geom_area(data = pred.6.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X162,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Predictor[162]", x = " ", y = " ")
    
    p.predictor <- (p.predictor.1 | p.predictor.2 | p.predictor.3)/(p.predictor.4 | p.predictor.5 | p.predictor.6) + 
      plot_layout(guides = "collect")
    
    save.figure(p.predictor, filename, path = path.to.storage, png = png, pdf = pdf)
  }
}

plot.alpha.inlabru.stan.compared <- function(
  res.inlabru, stan.predictor.df,
  path.to.storage, columns, pdf=T, png=F, a45=F, filename = "alpha_marginals_comparison") {
  
  if(!a45){
    pred.1.inlabru <- data.frame(res.inlabru$marginals.random$alpha$index.1)
    
    p.alpha.1 <- ggplot() + 
      geom_area(data = pred.1.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X1,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Alpha[1]", x = " ", y = " ")
    
    pred.3.inlabru <- data.frame(res.inlabru$marginals.random$alpha$index.3)
    
    p.alpha.3 <- ggplot() + 
      geom_area(data = pred.3.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X3,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Alpha[3]", x = " ", y = " ")
    
    pred.5.inlabru <- data.frame(res.inlabru$marginals.random$alpha$index.5)
    
    p.alpha.5 <- ggplot() + 
      geom_area(data = pred.5.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X5,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Alpha[5]", x = " ", y = " ")
    
    pred.7.inlabru <- data.frame(res.inlabru$marginals.random$alpha$index.7)
    
    p.alpha.7 <- ggplot() + 
      geom_area(data = pred.7.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X7,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Alpha[7]", x = " ", y = " ")
    
    pred.9.inlabru <- data.frame(res.inlabru$marginals.random$alpha$index.9)
    
    p.alpha.9 <- ggplot() + 
      geom_area(data = pred.9.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X9,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Alpha[9]", x = " ", y = " ")
    
    pred.11.inlabru <- data.frame(res.inlabru$marginals.random$alpha$index.11)
    
    p.alpha.11 <- ggplot() + 
      geom_area(data = pred.11.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X11,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Alpha[11]", x = " ", y = " ")
    
    pred.13.inlabru <- data.frame(res.inlabru$marginals.random$alpha$index.13)
    
    p.alpha.13 <- ggplot() + 
      geom_area(data = pred.13.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X13,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Alpha[13]", x = " ", y = " ")
    
    pred.15.inlabru <- data.frame(res.inlabru$marginals.random$alpha$index.15)
    
    p.alpha.15 <- ggplot() + 
      geom_area(data = pred.15.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X15,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Alpha[15]", x = " ", y = " ")
    
    pred.17.inlabru <- data.frame(res.inlabru$marginals.random$alpha$index.17)
    
    p.alpha.17 <- ggplot() + 
      geom_area(data = pred.17.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X17,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Alpha[17]", x = " ", y = " ")
    
    p.alpha <- (p.alpha.1 | p.alpha.3 | p.alpha.5) / (p.alpha.7 | p.alpha.9 | p.alpha.11) / (p.alpha.13 | p.alpha.15 | p.alpha.17) + 
      plot_layout(guides = "collect")
    
    save.figure(p.alpha, filename, path = path.to.storage, png = png, pdf = pdf)
  }
  if(a45){
    
    pred.1.inlabru <- data.frame(res.inlabru$marginals.random$alpha$index.1)
    
    p.alpha.1 <- ggplot() + 
      geom_area(data = pred.1.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X1,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Alpha[1]", x = " ", y = " ")
    
    pred.2.inlabru <- data.frame(res.inlabru$marginals.random$alpha$index.2)
    
    p.alpha.2 <- ggplot() + 
      geom_area(data = pred.2.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X2,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Alpha[2]", x = " ", y = " ")
    
    pred.3.inlabru <- data.frame(res.inlabru$marginals.random$alpha$index.3)
    
    p.alpha.3 <- ggplot() + 
      geom_area(data = pred.3.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X3,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Alpha[3]", x = " ", y = " ")
    
    pred.4.inlabru <- data.frame(res.inlabru$marginals.random$alpha$index.4)
    
    p.alpha.4 <- ggplot() + 
      geom_area(data = pred.4.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X4,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Alpha[4]", x = " ", y = " ")
    
    pred.5.inlabru <- data.frame(res.inlabru$marginals.random$alpha$index.5)
    
    p.alpha.5 <- ggplot() + 
      geom_area(data = pred.5.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X5,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Alpha[5]", x = " ", y = " ")
    
    pred.6.inlabru <- data.frame(res.inlabru$marginals.random$alpha$index.6)
    
    p.alpha.6 <- ggplot() + 
      geom_area(data = pred.6.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X6,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Alpha[6]", x = " ", y = " ")
    
    pred.7.inlabru <- data.frame(res.inlabru$marginals.random$alpha$index.7)
    
    p.alpha.7 <- ggplot() + 
      geom_area(data = pred.7.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X7,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Alpha[7]", x = " ", y = " ")
    
    pred.8.inlabru <- data.frame(res.inlabru$marginals.random$alpha$index.8)
    
    p.alpha.8 <- ggplot() + 
      geom_area(data = pred.8.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X8,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Alpha[8]", x = " ", y = " ")
    
    
    pred.9.inlabru <- data.frame(res.inlabru$marginals.random$alpha$index.9)
    
    p.alpha.9 <- ggplot() + 
      geom_area(data = pred.9.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X9,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Alpha[9]", x = " ", y = " ")
    
    
    p.alpha <- (p.alpha.1 | p.alpha.2 | p.alpha.3) / (p.alpha.4 | p.alpha.5 | p.alpha.6) / (p.alpha.7 | p.alpha.8 | p.alpha.9) + 
      plot_layout(guides = "collect")
    
    save.figure(p.alpha, filename, path = path.to.storage, png = png, pdf = pdf)
  }
}

plot.beta.inlabru.stan.compared <- function(
  res.inlabru, stan.predictor.df,
  path.to.storage, columns, pdf=T, png=F, a45=F, filename = "beta_marginals_comparison") {
  
  if(!a45){
    pred.1.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.1)
    
    p.beta.1 <- ggplot() + 
      geom_area(data = pred.1.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X1,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Beta[1]", x = " ", y = " ")
    
    pred.3.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.3)
    
    p.beta.3 <- ggplot() + 
      geom_area(data = pred.3.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X3,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Beta[3]", x = " ", y = " ")
    
    pred.5.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.5)
    
    p.beta.5 <- ggplot() + 
      geom_area(data = pred.5.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X5,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Beta[5]", x = " ", y = " ")
    
    pred.7.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.7)
    
    p.beta.7 <- ggplot() + 
      geom_area(data = pred.7.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X7,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Beta[7]", x = " ", y = " ")
    
    pred.9.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.9)
    
    p.beta.9 <- ggplot() + 
      geom_area(data = pred.9.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X9,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Beta[9]", x = " ", y = " ")
    
    pred.11.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.11)
    
    p.beta.11 <- ggplot() + 
      geom_area(data = pred.11.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X11,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Beta[11]", x = " ", y = " ")
    
    pred.13.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.13)
    
    p.beta.13 <- ggplot() + 
      geom_area(data = pred.13.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X13,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Beta[13]", x = " ", y = " ")
    
    pred.15.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.15)
    
    p.beta.15 <- ggplot() + 
      geom_area(data = pred.15.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X15,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Beta[15]", x = " ", y = " ")
    
    pred.17.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.17)
    
    p.beta.17 <- ggplot() + 
      geom_area(data = pred.17.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X17,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Beta[17]", x = " ", y = " ")
    
    p.beta <- (p.beta.1 | p.beta.3 | p.beta.5) / (p.beta.7 | p.beta.9 | p.beta.11) / (p.beta.13 | p.beta.15 | p.beta.17) + 
      plot_layout(guides = "collect")
    
    save.figure(p.beta, filename, path = path.to.storage, png = png, pdf = pdf)
  }
  if(a45){
    
    pred.1.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.1)
    
    p.beta.1 <- ggplot() + 
      geom_area(data = pred.1.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X1,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Beta[1]", x = " ", y = " ")
    
    pred.2.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.2)
    
    p.beta.2 <- ggplot() + 
      geom_area(data = pred.2.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X2,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Beta[2]", x = " ", y = " ")
    
    pred.3.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.3)
    
    p.beta.3 <- ggplot() + 
      geom_area(data = pred.3.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X3,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Beta[3]", x = " ", y = " ")
    
    pred.4.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.4)
    
    p.beta.4 <- ggplot() + 
      geom_area(data = pred.4.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X4,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Beta[4]", x = " ", y = " ")
    
    pred.5.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.5)
    
    p.beta.5 <- ggplot() + 
      geom_area(data = pred.5.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X5,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Beta[5]", x = " ", y = " ")
    
    pred.6.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.6)
    
    p.beta.6 <- ggplot() + 
      geom_area(data = pred.6.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X6,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Beta[6]", x = " ", y = " ")
    
    pred.7.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.7)
    
    p.beta.7 <- ggplot() + 
      geom_area(data = pred.7.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X7,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Beta[7]", x = " ", y = " ")
    
    pred.8.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.8)
    
    p.beta.8 <- ggplot() + 
      geom_area(data = pred.8.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X8,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Beta[8]", x = " ", y = " ")

    
    pred.9.inlabru <- data.frame(res.inlabru$marginals.random$beta$index.9)
    
    p.beta.9 <- ggplot() + 
      geom_area(data = pred.9.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.predictor.df, aes(x = X9,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Beta[9]", x = " ", y = " ")
    
    
    p.beta <- (p.beta.1 | p.beta.2 | p.beta.3) / (p.beta.4 | p.beta.5 | p.beta.6) / (p.beta.7 | p.beta.8 | p.beta.9) + 
      plot_layout(guides = "collect")
    
    save.figure(p.beta, filename, path = path.to.storage, png = png, pdf = pdf)
  }
}

plot.kappa.inlabru.stan.compared <- function(
  res.inlabru, stan.predictor.df,
  path.to.storage, columns, pdf=T, png=F, filename = "kappa_marginals_comparison") {
  
  pred.1.inlabru <- data.frame(res.inlabru$marginals.random$kappa$index.1)
  
  p.kappa.1 <- ggplot() + 
    geom_area(data = pred.1.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_density(data = stan.predictor.df, aes(x = X1,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Kappa[1]", x = " ", y = " ")
  
  pred.3.inlabru <- data.frame(res.inlabru$marginals.random$kappa$index.3)
  
  p.kappa.3 <- ggplot() + 
    geom_area(data = pred.3.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_density(data = stan.predictor.df, aes(x = X3,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Kappa[3]", x = " ", y = " ")
  
  pred.5.inlabru <- data.frame(res.inlabru$marginals.random$kappa$index.5)
  
  p.kappa.5 <- ggplot() + 
    geom_area(data = pred.5.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_density(data = stan.predictor.df, aes(x = X5,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Kappa[5]", x = " ", y = " ")
  
  pred.7.inlabru <- data.frame(res.inlabru$marginals.random$kappa$index.7)
  
  p.kappa.7 <- ggplot() + 
    geom_area(data = pred.7.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_density(data = stan.predictor.df, aes(x = X7,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Kappa[7]", x = " ", y = " ")
  
  
  pred.9.inlabru <- data.frame(res.inlabru$marginals.random$kappa$index.9)
  
  p.kappa.9 <- ggplot() + 
    geom_area(data = pred.9.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_density(data = stan.predictor.df, aes(x = X9,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Kappa[9]", x = " ", y = " ")
  
  pred.11.inlabru <- data.frame(res.inlabru$marginals.random$kappa$index.11)
  
  p.kappa.11 <- ggplot() + 
    geom_area(data = pred.11.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_density(data = stan.predictor.df, aes(x = X11,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Kappa[11]", x = " ", y = " ")
  
  pred.13.inlabru <- data.frame(res.inlabru$marginals.random$kappa$index.13)
  
  p.kappa.13 <- ggplot() + 
    geom_area(data = pred.13.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_density(data = stan.predictor.df, aes(x = X13,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Kappa[13]", x = " ", y = " ")
  
  pred.15.inlabru <- data.frame(res.inlabru$marginals.random$kappa$index.15)
  
  p.kappa.15 <- ggplot() + 
    geom_area(data = pred.15.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_density(data = stan.predictor.df, aes(x = X15,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Kappa[15]", x = " ", y = " ")
  
  pred.17.inlabru <- data.frame(res.inlabru$marginals.random$kappa$index.17)
  
  p.kappa.17 <- ggplot() + 
    geom_area(data = pred.17.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
    geom_density(data = stan.predictor.df, aes(x = X17,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
    theme_classic() + 
    scale_color_manual(name = "", values = palette) + 
    scale_fill_manual(name = "", values = palette) + 
    labs(title = "Kappa[17]", x = " ", y = " ")
  
  p.kappa <- (p.kappa.1 | p.kappa.3 | p.kappa.5) / (p.kappa.7 | p.kappa.9 | p.kappa.11) / (p.kappa.13 | p.kappa.15 | p.kappa.17) + 
    plot_layout(guides = "collect")
  
  save.figure(p.kappa, filename, path = path.to.storage, png = png, pdf = pdf)
}

plot.epsilon.inlabru.stan.compared <- function(
  res.inlabru, stan.epsilon.df,
  path.to.storage, pdf=T, png=F, a45=F, filename = "epsilon_marginals_comparison") {
  
  if(!a45){
    pred.1.inlabru <- data.frame(res.inlabru$marginals.random$epsilon$index.1)
    
    p.epsilon.1 <- ggplot() + 
      geom_area(data = pred.1.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.epsilon.df, aes(x = X1,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Epsilon[1]", x = " ", y = " ")
    
    pred.2.inlabru <- data.frame(res.inlabru$marginals.random$epsilon$index.64)
    
    p.epsilon.2 <- ggplot() + 
      geom_area(data = pred.2.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.epsilon.df, aes(x = X64,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Epsilon[64]", x = " ", y = " ")
    
    pred.3.inlabru <- data.frame(res.inlabru$marginals.random$epsilon$index.128)
    
    p.epsilon.3 <- ggplot() + 
      geom_area(data = pred.3.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.epsilon.df, aes(x = X128,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Epsilon[128]", x = " ", y = " ")
    
    pred.4.inlabru <- data.frame(res.inlabru$marginals.random$epsilon$index.192)
    
    p.epsilon.4 <- ggplot() + 
      geom_area(data = pred.4.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.epsilon.df, aes(x = X192,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Epsilon[192]", x = " ", y = " ")
    
    pred.5.inlabru <- data.frame(res.inlabru$marginals.random$epsilon$index.256)
    
    p.epsilon.5 <- ggplot() + 
      geom_area(data = pred.5.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.epsilon.df, aes(x = X256,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Epsilon[256]", x = " ", y = " ")
    
    pred.6.inlabru <- data.frame(res.inlabru$marginals.random$epsilon$index.324)
    
    p.epsilon.6 <- ggplot() + 
      geom_area(data = pred.6.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.epsilon.df, aes(x = X324,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Epsilon[324]", x = " ", y = " ")
    
    p.epsilon <- (p.epsilon.1 | p.epsilon.2 | p.epsilon.3)/(p.epsilon.4 | p.epsilon.5 | p.epsilon.6) + 
      plot_layout(guides = "collect")
    
    save.figure(p.epsilon, filename, path = path.to.storage, png = png, pdf = pdf)
  }
  if(a45){
    
    pred.1.inlabru <- data.frame(res.inlabru$marginals.random$epsilon$index.1)
    
    p.epsilon.1 <- ggplot() + 
      geom_area(data = pred.1.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.epsilon.df, aes(x = X1,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Epsilon[1]", x = " ", y = " ")
    
    pred.2.inlabru <- data.frame(res.inlabru$marginals.random$epsilon$index.32)
    
    p.epsilon.2 <- ggplot() + 
      geom_area(data = pred.2.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.epsilon.df, aes(x = X32,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Epsilon[32]", x = " ", y = " ")
    
    pred.3.inlabru <- data.frame(res.inlabru$marginals.random$epsilon$index.64)
    
    p.epsilon.3 <- ggplot() + 
      geom_area(data = pred.3.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.epsilon.df, aes(x = X64,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Epsilon[64]", x = " ", y = " ")
    
    pred.4.inlabru <- data.frame(res.inlabru$marginals.random$epsilon$index.96)
    
    p.epsilon.4 <- ggplot() + 
      geom_area(data = pred.4.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.epsilon.df, aes(x = X96,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Epsilon[96]", x = " ", y = " ")
    
    pred.5.inlabru <- data.frame(res.inlabru$marginals.random$epsilon$index.128)
    
    p.epsilon.5 <- ggplot() + 
      geom_area(data = pred.5.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.epsilon.df, aes(x = X128,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Epsilon[128]", x = " ", y = " ")
    
    pred.6.inlabru <- data.frame(res.inlabru$marginals.random$epsilon$index.162)
    
    p.epsilon.6 <- ggplot() + 
      geom_area(data = pred.6.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_density(data = stan.epsilon.df, aes(x = X162,   color = "Stan", fill = "Stan"),  alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Epsilon[162]", x = " ", y = " ")
    
    p.epsilon <- (p.epsilon.1 | p.epsilon.2 | p.epsilon.3)/(p.epsilon.4 | p.epsilon.5 | p.epsilon.6) + 
      plot_layout(guides = "collect")
    
    save.figure(p.epsilon, filename, path = path.to.storage, png = png, pdf = pdf)
  }
}