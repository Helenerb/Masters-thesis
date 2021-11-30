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

plot.predictor.inlabru.stan.compared.old <- function(
  inlabru.predictor.df, stan.predictor.df,
  path.to.storage, columns, pdf=T, png=F, a45=F) {
  
  if(!a45){
    p.predictor.1 <- ggplot() + 
      geom_histogram(data = inlabru.predictor.df, aes(x = X1, y = after_stat(density), color = "Inlabru", fill = "Inlabru"), alpha = 0.5, bins = 50) + 
      geom_histogram(data = stan.predictor.df, aes(x = X1, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins=50) + 
      scale_color_manual(name = " ", values = palette) + 
      scale_fill_manual(name = " ", values = palette) +
      theme_classic() + 
      labs(x = "Value of predictor", y = " ", title = "Marginal of eta[1]")
    
    p.predictor.2 <- ggplot(data = data.frame(), aes(x = X64, y = after_stat(density))) + 
      geom_histogram(data = inlabru.predictor.df, aes(color = "Inlabru", fill = "Inlabru"), alpha = 0.5, bins = 50) + 
      geom_histogram(data = stan.predictor.df, aes(color = "Stan", fill = "Stan"), alpha = 0.5, bins = 50) + 
      scale_color_manual(name = " ", values = palette) + 
      scale_fill_manual(name = " ", values = palette) +
      theme_classic() + 
      labs(x = "Value of predictor", y = " ", title = "Marginal of eta[64]")
    
    p.predictor.3 <- ggplot(data = data.frame(), aes(x = X128, y = after_stat(density))) + 
      geom_histogram(data = inlabru.predictor.df, aes(color = "Inlabru", fill = "Inlabru"), alpha = 0.5, bins = 50) + 
      geom_histogram(data = stan.predictor.df, aes(color = "Stan", fill = "Stan"), alpha = 0.5, bins = 50) + 
      scale_color_manual(name = " ", values = palette) + 
      scale_fill_manual(name = " ", values = palette) +
      theme_classic() + 
      labs(x = "Value of predictor", y = " ", title = "Marginal of eta[128]")
    
    p.predictor.4 <- ggplot(data = data.frame(), aes(x = X192, y = after_stat(density))) + 
      geom_histogram(data = inlabru.predictor.df, aes(color = "Inlabru", fill = "Inlabru"), alpha = 0.5, bins = 50) + 
      geom_histogram(data = stan.predictor.df, aes(color = "Stan", fill = "Stan"), alpha = 0.5, bins = 50) + 
      scale_color_manual(name = " ", values = palette) + 
      scale_fill_manual(name = " ", values = palette) +
      theme_classic() + 
      labs(x = "Value of predictor", y = " ", title = "Marginal of eta[192]")
    
    p.predictor.5 <- ggplot(data = data.frame(), aes(x = X256, y = after_stat(density))) + 
      geom_histogram(data = inlabru.predictor.df, aes(color = "Inlabru", fill = "Inlabru"), alpha = 0.5, bins = 50) + 
      geom_histogram(data = stan.predictor.df, aes(color = "Stan", fill = "Stan"), alpha = 0.5, bins = 50) + 
      scale_color_manual(name = " ", values = palette) + 
      scale_fill_manual(name = " ", values = palette) +
      theme_classic() + 
      labs(x = "Value of predictor", y = " ", title = "Marginal of eta[256]")
    
    p.predictor.6 <- ggplot(data = data.frame(), aes(x = X324, y = after_stat(density))) + 
      geom_histogram(data = inlabru.predictor.df, aes(color = "Inlabru", fill = "Inlabru"), alpha = 0.5, bins = 50) + 
      geom_histogram(data = stan.predictor.df, aes(color = "Stan", fill = "Stan"), alpha = 0.5, bins = 50) + 
      scale_color_manual(name = " ", values = palette) + 
      scale_fill_manual(name = " ", values = palette) +
      theme_classic() + 
      labs(x = "Value of predictor", y = " ", title = "Marginal of eta[324]")
    
    p.predictor <- (p.predictor.1 | p.predictor.2 | p.predictor.3)/(p.predictor.4 | p.predictor.5 | p.predictor.6) + 
      plot_layout(guides = "collect")
    
    save.figure(p.predictor, "predictor_marginals_comparison", path = path.to.storage, png = png, pdf = pdf)
  }
  if(a45){
    
    p.predictor.1 <- ggplot() + 
      geom_histogram(data = inlabru.predictor.df, aes(x = X1, y = after_stat(density), color = "Inlabru", fill = "Inlabru"), alpha = 0.5, bins = 50) + 
      geom_histogram(data = stan.predictor.df, aes(x = X1, y = after_stat(density), color = "Stan", fill = "Stan"), alpha = 0.5, bins=50) + 
      scale_color_manual(name = " ", values = palette) + 
      scale_fill_manual(name = " ", values = palette) +
      theme_classic() + 
      labs(x = "Value of predictor", y = " ", title = "Marginal of eta[1]")
    
    p.predictor.2 <- ggplot(data = data.frame(), aes(x = X32, y = after_stat(density))) + 
      geom_histogram(data = inlabru.predictor.df, aes(color = "Inlabru", fill = "Inlabru"), alpha = 0.5, bins = 50) + 
      geom_histogram(data = stan.predictor.df, aes(color = "Stan", fill = "Stan"), alpha = 0.5, bins = 50) + 
      scale_color_manual(name = " ", values = palette) + 
      scale_fill_manual(name = " ", values = palette) +
      theme_classic() + 
      labs(x = "Value of predictor", y = " ", title = "Marginal of eta[32]")
    
    p.predictor.3 <- ggplot(data = data.frame(), aes(x = X64, y = after_stat(density))) + 
      geom_histogram(data = inlabru.predictor.df, aes(color = "Inlabru", fill = "Inlabru"), alpha = 0.5, bins = 50) + 
      geom_histogram(data = stan.predictor.df, aes(color = "Stan", fill = "Stan"), alpha = 0.5, bins = 50) + 
      scale_color_manual(name = " ", values = palette) + 
      scale_fill_manual(name = " ", values = palette) +
      theme_classic() + 
      labs(x = "Value of predictor", y = " ", title = "Marginal of eta[64]")
    
    p.predictor.4 <- ggplot(data = data.frame(), aes(x = X96, y = after_stat(density))) + 
      geom_histogram(data = inlabru.predictor.df, aes(color = "Inlabru", fill = "Inlabru"), alpha = 0.5, bins = 50) + 
      geom_histogram(data = stan.predictor.df, aes(color = "Stan", fill = "Stan"), alpha = 0.5, bins = 50) + 
      scale_color_manual(name = " ", values = palette) + 
      scale_fill_manual(name = " ", values = palette) +
      theme_classic() + 
      labs(x = "Value of predictor", y = " ", title = "Marginal of eta[96]")
    
    p.predictor.5 <- ggplot(data = data.frame(), aes(x = X128, y = after_stat(density))) + 
      geom_histogram(data = inlabru.predictor.df, aes(color = "Inlabru", fill = "Inlabru"), alpha = 0.5, bins = 50) + 
      geom_histogram(data = stan.predictor.df, aes(color = "Stan", fill = "Stan"), alpha = 0.5, bins = 50) + 
      scale_color_manual(name = " ", values = palette) + 
      scale_fill_manual(name = " ", values = palette) +
      theme_classic() + 
      labs(x = "Value of predictor", y = " ", title = "Marginal of eta[128]")
    
    p.predictor.6 <- ggplot(data = data.frame(), aes(x = X162, y = after_stat(density))) + 
      geom_histogram(data = inlabru.predictor.df, aes(color = "Inlabru", fill = "Inlabru"), alpha = 0.5, bins = 50) + 
      geom_histogram(data = stan.predictor.df, aes(color = "Stan", fill = "Stan"), alpha = 0.5, bins = 50) + 
      scale_color_manual(name = " ", values = palette) + 
      scale_fill_manual(name = " ", values = palette) +
      theme_classic() + 
      labs(x = "Value of predictor", y = " ", title = "Marginal of eta[162]")
    
    p.predictor <- (p.predictor.1 | p.predictor.2 | p.predictor.3)/(p.predictor.4 | p.predictor.5 | p.predictor.6) + 
      plot_layout(guides = "collect")
    
    save.figure(p.predictor, "predictor_marginals_comparison", path = path.to.storage, png = png, pdf = pdf)
  }
}

plot.predictor.inlabru.stan.compared <- function(
  res.inlabru, stan.predictor.df,
  path.to.storage, columns, pdf=T, png=F, a45=F) {
  
  if(!a45){
    pred.1.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.001)
  
    p.predictor.1 <- ggplot() + 
      geom_area(data = pred.1.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_histogram(data = stan.predictor.df, aes(x = X1, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Predictor[1], Inlabru and stan", x = " ", y = " ")
    
    pred.2.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.064)
    
    p.predictor.2 <- ggplot() + 
      geom_area(data = pred.2.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_histogram(data = stan.predictor.df, aes(x = X64, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Predictor[64], Inlabru and stan", x = " ", y = " ")
    
    pred.3.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.128)
    
    p.predictor.3 <- ggplot() + 
      geom_area(data = pred.3.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_histogram(data = stan.predictor.df, aes(x = X128, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Predictor[128], Inlabru and stan", x = " ", y = " ")
    
    pred.4.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.192)
    
    p.predictor.4 <- ggplot() + 
      geom_area(data = pred.4.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_histogram(data = stan.predictor.df, aes(x = X192, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Predictor[192], Inlabru and stan", x = " ", y = " ")
    
    pred.5.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.256)
    
    p.predictor.5 <- ggplot() + 
      geom_area(data = pred.5.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_histogram(data = stan.predictor.df, aes(x = X256, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Predictor[256], Inlabru and stan", x = " ", y = " ")
    
    pred.6.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.324)
    
    p.predictor.6 <- ggplot() + 
      geom_area(data = pred.6.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_histogram(data = stan.predictor.df, aes(x = X324, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Predictor[324], Inlabru and stan", x = " ", y = " ")
    
    p.predictor <- (p.predictor.1 | p.predictor.2 | p.predictor.3)/(p.predictor.4 | p.predictor.5 | p.predictor.6) + 
      plot_layout(guides = "collect")
    
    save.figure(p.predictor, "predictor_marginals_comparison", path = path.to.storage, png = png, pdf = pdf)
  }
  if(a45){
    
    pred.1.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.001)
    
    p.predictor.1 <- ggplot() + 
      geom_area(data = pred.1.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_histogram(data = stan.predictor.df, aes(x = X1, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Predictor[1], Inlabru and stan", x = " ", y = " ")
    
    pred.2.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.032)
    
    p.predictor.2 <- ggplot() + 
      geom_area(data = pred.2.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_histogram(data = stan.predictor.df, aes(x = X32, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Predictor[32], Inlabru and stan", x = " ", y = " ")
    
    pred.3.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.064)
    
    p.predictor.3 <- ggplot() + 
      geom_area(data = pred.3.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_histogram(data = stan.predictor.df, aes(x = X64, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Predictor[64], Inlabru and stan", x = " ", y = " ")
    
    pred.4.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.096)
    
    p.predictor.4 <- ggplot() + 
      geom_area(data = pred.4.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_histogram(data = stan.predictor.df, aes(x = X96, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Predictor[96], Inlabru and stan", x = " ", y = " ")
    
    pred.5.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.128)
    
    p.predictor.5 <- ggplot() + 
      geom_area(data = pred.5.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_histogram(data = stan.predictor.df, aes(x = X128, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Predictor[128], Inlabru and stan", x = " ", y = " ")
    
    pred.6.inlabru <- data.frame(res.inlabru$marginals.linear.predictor$APredictor.162)
    
    p.predictor.6 <- ggplot() + 
      geom_area(data = pred.6.inlabru, aes(x = x, y = y, fill = "Inlabru", color = "Inlabru"), alpha = 0.5) + 
      geom_histogram(data = stan.predictor.df, aes(x = X162, y = after_stat(density), color = "Stan", fill = "Stan"), bins = 50, alpha = 0.5) + 
      theme_classic() + 
      scale_color_manual(name = "", values = palette) + 
      scale_fill_manual(name = "", values = palette) + 
      labs(title = "Predictor[162], Inlabru and stan", x = " ", y = " ")
    
    p.predictor <- (p.predictor.1 | p.predictor.2 | p.predictor.3)/(p.predictor.4 | p.predictor.5 | p.predictor.6) + 
      plot_layout(guides = "collect")
    
    save.figure(p.predictor, "predictor_marginals_comparison", path = path.to.storage, png = png, pdf = pdf)
  }
}