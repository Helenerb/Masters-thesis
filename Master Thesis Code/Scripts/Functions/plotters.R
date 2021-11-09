# script containing functions for plotting. 

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