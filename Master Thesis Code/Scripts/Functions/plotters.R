# script containing functions for plotting. 

save.figure <- function(plot, name, path){
  #'@param plot gg object
  #'@param name string:  on the format '<name>.png'
  #'@param path string: path to location where figure should be stored
  # ggsave(paste(name, '.png', sep=""),
  #        plot = plot,
  #        device = "png",
  #        path = path,
  #        height = 5, width = 8, 
  #        dpi = "retina"
  # )
  ggsave(paste(name, '.pdf', sep=""),
         plot = plot,
         device = "pdf",
         path = path,
         height = 5, width = 8, 
         dpi = "retina"
  )
}

confidence_bound_male_female <- function(param_name){
  plot <- ggplot() + 
    geom_ribbon(data = res.abKg$summary.random$alpha0, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = 'Male'), alpha = 0.4) +
    geom_ribbon(data = res.abKg$summary.random$alpha1, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = 'Female'), alpha = 0.4) +
    geom_point(data = res.abKg$summary.random$alpha0, aes(x = ID, y = mean, color = "Male")) + 
    geom_point(data = res.abKg$summary.random$alpha1, aes(x = ID, y = mean, color = "Female")) + 
    scale_color_manual(name = " ", values = palette.basis) + 
    scale_fill_manual(name = " ", values = palette.basis) +
    labs(x = "x", y = "alpha", title = "Alpha")
  
  return(plot)
}

