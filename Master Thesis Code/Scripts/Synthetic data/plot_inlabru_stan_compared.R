# script contatining functions for comparing inlabru and stan results

plot.inlabru.stan.compared.lc <- function(stan.summaries, inlabru.summaries, underlying.effects){
  obs <- underlying.effects$obs
  
  #  ----   intercept   ----
  p.intercept <- ggplot() + 
    geom_area(data=inlabru.summaries$data.fixed, aes(x = Int.x, y = Int.y, fill = "Inlabru"), alpha = 0.4) + 
    geom_vline(data=inlabru.summaries$data.fixed, aes(xintercept = mean[1], color = "Estimated", fill="Inlabru")) + 
    geom_vline(data=stan.summaries$summary.fixed, aes(xintercept = mean[2], color = "Stan")) + 
    geom_vline(aes(xintercept = `2.5%`[2], color = "Stan"), alpha = 0.5) +
    geom_vline(aes(xintercept = `97.5%`[2], color = "Stan"), alpha = 0.5) + 
    geom_vline(aes(xintercept = underlying.effects$age.intercept.true, color="True", fill="True")) +
    scale_color_manual(name = " ", values = palette) + 
    scale_fill_manual(name = " ", values = palette) +
    labs(x = "Value of intercept", y = " ", title = "Intercept")
  
  # ---   alpha   ----
  p.alpha <- ggplot() + 
    geom_ribbon(data=inlabru.summaries$data.alpha, aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, fill = "Inlabru"), alpha = 0.4) + 
    geom_point(data=inlabru.summaries$data.alpha, aes(x = ID, y = mean, color = "Inlabru", fill = "Inlabru"), size = 0.5) + 
    
    geom_point(data=stan.summaries$summary.alpha, aes(x=index, y=mean, color="Stan")) + 
    geom_line(data=stan.summaries$summary.alpha, aes(x=index, y=`2.5%`, color="Stan"), alpha=0.5) + 
    geom_line(data=stan.summaries$summary.alpha, aes(x=index, y=`97.5%`, color="Stan"), alpha=0.5) +
    
    geom_point(aes(x = ID, y = underlying.effects$alpha.true, color = "True", fill = "True"), size = 0.5) + 
    scale_color_manual(name = "",
                       values = palette ) +
    scale_fill_manual(name = "",
                      values = palette ) +
    labs(title="Alpha", x = "x", y='')
}

plot.inlabru.stan.compared.cohort <- function(stan.summaries, res.inlabru, underlying.effects){
  obs <- underlying.effects$obs
}
