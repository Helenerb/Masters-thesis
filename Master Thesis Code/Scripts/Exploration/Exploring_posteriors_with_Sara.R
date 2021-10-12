# Dumping of code from exploring finding posterio distributions together with Sara


samps = inla.posterior.sample(res.inlabru.lc.1, n = 1000)

phi.plus.kappa <- function(t_max){
  t = 0:t_max
  res = kappa + phi*t
  return(res)
}

phi.plus.kappa.99 <- function(){
  t = 0:99
  res = kappa + phi*t
  return(res)
}

posterior <- inla.posterior.sample.eval(fun = phi.plus.kappa.99, samples=samps)

posterior.df <- data.frame(t = 1:100,
                           mean = apply(posterior, 1, mean),
                           q1 = apply(posterior, 1, quantile, 0.025),
                           q2 = apply(posterior, 1, quantile, 0.975)) %>%
  mutate(kappa = underlying.effects.lc$kappa.true[t]) %>%
  mutate(phi.t = underlying.effects.lc$phi.true*(t-1)) %>%
  mutate(kappa.phi = kappa + phi.t)

gg.posterior <- ggplot(data = posterior.df) +
  geom_ribbon(aes(x = t, ymin = q1, ymax = q2, color = "Inlabru", fill = "Inlabru"), alpha = 0.5) + 
  geom_point(aes(x = t, y = mean, color = "Inlabru", fill = "Inlabru")) +
  geom_point(aes(x = t, y = kappa.phi, color = "True", fill = "True")) +
  scale_color_manual(name = " ", values = palette.basis) + 
  scale_fill_manual(name = " ", values = palette.basis) + 
  labs(title = "Phi*t + kappa", x = "t", y = "")

gg.posterior

phi.kappa.samp <- inla.posterior.sample.eval(samples=samps, c("phi", "kappa"))

samp.matrix <- matrix(rep(1:20, 100), byrow = FALSE, nrow= 20)
phi.matrix <- matrix(rep(phi.kappa.samp[1,], 20), nrow = 20, byrow=TRUE)
samp.matrix <- samp.matrix * phi.matrix + phi.kappa.samp[-1,]

data.frame(t = 1:20, mean = apply(samp.matrix, 1, mean),
           q1 = apply(samp.matrix, 1, quantile, 0.025),
           q2 = apply(samp.matrix, 1, quantile, 0.975)) %>%
  ggplot() + 
  geom_line(aes(x = t, y= mean)) + 
  geom_line(aes(x = t, y = q1), alpha = 0.5) +
  geom_line(aes(x = t, y = q2), alpha = 0.5)

samples[1]

ggplot(data = posterior.phi) + geom_point(aes(x = t, y = mean)) + geom_line(aes(x = t, y =q0.025)) + geom_line(aes(x = t, y = q0.975))
