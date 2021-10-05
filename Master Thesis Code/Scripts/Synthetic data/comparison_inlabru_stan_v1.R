source("configuration_v1.R")

# configuration without cohort effect
underlying.effects.lc <- configuration.v5()

obs.lc <- underlying.effects.lc$obs

# cohfiguration with cohort effect
underlying.effects.lc.cohort <- configuration.v7()

obs.lc.cohort <- underlying.effects.lc.cohort$obs

source("inlabru_analyses.R")

res.inlabru.lc.1 <- inlabru.lc.1(obs.lc)

res.inlabru.lc.cohort.1 <- inlabru.lc.cohort.1(obs.lc.cohort)

#   ----   plotting results of inlabru fit   ----
source("plot_inlabru_vs_underlying.R")

# plotting results from run with cohort effects:
plots.lc.cohort <- plot.inlabru.vs.underlying.v1(res.inlabru.lc.cohort.1,
                                                 underlying.effects.lc.cohort)
plots.lc.cohort$p.alpha
plots.lc.cohort$p.beta
plots.lc.cohort$p.phi
plots.lc.cohort$p.intercept
plots.lc.cohort$p.kappa
plots.lc.cohort$p.eta
plots.lc.cohort$p.eta.2
plots.lc.cohort$p.eta.t
plots.lc.cohort$p.eta.x
plots.lc.cohort$p.gamma

# plotting results from run with cohort effects:
plots <- plot.inlabru.vs.underlying.v5(res.inlabru.lc.1, underlying.effects.lc)
plots$p.alpha
plots$p.beta
plots$p.phi
plots$p.intercept
plots$p.kappa
plots$p.eta
plots$p.eta.2
plots$p.eta.t
plots$p.eta.x


#    ----   Old code: don´t know if you still need it yet.   ----

obs <- underlying.effects$obs

obs.pred = data.frame(x = obs$x, t = obs$t, xt = obs$xt, c = obs$c,
                      x.c = obs$x.c, t.c = obs$t.c, E = obs$E, Y = obs$Y,
                      age.intercept.obs = obs$age.intercept,
                      alpha.obs = obs$alpha, beta.obs = obs$beta,
                      phi.t.obs = obs$phi.t, kappa.obs = obs$kappa,
                      gamma.obs = obs$gamma, eta.obs = obs$eta,
                      epsilon.obs = obs$epsilon
                      )

# posterior distribution of beta*phi*t
#beta.phi.posterior <- predict(res.inlabru, data = NULL, formula = ~ beta*phi*t)
beta.phi.posterior <- predict(res.inlabru, data = obs.pred , formula = ~ beta*phi + beta*kappa)

ggplot(data = beta.phi.posterior) +
  geom_line(aes(x = xt, y = mean, color = "Predict")) + 
  geom_line(aes(x = xt, y = beta.obs*phi.t.obs + beta.obs*kappa.obs, color = "Observed"))

alpha.beta.kappa.posterior <- predict(res.inlabru, data = obs.pred , formula = ~ alpha + beta*phi + beta*kappa)

ggplot(data = alpha.beta.kappa.posterior) +
  geom_point(aes(x = xt, y = mean, color = "Predict"), size = 1) + 
  geom_point(aes(x = xt, y = alpha.obs + beta.obs*phi.t.obs + beta.obs*kappa.obs, color = "Observed"), size = 1)

alpha.beta.kappa.gamma.posterior <- predict(res.inlabru, data = obs.pred , formula = ~ alpha + beta*phi + beta*kappa + gamma)

ggplot(data = alpha.beta.kappa.gamma.posterior) +
  geom_point(aes(x = xt, y = mean, color = "Predict")) + 
  geom_point(aes(x = xt, y = alpha.obs + beta.obs*phi.t.obs + beta.obs*kappa.obs + gamma.obs, color = "Observed"))

# observe shioft for this configuration:
int.alpha.beta.kappa.gamma.posterior <- predict(res.inlabru, data = obs.pred , formula = ~ Int + alpha + beta*phi + beta*kappa + gamma)

ggplot(data = int.alpha.beta.kappa.gamma.posterior) +
  geom_point(aes(x = xt, y = mean, color = "Predict")) + 
  geom_point(aes(x = xt, y = age.intercept.obs + alpha.obs + beta.obs*phi.t.obs + beta.obs*kappa.obs + gamma.obs, color = "Observed"))

int.alpha.beta.kappa.gamma.epsilon.posterior <- predict(res.inlabru, data = obs.pred , formula = ~ Int + alpha + beta*phi + beta*kappa + gamma + epsilon)

ggplot(data = int.alpha.beta.kappa.gamma.epsilon.posterior) +
  geom_point(aes(x = xt, y = mean, color = "Predict")) + 
  geom_point(aes(x = xt, y = age.intercept.obs + alpha.obs + beta.obs*phi.t.obs + beta.obs*kappa.obs + gamma.obs  + epsilon.obs, color = "Observed"))

# attempt removing intercept from calculation
alpha.beta.kappa.gamma.epsilon.posterior <- predict(res.inlabru, data = obs.pred , formula = ~ alpha + beta*phi + beta*kappa + gamma + epsilon)

ggplot(data = alpha.beta.kappa.gamma.epsilon.posterior) +
  #geom_line(aes(x = xt, y = mean, color = "Predict wo Int")) + 
  geom_line(data = int.alpha.beta.kappa.gamma.epsilon.posterior,
             aes(x = xt, y = mean, color = "Predict w Int")) +
  # geom_line(aes(x = xt,
  #                y = age.intercept.obs + alpha.obs + beta.obs*phi.t.obs +
  #                  beta.obs*kappa.obs + gamma.obs  + epsilon.obs, color = "Observed")) +
  # geom_line(aes(x = xt,
  #                y = alpha.obs + beta.obs*phi.t.obs + beta.obs*kappa.obs +
  #                  gamma.obs  + epsilon.obs, color = "Observed wo age intercept")) +
  geom_line(aes(x = xt,
                 y = age.intercept.obs + alpha.obs + beta.obs*phi.t.obs +
                   beta.obs*kappa.obs + gamma.obs  + epsilon.obs -
                   mean(beta.obs*phi.t.obs),
                 color = "Observed mean(period effect extracted)")) +
  geom_line(aes(x = xt,
                y = age.intercept.obs + alpha.obs + beta.obs*phi.t.obs +
                  beta.obs*kappa.obs + gamma.obs  + epsilon.obs + 0.6296,
                color = "Observed mean(mean(period effect extracted))"))
  
all.posteroops <- predict(res.inlabru, data = obs.pred , formula = NULL)

ggplot() +
  geom_point(data = data.frame(mean = all.posteroops$Predictor$mean[1:8000]),
             aes(x = obs.pred$xt, y = mean, color = "Predict")) + 
  geom_point(data = obs.pred, aes(x = xt, y = exp(eta.obs), color = "Observed"))


