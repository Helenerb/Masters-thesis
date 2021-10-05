source("configuration_v1.R")

underlying.effects <- configuration.v7()

obs <- underlying.effects$obs

nx = length(unique(obs$x))
nt = length(unique(obs$t))
nc = length(unique(obs$c))

#   ----   Inlabru model fitting   ----

A.mat = matrix(1, nrow = 1, ncol = nx)  #  helper values for constraining of beta
e.vec = 1  #  helper values for constraining of beta

pc.prior <- list(prec = list(prior = "pc.prec", param = c(1,0.05)))

comp = ~ -1 +
  Int(1) +
  alpha(x, model = "rw1", values=unique(obs$x), hyper = pc.prior, constr = TRUE) +
  phi(t, model = "linear", prec.linear = 1) +
  beta(x.c, model = "iid", extraconstr = list(A = A.mat, e = e.vec), hyper = pc.prior) +
  kappa(t.c, model = "rw1", values = unique(obs$t), constr = TRUE, hyper = pc.prior) +
  gamma(c, model = "rw1", values = unique(obs$c), constr = TRUE, hyper = pc.prior) +
  epsilon(xt, model = "iid", hyper = pc.prior)

formula = Y ~  Int + alpha + beta*phi + beta*kappa + gamma + epsilon

likelihood = like(formula = formula, family = "poisson", data = obs, E = obs$E)

c.c <- list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE)  # control.compute

res.inlabru = bru(components = comp,
                  likelihood, 
                  options = list(verbose = F,
                                 bru_verbose = 1, 
                                 num.threads = "1:1",
                                 control.compute = c.c,
                                 bru_max_iter=30,
                                 control.predictor = list(link = 1)
                  ))

#   ----   plotting results of inlabru fit   ----
source("plot_inlabru_vs_underlying.R")

plots <- plot.inlabru.vs.underlying.v1(res.inlabru, underlying.effects)
plots$p.alpha
plots$p.beta
plots$p.phi
plots$p.intercept
plots$p.kappa
plots$p.eta
plots$p.eta.2
plots$p.eta.t
plots$p.eta.x
plots$p.gamma

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


