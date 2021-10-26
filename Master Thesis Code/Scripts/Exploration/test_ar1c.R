# try to fit example script for ar1c process with inlabru

library(inlabru)

n = 500
N = n+1
phi = 0.98
sd.h = 0.4
prec = (1/sd.h)^2
prec.prime = prec / (1-phi^2)
beta = 0.05
h = numeric(n)
y = numeric(n)
z = numeric(n)
s = 0.01
z[n] = NA ## not used
h[1] = rnorm(1, sd = sqrt(1/prec))
y[1] = rnorm(1, sd = sd.h) + rnorm(1, sd = s)
for(i in 2:n) {
  z[i-1] = rnorm(1)
  h[i] = phi * h[i-1] + beta * z[i-1] + rnorm(1, sd = sqrt(1/prec.prime))
  y[i] = h[i] + rnorm(1, sd = s)
}
idx = 1:n

components = ~ -1 +
  kappa(idx, model="ar1c", values = 1:N, args.ar1c = list(Z = cbind(z), Q.beta = matrix(1, 1, 1)))

formula = y ~ kappa

likelihood = like(formula = formula, family = "gaussian", data = data.frame(y, idx))

control.family = list(hyper = list(prec = list(initial = log(1/s^2),fixed=TRUE)))

res.inlabru = bru(components = components,
                  likelihood, 
                  options = list(verbose = T,
                                 bru_verbose = 1, 
                                 num.threads = "1:1",
                                 bru_max_iter=30,
                                 control.family = control.family
                  ))

plot(y, res.inlabru$summary.random$kappa$mean[1:500])
plot(y)

par(mfrow=c(2, 2))
plot(idx, y, type="l", main = "data")
plot(inla.tmarginal(function(x) sqrt(1/exp(x)),
                    r$internal.marginals.hyperpar[[1]]),
     type = "l", lwd=2, main = "sd(h)")
abline(v = sd.h)
plot(inla.smarginal(r$marginals.hyperpar[[2]]),
     type = "l", lwd=2, main = "phi")
abline(v = phi)
## the N’th element is ’beta’
plot(inla.smarginal(r$marginals.random$idx[[N]]), type="l", lwd=2,
     main = "beta")
abline(v = beta)
