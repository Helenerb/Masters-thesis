# script for readying synthetic lung cancer above a45  - v7

# note: assumes working directory at .../Master Thesis Code

synthetic.male.lung.a45.v7 <- function(){
  obs <- read.csv("Data/synthetic_male_lung_7.csv")
  obs <- obs %>% mutate(x.old = x, x = x - 9)
  
  underlying.effects <- list(obs = obs, nx = 9, nt = 18,
                               alpha.true = {obs %>% filter(t == 0)}$alpha,
                               beta.true = {obs %>% filter(t == 0)}$beta,
                               kappa.true = {obs %>% filter(x == 9)}$kappa,
                               intercept = unique(obs$intercept),
                               tau.alpha.true = unique(obs$tau.alpha),
                               tau.beta.true = unique(obs$tau.beta),
                               tau.kappa.true = unique(obs$tau.kappa),
                               tau.epsilon.true = unique(obs$tau.epsilon))
  return(list(obs = obs, underlying.effects = underlying.effects))
}