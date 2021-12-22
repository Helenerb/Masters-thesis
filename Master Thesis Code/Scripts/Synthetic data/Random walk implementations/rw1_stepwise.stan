//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

data {
  real y[100];
}

parameters {
  vector[100] eta;
  real<lower=0> tau_eta;
  real<lower=0> tau_y;
}

model {
  tau_eta ~ gamma(1, 0.00005);
  tau_y ~ gamma(1, 0.00005);
  
  // random walk impplementation, stepwise 
  eta[1] ~  normal(0, 100);
  eta[2:100] ~ normal(eta[1:99], 1/sqrt(tau_eta));
  
  sum(eta) ~ normal(0, 0.001*100);
  
  y ~ normal(eta, 1/sqrt(tau_y));
}

