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
}

transformed parameters {
  real tau_eta = 7500;
  real tau_y = 4.5;
}

model {
  
  // random walk of order two implementation, stepwise
  eta[1] ~  normal(0, 100);
  eta[2] ~ normal(eta[1], 100);
  eta[3:100] ~ normal(2*eta[2:99] - eta[1:98], 1/sqrt(tau_eta));
  
  sum(eta) ~ normal(0, 0.001*100);
  
  y ~ normal(eta, 1/sqrt(tau_y));
}

