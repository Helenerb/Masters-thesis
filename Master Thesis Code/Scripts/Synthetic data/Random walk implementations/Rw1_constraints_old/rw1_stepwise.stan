//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

functions {
  real log_gamma_lpdf(real theta, real a, real b){
    return a*log(b) - lgamma(a) + theta*a - b*exp(theta);
  }
}

data {
  real y[100];
}

parameters {
  vector[100] eta;
  real theta_eta;
  real theta_y;
  
  real intercept;
}

transformed parameters{
  real tau_eta = exp(theta_eta);
  real tau_y = exp(theta_y);
}

model {
  theta_eta ~ log_gamma(1, 0.00005);
  theta_y ~ log_gamma(1, 0.00005);
  
  intercept ~ normal(0, 1/sqrt(0.001));
  
  // random walk impplementation, stepwise 
  eta[1] ~  normal(0, 100);
  eta[2:100] ~ normal(eta[1:99], sqrt(6.0/99.0) * 1/sqrt(tau_eta));
  
  sum(eta) ~ normal(0, 0.001*100);
  
  y ~ normal(eta, 1/sqrt(tau_y));
}

