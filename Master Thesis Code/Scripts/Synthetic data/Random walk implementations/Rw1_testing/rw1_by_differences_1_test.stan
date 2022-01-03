// NOTE: This file is for testing purposes!!!

functions {
  real log_gamma_lpdf(real theta, real a, real b){
    //return a*log(b) - lgamma(a) + theta*(a - 1) - b*exp(theta);
    return a*log(b) - lgamma(a) + theta*a - b*exp(theta);
  }
}

data {
  real y[100];
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  //real<lower=0> tau_eta;
  real theta_eta;
  
  real<lower=0> tau_y;
  vector[100] eta;
  
  //real intercept;
}

transformed parameters{
  //real theta_eta = 2;
  real tau_eta = exp(theta_eta);
  vector[100 - 1] eta_diff = eta[2:100] - eta[1:99];
}
// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  theta_eta ~ log_gamma(1, 0.00005);
  
  tau_y ~ gamma(1, 0.00005);
  
  //intercept ~ normal(0, 1/sqrt(0.001));
  
  //eta_diff ~ normal(0, 1/sqrt(tau_eta));
  target += normal_lupdf(eta_diff | 0, 1/sqrt(tau_eta));
  
  //sum(eta) ~ normal(0, 0.001*100);
  
  //y ~ normal(rep_vector(intercept, 100) + eta, 1/sqrt(tau_y));
  y ~ normal(eta, 1/sqrt(tau_y));
}

