// For testing purposes!!!
functions {
  real log_gamma_lpdf(real theta, real a, real b){
    //return a*log(b) - lgamma(a) + theta*(a - 1) - b*exp(theta);
    return a*log(b) - lgamma(a) + theta*a - b*exp(theta);
  }
}

data {
  real y[100];
}

parameters {
  //real<lower=0> tau_eta;
  real theta_eta;
  real<lower=0> tau_y;
  
  vector[100] eta;
}

transformed parameters{
  tau_eta = exp(theta_eta);
}

model {
  //tau_eta ~ gamma(1, 0.00005);
  theta_eta ~ log_gamma(1, 0.00005);
  tau_y ~ gamma(1, 0.00005);
  
  //target += -0.5*dot_self(eta_scaled[2:100] - eta_scaled[1:99]);
  target += 99.0*0.5*log(tau_eta)-0.5*tau_eta*dot_self(eta[2:100]  - eta[1:99]);
  
  // sum-to-zero constraint
  sum(eta) ~ normal(0, 0.001*100);
  
  y ~ normal(eta, 1/sqrt(tau_y));
}

