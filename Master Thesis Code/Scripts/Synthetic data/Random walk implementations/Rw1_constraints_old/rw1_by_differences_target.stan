//
functions {
  real log_gamma_lpdf(real theta, real a, real b){
    return a*log(b) - lgamma(a) + theta*a - b*exp(theta);
  }
}

data {
  real y[100];
}

parameters {
  real theta_eta;
  real theta_y;
  
  vector[100] eta_scaled;
}

transformed parameters {
  real tau_eta = exp(theta_eta);
  real tau_y = exp(theta_y);
}

model {
  theta_eta ~ log_gamma(1, 0.00005);
  theta_y ~ log_gamma(1, 0.00005);
  
  target += -0.5*dot_self(eta_scaled[2:100] - eta_scaled[1:99]);
  
  // sum-to-zero constraint
  //sum(eta_scaled) ~ normal(0, 0.001*100);
  
  y ~ normal(eta_scaled*(sqrt(6.0/99.0) * 1/sqrt(tau_eta)), 1/sqrt(tau_y));
}

generated quantities {
  vector[100] eta = eta_scaled*(sqrt(6.0/99.0) * 1/sqrt(tau_eta));
}

