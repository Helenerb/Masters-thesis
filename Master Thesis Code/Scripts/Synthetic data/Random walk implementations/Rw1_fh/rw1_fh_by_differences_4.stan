//

data {
  real y[100];
}

parameters {
  vector[100] eta_scaled;
}

transformed parameters {
  real theta_eta = log(7500);
  real theta_y = log(4.5);
  
  real tau_eta = exp(theta_eta);
  real tau_y = exp(theta_y);
}

model {
  //target += -0.5*dot_self(eta_scaled[2:100] - eta_scaled[1:99]);
  target += -0.5*dot_self(eta_scaled[2:100] - eta_scaled[1:99]);
  
  // sum-to-zero constraint
  sum(eta_scaled) ~ normal(0, 0.001*100);
  
  y ~ normal(eta_scaled*1/sqrt(tau_eta), 1/sqrt(tau_y));
}

generated quantities {
  vector[100] eta = eta_scaled*1/sqrt(tau_eta);
}
