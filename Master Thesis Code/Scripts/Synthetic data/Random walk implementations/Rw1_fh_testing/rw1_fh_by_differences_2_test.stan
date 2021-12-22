//

data {
  real y[100];
}

parameters {
  vector[100] eta_scaled;
}

transformed parameters {
  real tau_eta = 75;
  real tau_y = 4.5;
}

model {
  target += -0.5*dot_self(eta_scaled[2:100] - eta_scaled[1:99]);
  
  // sum-to-zero constraint
  sum(eta_scaled) ~ normal(0, 0.001*100);
  
  y ~ normal(eta_scaled*1/sqrt(tau_eta), 1/sqrt(tau_y));
}

generated quantities {
  vector[100] eta = eta_scaled*1/sqrt(tau_eta);
}

