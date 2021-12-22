//

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
  //target += -0.5*dot_self(eta_scaled[2:100] - eta_scaled[1:99]);
  target += 98/2*log(tau_eta)-0.5*tau_eta*dot_self(eta[1:98] - 2*eta[2:99] + eta[3:100]);
  
  // sum-to-zero constraint
  sum(eta) ~ normal(0, 0.001*100);
  
  y ~ normal(eta, 1/sqrt(tau_y));
}

