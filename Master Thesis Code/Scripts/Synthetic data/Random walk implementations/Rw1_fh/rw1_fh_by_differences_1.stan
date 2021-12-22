
data {
  real y[100];
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[100] eta;
}

transformed parameters{
  real tau_eta = 7500;
  real tau_y = 4.5;
  
  vector[100 - 1] eta_diff = eta[2:100] - eta[1:99];
}

model {
  
  eta_diff ~ normal(0, 1/sqrt(tau_eta));
  
  sum(eta) ~ normal(0, 0.001*100);
  
  y ~ normal(eta, 1/sqrt(tau_y));
}

