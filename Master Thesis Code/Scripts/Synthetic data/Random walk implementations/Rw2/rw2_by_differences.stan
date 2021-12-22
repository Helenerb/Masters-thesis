
data {
  real y[100];
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=0> tau_eta;
  real<lower=0> tau_y;
  vector[100] eta;
}

transformed parameters{
  vector[100 - 2] eta_diff = eta[1:98] - 2*eta[2:99] + eta[3:100];
}

model {
  tau_eta ~ gamma(1, 0.00005);
  tau_y ~ gamma(1, 0.00005);
  
  eta_diff ~ normal(0, 1/sqrt(tau_eta));
  
  sum(eta) ~ normal(0, 0.001*100);
  
  y ~ normal(eta, 1/sqrt(tau_y));
}

