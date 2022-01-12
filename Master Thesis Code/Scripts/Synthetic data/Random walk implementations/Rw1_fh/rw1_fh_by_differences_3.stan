//

data {
  real y[100];
}

parameters {
  vector[100] eta;
}

transformed parameters {
  real tau_eta = exp(2);
  real tau_y = 4;
}

model {
  // Note: There was an error in this implementation earlier! (98/2) (which is the case for rw2 not rw1)
  //target += 98/2*log(tau_eta)-0.5*tau_eta*dot_self(eta[2:100] - eta[1:99]);
  target += -99.0*0.5*log(6.0/99.0) + 99.0*0.5*log(tau_eta) - 0.5*(99.0/6.0)*tau_eta*dot_self(eta[2:100]  - eta[1:99]);
  
  // sum-to-zero constraint
  //sum(eta) ~ normal(0, 0.001*100);
  
  y ~ normal(eta, 1/sqrt(tau_y));
}

