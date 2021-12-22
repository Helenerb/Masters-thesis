// Implements rw1 with target+= with "constant" term

// data {
//   real y[100];
// }
// 
// parameters {
//   vector[100] eta;
// }
// 
// transformed parameters {
//   real tau_eta = 75;
//   real tau_y = 4.5;
// }
// 
// model {
//   target += 98/2*log(tau_eta)-0.5*tau_eta*dot_self(eta[2:100] - eta[1:99]);
//   
//   // sum-to-zero constraint
//   sum(eta) ~ normal(0, 0.001*100);
//   
//   y ~ normal(eta, 1/sqrt(tau_y));
// }

data {
  int y[100];
  vector[100] E;
}

parameters {
  vector[100] eta;
}

transformed parameters {
  real tau_eta = 75;
  //real tau_y = 4.5;
}

model {
  target += 98/2*log(tau_eta)-0.5*tau_eta*dot_self(eta[2:100] - eta[1:99]);
  
  // sum-to-zero constraint
  sum(eta) ~ normal(0, 0.001*100);
  
  y ~ poisson_log(log(E) + eta);
}

