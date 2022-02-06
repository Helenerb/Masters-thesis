// This stan program is for testing purposes!

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
  vector[100] eta;
  
  real theta_eta;
  real theta_y;
  real intercept;
}

transformed parameters{
  real tau_eta = exp(theta_eta);
  real tau_y = exp(theta_y);
}

model {
  //tau_eta ~ gamma(1, 0.00005);
  theta_eta ~ log_gamma(1, 0.00005);
  theta_y ~ log_gamma(1, 0.00005);
  
  intercept ~ normal(0, 1/sqrt(0.001));
  
  // random walk impplementation, stepwise 
  eta[1] ~  normal(0, 100);
  eta[2:100] ~ normal(eta[1:99], sqrt(0.06060606) * 1/sqrt(tau_eta));
  
  sum(eta) ~ normal(0, 0.001*100);
  
  y ~ normal(rep_vector(intercept, 100) + eta, 1/sqrt(tau_y));
}
