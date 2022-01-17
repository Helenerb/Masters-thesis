// This STAN model implements the traditional lee-carter model, with the following
// model choices:

// gaus: Gaussian Lee-Carter model
// lin: Linear model (no beta term))
// fh: fixed hypers
// rw1: alpha and kappa as rw1

functions {
  real log_gamma_lpdf(real theta, real a, real b){
    return a*log(b) - lgamma(a) + theta*a - b*exp(theta);
  }
}


data {
  int<lower=0> T;  // number of periods
  int<lower=0> X;  // number of age groups
  
  int x[X*T];   // age group indices in the data
  int t[X*T];   // period indices in the data
  
  real log_mr[X*T];
  
  real nx;  //  number of age steps as a float
  real nt;  // number of period steps as a float
}

parameters {
  // hyperparameters
  real theta_alpha;
  real theta_beta;
  real theta_kappa;
  real theta_epsilon;
  
  // random effects
  vector[X] alpha;
  vector[X] beta;
  vector[T] kappa;
}

transformed parameters {
  real tau_alpha = exp(theta_alpha);
  real tau_beta = exp(theta_beta);
  real tau_kappa = exp(theta_kappa);
  real tau_epsilon = exp(theta_epsilon); 
  
  vector[X*T] eta = alpha[x] + beta[x].*kappa[t];
}

model {
  theta_alpha ~ log_gamma(1, 0.00005);
  theta_beta ~ log_gamma(1, 0.00005);
  theta_kappa ~ log_gamma(1, 0.005);  // prior for higher variance
  theta_epsilon ~ log_gamma(1, 0.00005);
  
  alpha[1] ~ normal(0, 100);  // flat prior
  alpha[2:X] ~ normal(alpha[1:X-1], 1/sqrt(tau_alpha));
  
  beta ~ normal(0, 1/sqrt(tau_beta));
  sum(beta) ~ normal(1, 0.001*nx);
  
  // simpler implementation - kappa as rw1
  kappa[1] ~ normal(0, 100);
  kappa[2:T] ~ normal(kappa[1:T-1], 1/sqrt(tau_kappa));
  
  sum(kappa) ~ normal(0, 0.001*nt);

  // send exp_mr as response in the first place - to get identical model to inlabru
  
  log_mr ~ normal(eta, 1/sqrt(tau_epsilon));
}

