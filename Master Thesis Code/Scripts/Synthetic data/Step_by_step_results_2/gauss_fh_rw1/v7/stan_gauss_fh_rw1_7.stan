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
  //real intercept;  // included intercept
  
  vector[X] alpha;
  vector[X] beta;
  vector[T] kappa;
}

transformed parameters {
  real tau_alpha = 1.96;
  real tau_beta = 202;
  //real tau_kappa = 336;
  real tau_kappa = 30;
  real tau_epsilon = 420; 
  
  vector[X*T] eta = alpha[x] + beta[x].*kappa[t];
}

model {
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

