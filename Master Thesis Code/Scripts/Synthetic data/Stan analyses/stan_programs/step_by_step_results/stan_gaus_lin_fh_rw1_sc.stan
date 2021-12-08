// This STAN model implements the traditional lee-carter model

// The input data is a vector 'y' of length 'N'.
// functions {
//   real log_gamma_lpdf(real theta, real a, real b){
//     //return a*log(b) - lgamma(a) + theta*(a - 1) - b*exp(theta);
//     return a*log(b) - lgamma(a) + theta*a - b*exp(theta);
//   }
// }


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
  real intercept;  // included intercept
  
  vector[X] alpha;  //  before sum-to-zero is implemented
  //vector[X] beta;  //  before sum-to-unit is implemented
  vector[T] kappa;  //  before sum-to-zero is implemented
  
  //vector[X*T] epsilon;  // unconstrained
}

transformed parameters {
  //vector[X] alpha = append_row(alpha_raw, 0 - sum(alpha_raw));
  //vector[X] beta = append_row(beta_raw, 1 - sum(beta_raw));
  real tau_alpha = 1.96;
  //real tau_beta = 202;
  real tau_kappa = 336;
  real tau_epsilon = 420;  // precision of normal distribution of epsilon - error term
  
  // split kappa into driftless rw and linear term
  //vector[X*T] eta = rep_vector(intercept, X*T) + alpha[x] + beta[x].*kappa[t] + epsilon;
  vector[X*T] eta = rep_vector(intercept, X*T) + alpha[x] + kappa[t];
}

model {
  intercept ~ normal(0, 1/sqrt(0.001));
  
  
  alpha[1] ~ normal(0, 100);  // flat prior
  alpha[2:X] ~ normal(alpha[1:X-1], 1/sqrt(tau_alpha));
  
  sum(alpha) ~ normal(0, 0.001*nx);  // soft sum-to-zero
  
  //beta ~ normal(0, 1/sqrt(tau_beta));
  
  //sum(beta) ~ normal(1, 0.001*nx);
  
  // random walk of order two
  //kappa_raw[1] ~ normal(0, 1/sqrt(tau_kappa));
  //kappa_raw[2] ~ normal(kappa_raw[1], 1/sqrt(tau_kappa));
  //kappa_raw[3:T] ~ normal(2*kappa_raw[2:T-1] - kappa_raw[1:T-2], 1/sqrt(tau_kappa));
  
  // simpler implementation - kappa as rw1
  kappa[1] ~ normal(0, 100);
  kappa[2:T] ~ normal(kappa[1:T-1], 1/sqrt(tau_kappa));
  
  sum(kappa) ~ normal(0, 0.001*nt);
  
  //epsilon ~ normal(0, 1/sqrt(tau_epsilon));
  
  // send exp_mr as response in the first place - to get identical model to inlabru
  
  log_mr ~ normal(eta, 1/sqrt(tau_epsilon));
}

