// This is the second attempt at successful model implementatioj in STAN. 
// This file implements a Poisson-version of the Lee-Carter model, such as 
// the one described in the project thesis as model 2.1. 
// This STAN program is run in the R script explore_STAN_v2.R. 

// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// OK
data {
  int<lower=0> T;  // number of periods
  int<lower=0> X;  // number of age groups
  int x[X*T];   // age group indices in the data
  int t[X*T];   // period indices in the data
  vector[X*T] ts; // vector containing period indices
  int Y[X*T];   // observed deaths -> response variable
  vector[X*T] E;   // observed population at risk
}

// The parameters accepted by the model. Our model
// accepts the parameters tau_alpha, tau_beta, tau_kappa, tau_epsilon, phi, 
// alpha, beta_raw, kappa_raw, epsilon. 

// OK
parameters {
  real<lower=0> tau_alpha;  // precision of random walk, alpha
  real<lower=0> tau_beta;  // precision of normal distribution of beta
  real<lower=0> tau_kappa;  // precision of random walk, kappa
  real<lower=0> tau_epsilon;  // precision of normal distribution of epsilon
  
  real phi;        // drift of random walk, kappa
  
  real intercept;  // included intercept
  
  vector[X] alpha_raw;  //  before sum-to-zero is implemented
  vector[X] beta_raw;  //  before sum-to-unit is implemented
  
  vector[T] kappa_raw;  //  before sum-to-zero is implemented
  
  vector[X*T] epsilon;  // unconstrained
}


// OK
transformed parameters {
  //vector[X] alpha = append_row(alpha_raw, 0 - sum(alpha_raw));
  //vector[X] beta = append_row(beta_raw, 1 - sum(beta_raw));
  
  vector[X] alpha = alpha_raw - mean(alpha_raw);
  vector[X] beta = beta_raw - mean(beta_raw) + 1/80.0;
  vector[T] kappa = kappa_raw - mean(kappa_raw);
  
  // split kappa into driftless rw and linear term
  //vector[X*T] eta = rep_vector(intercept, X*T) + alpha[x] + beta[x].*kappa[t] + epsilon;
  vector[X*T] eta = rep_vector(intercept, X*T) + alpha[x] + beta[x].*kappa[t] + (phi*ts).*beta[x] + epsilon;
}

// OK
model {
  // hyperpriors: attempt with exponential priors on the standard deviations
  
  target += gamma_lpdf(tau_alpha | 1, 0.00005);
  target += gamma_lpdf(tau_beta | 1, 0.00005);
  target += gamma_lpdf(tau_kappa | 1, 0.00005);
  target += gamma_lpdf(tau_epsilon | 1, 0.00005);
  
  // prior distributions
  //target += normal_lpdf(alpha_raw[2:X-1]| alpha_raw[1:X-2], 1/sqrt(tau_alpha));
  //target += normal_lpdf(beta_raw | 0, 1/sqrt(tau_beta));
  
  target += normal_lpdf(alpha_raw[2:X] | alpha_raw[1:X-1], 1/sqrt(tau_alpha));
  target += normal_lpdf(beta_raw | 0, 1/sqrt(tau_beta));
  
  target += normal_lpdf(phi | 0, 1);  //  normal prior for linear term
  
  // add starting point for kappa
  // target += normal_lpdf(kappa[1] | 0, 1/sqrt(tau_kappa));
  
  // change to driftless random walk and linear term
  //target += normal_lpdf(kappa[2:T] | phi + kappa[1:T-1], 1/sqrt(tau_kappa));
  
  // driftless kappa
  target += normal_lpdf(kappa_raw[2:T] | kappa_raw[1:T-1], 1/sqrt(tau_kappa));
  
  //epsilon ~ normal(0, 1/sqrt(tau_epsilon));
  target += normal_lpdf(epsilon | 0, 1/sqrt(tau_epsilon));
  
  // We might have to explicitly handle the matrix format
  //Y ~ poisson_log(E .* eta);  // elementwise multiplication
  target += poisson_log_lpmf(Y | log(E) + eta);
}
