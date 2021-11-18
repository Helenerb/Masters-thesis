// This STAN model implements the traditional lee-carter model

// The input data is a vector 'y' of length 'N'.
functions {
  real log_gamma_lpdf(real theta, real a, real b){
    return a*log(b) - lgamma(a) + theta*(a - 1) - b*exp(theta);
  }
}


data {
  int<lower=0> T;  // number of periods
  int<lower=0> X;  // number of age groups
  
  int x[X*T];   // age group indices in the data
  int t[X*T];   // period indices in the data
  
  real exp_mortality_rate[X*T];   // observed mortality - response variable
  
  real nx;  //  number of age steps as a float
  real nt;  // number of period steps as a float
}

parameters {
  real theta_alpha;  // log precision of alpha
  real theta_beta;  // log-precision of beta
  real theta_kappa;  // log-precision of kappa
  real theta_epsilon;  // log-precision of epsilon
  
  real intercept;  // included intercept
  
  vector[X] alpha_raw;  //  before sum-to-zero is implemented
  vector[X] beta_raw;  //  before sum-to-unit is implemented
  
  vector[T] kappa_raw;  //  before sum-to-zero is implemented
  
  //vector[X*T] epsilon;  // unconstrained
}

transformed parameters {
  //vector[X] alpha = append_row(alpha_raw, 0 - sum(alpha_raw));
  //vector[X] beta = append_row(beta_raw, 1 - sum(beta_raw));
  real tau_alpha = exp(theta_alpha);  // precision of random walk, alpha
  real tau_beta = exp(theta_beta);  // precision of normal distribution of beta
  real tau_kappa  = exp(theta_kappa);  // precision of random walk, kappa
  real tau_epsilon = exp(theta_epsilon);  // precision of normal distribution of epsilon - error term
  
  vector[X] alpha = alpha_raw - mean(alpha_raw);
  vector[X] beta = beta_raw - mean(beta_raw) + 1/nx;
  vector[T] kappa = kappa_raw - mean(kappa_raw);
  
  // split kappa into driftless rw and linear term
  //vector[X*T] eta = rep_vector(intercept, X*T) + alpha[x] + beta[x].*kappa[t] + epsilon;
  vector[X*T] eta = rep_vector(intercept, X*T) + alpha[x] + beta[x].*kappa[t];
}

model {
  // hyperpriors: attempt with exponential priors on the standard deviations
  
  //target += gamma_lpdf(tau_alpha | 1, 0.00005);
  //tau_alpha ~ gamma(1, 0.00005);
  theta_alpha ~ log_gamma(1, 0.00005);
  
  //target += gamma_lpdf(tau_beta | 1, 0.00005);
  //tau_beta ~ gamma(1, 0.00005);
  theta_beta ~ log_gamma(1, 0.00005);
  
  //target += gamma_lpdf(tau_kappa | 1, 0.005);  // more informative prior on tau_kappa
  //tau_kappa ~ gamma(1, 0.005);
  theta_kappa ~ log_gamma(1, 0.005);
  
  //target += gamma_lpdf(tau_epsilon | 1, 0.00005);
  //tau_epsilon ~ gamma(1, 0.00005);
  theta_epsilon ~ log_gamma(1, 0.00005);
  
  //target += normal_lpdf(intercept | 0, 1/sqrt(0.001));
  intercept ~ normal(0, 1/sqrt(0.001));
  
  // prior distributions
  //target += normal_lpdf(alpha_raw[2:X-1]| alpha_raw[1:X-2], 1/sqrt(tau_alpha));
  //target += normal_lpdf(beta_raw | 0, 1/sqrt(tau_beta));
  
  //target += normal_lpdf(alpha_raw[2:X] | alpha_raw[1:X-1], 1/sqrt(tau_alpha));
  
  alpha_raw[1] ~ normal(0, 1/sqrt(tau_alpha));
  alpha_raw[2:X] ~ normal(alpha_raw[1:X-1], 1/sqrt(tau_alpha));
  
  //target += normal_lpdf(beta_raw | 0, 1/sqrt(tau_beta));
  beta_raw ~ normal(0, 1/sqrt(tau_beta));  // Should I explicitly say that this is a vector?
  
  // random walk of order two
  //target += normal_lpdf(kappa_raw[1] | 0, 1/sqrt(tau_kappa));
  kappa_raw[1] ~ normal(0, 1/sqrt(tau_kappa));
  
  //target += normal_lpdf(kappa_raw[2] | kappa_raw[1], 1/sqrt(tau_kappa));
  kappa_raw[2] ~ normal(kappa_raw[1], 1/sqrt(tau_kappa));
  
  //target += normal_lpdf(kappa_raw[3:T] | 2*kappa_raw[2:T-1] - kappa_raw[1:T-2], 1/sqrt(tau_kappa));
  kappa_raw[3:T] ~ normal(2*kappa_raw[2:T-1] - kappa_raw[1:T-2], 1/sqrt(tau_kappa));
  
  //epsilon ~ normal(0, 1/sqrt(tau_epsilon));
  //target += normal_lpdf(epsilon | 0, 1/sqrt(tau_epsilon));
  
  // We might have to explicitly handle the matrix format
  //Y ~ poisson_log(E .* eta);  // elementwise multiplication
  //target += normal_lupdf(exp_mortality_rate | eta, rep_vector(1/sqrt(tau_epsilon), X*T);
  exp_mortality_rate ~ normal(eta, 1/sqrt(tau_epsilon));
}

