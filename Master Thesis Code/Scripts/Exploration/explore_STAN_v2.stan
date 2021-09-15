// This is the second attempt at successful model implementatioj in STAN. 
// This file implements a Poisson-version of the Lee-Carter model, such as 
// the one described in the project thesis as model 2.1. 
// This STAN program is run in the R script explore_STAN_v2.R. 

// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

data {
  int<lower=0> T;  // number of periods
  int<lower=0> X;  // number of age groups
  int x[X*T];   // age group indices in the data
  int t[X*T];   // period indices in the data
  int Y[X*T];   // observed deaths -> response variable
  vector[X*T] E;   // observed population at risk
}

// The parameters accepted by the model. Our model
// accepts the parameters tau_alpha, tau_beta, tau_kappa, tau_epsilon, phi, 
// alpha, beta_raw, kappa_raw, epsilon. 
parameters {
  //real tau_alpha;  // precision of random walk, alpha
  real<lower=0> sigma_alpha;  // sd of random walk, alpha
  //real tau_beta;   // precision of normal distribution of beta
  real<lower=0> sigma_beta;  // sd of normal distribution of beta
  //real tau_kappa;  // precision of random walk, kappa
  real<lower=0> sigma_kappa;  // sd of random walk, kappa
  // real tau_gamma;
  //real tau_epsilon;
  real<lower=0> sigma_epsilon;  // sd of normal distribution of epsilon
  real phi;        // drift of random walk, kappa
  
  real intercept;
  
  vector[X-1] alpha_raw;
  vector[X-1] beta_raw;  // one of the methods for sum-to-zero/unit
  vector[T] kappa;  // unconstrained - random walk with drift
  // vector[T*X -1 ] gamma_raw;  // will sum-to-zero
  vector[X*T] epsilon;
}

transformed parameters {
  vector[X] alpha = append_row(alpha_raw, 0 - sum(alpha_raw));
  vector[X] beta = append_row(beta_raw, 1 - sum(beta_raw));
  //vector[T] kappa = append_row(kappa_raw, 0 - sum(kappa_raw));
  // vector[X*T] gamma = append_row(gamma_raw, -sum(gamma_raw));
  //vector[T*X] eta = alpha + beta*kappa + epsilon; // unsure if this is correct implementation
  
  vector[X*T] eta = rep_vector(intercept, X*T) + alpha[x] + beta[x].*kappa[t] + epsilon;
}

model {
  // hyperpriors: attempt with exponential priors on the standard deviations
  target += exponential_lpdf(sigma_alpha | 10);
  target += exponential_lpdf(sigma_beta | 20);
  target += exponential_lpdf(sigma_kappa | 10);
  target += exponential_lpdf(sigma_epsilon | 100);
  
  // prior distributions
  //alpha[2:X] ~ normal(alpha[1:X-1], 1/sqrt(tau_alpha));
  target += normal_lpdf(alpha[2:X]| alpha[1:X-1], sigma_alpha);
  //beta ~ normal(0, 1/sqrt(tau_beta));
  target += normal_lpdf(beta | 0, sigma_beta);
  //kappa[2:T] ~ normal(phi + kappa[1:T-1], 1/sqrt(tau_kappa)); 
  target += normal_lpdf(kappa[2:T] | phi + kappa[1:T-1], sigma_beta);
  
  //epsilon ~ normal(0, 1/sqrt(tau_epsilon));
  target += normal_lpdf(epsilon | 0, sigma_epsilon);
  
  // We might have to explicitly handle the matrix format
  //Y ~ poisson_log(E .* eta);  // elementwise multiplication
  target += poisson_log_lpmf(Y | log(E) + eta);
}
