// poiss: poisson lee-carter
// lin: linear (no beta term)
// gp: gamma priors for hyperparameters
// rw1: alpha and kappa as random walks

// The input data is a vector 'y' of length 'N'.
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
  
  real nx;  //  number of age steps as a float
  real nt;  // number of period steps as a float
  
  int Y[X*T];   // observed deaths -> response variable
  vector[X*T] E;   // observed population at risk
}

parameters {
  // log-precisions:
  real theta_alpha;
  // real theta_beta;
  real theta_kappa;
  real theta_epsilon;
  
  //real intercept;  // included intercept
  
  vector[X] alpha;
  //vector[X] beta;  
  vector[T] kappa;
  vector[X*T] epsilon;  // unconstrained
  
}

transformed parameters {
  real tau_alpha = exp(theta_alpha);
  //real tau_beta = exp(theta_beta);
  real tau_kappa = exp(theta_kappa);
  real tau_epsilon = exp(theta_epsilon);  // precision of normal distribution of epsilon - error term
  
  vector[X*T] eta = alpha[x] + kappa[t] + epsilon;
}

model {
  theta_alpha ~ log_gamma(1, 0.00005);
  // theta_beta ~ log_gamma(1, 0.00005);
  theta_kappa ~ log_gamma(1, 0.005);
  theta_epsilon ~ log_gamma(1, 0.00005);
  
  //intercept ~ normal(0, 1/sqrt(0.001));
  
  alpha[1] ~ normal(0, 100);  // flat prior
  alpha[2:X] ~ normal(alpha[1:X-1], 1/sqrt(tau_alpha));
  
  //sum(alpha) ~ normal(0, 0.001*nx);  // soft sum-to-zero
  
  // simpler implementation - kappa as rw1
  kappa[1] ~ normal(0, 100);
  kappa[2:T] ~ normal(kappa[1:T-1], 1/sqrt(tau_kappa));
  
  sum(kappa) ~ normal(0, 0.001*nt);
  
  epsilon ~ normal(0, 1/sqrt(tau_epsilon));
  
  Y ~ poisson_log(log(E) + eta);
}

