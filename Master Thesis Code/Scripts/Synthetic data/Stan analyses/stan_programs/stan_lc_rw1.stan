 // Stan implementation of lee-carter model (without cohort), where the period
// effect is modelled as a first order random walk. sum-to-zero constraints
// are used for the period effect. 

// define our own log-gamma 
functions {
  real log_gamma_lpdf(real theta, real a, real b){
    //return a*log(b) - lgamma(a) + theta*(a - 1) - b*exp(theta);
    return a*log(b) - lgamma(a) + theta*a - b*exp(theta);
  }
}


data {
  int<lower=0> T;  // number of periods
  int<lower=0> X;  // number of age groups
  int x[X*T];   // age group indices in the data
  int t[X*T];   // period indices in the data
  vector[X*T] ts; // vector containing period indices
  int Y[X*T];   // observed deaths -> response variable
  vector[X*T] E;   // observed population at risk
  real nx;  //  number of age steps as a float
  real nt;  // number of period steps as a float
}

// The parameters accepted by the model. Our model
// accepts the parameters theta_alpha, theta_beta, theta_kappa, theta_epsilon,
// intercept, 
// alpha_raw, beta_raw, kappa_raw, epsilon. 

parameters {
  // log-precisions for hyperparameters
  real theta_alpha;  //  log-precision of alpha
  real theta_beta;  // log-precision of beta
  real theta_kappa;  // log-precision of kappa
  real theta_epsilon;  // log-precision of epsilon
  
  real intercept;  // included intercept
  
  // soft sum-to-zero
  vector[X] alpha;  //  before sum-to-zero is implemented
  vector[X] beta;  //  before sum-to-unit is implemented
  vector[T] kappa;  //  before sum-to-zero is implemented
  
  vector[X*T] epsilon;  // unconstrained
}


transformed parameters {
  // transform log-precisions to precisions
  real tau_alpha = exp(theta_alpha);
  real tau_beta = exp(theta_beta);
  real tau_kappa = exp(theta_kappa);
  real tau_epsilon = exp(theta_epsilon);
  
  // impose constraints
  //vector[X] alpha = alpha_raw - mean(alpha_raw);
  //vector[X] beta = beta_raw - mean(beta_raw) + 1/nx;
  //vector[T] kappa = kappa_raw - mean(kappa_raw);
  
  // construct linear predictor
  vector[X*T] eta = rep_vector(intercept, X*T) + alpha[x] + beta[x].*kappa[t] + epsilon;
}

model {
  // log-gamma priors on log-precisions
  theta_alpha ~ log_gamma(1, 0.00005);
  theta_beta ~ log_gamma(1, 0.00005);
  theta_kappa ~ log_gamma(1, 0.005);
  theta_epsilon ~ log_gamma(1, 0.00005);
  
  intercept ~ normal(0, 1/sqrt(0.001));
  
  // TODO: add prior on first element of random walk
  alpha[1] ~ normal(0, 100);
  alpha[2:X] ~ normal(alpha[1:X-1], 1/sqrt(tau_alpha));
  
  sum(alpha) ~ normal(0, 0.001*nx);
  
  beta ~ normal(0, 1/sqrt(tau_beta));
  sum(beta) ~ normal(1, 0.001*nx);
  
  // random walk of order two
  kappa[1] ~ normal(0, 100);
  kappa[2:T] ~ normal(kappa[1:T-1], 1/sqrt(tau_kappa));
  
  sum(kappa) ~ normal(0, 0.001*nt);
  
  epsilon ~ normal(0, 1/sqrt(tau_epsilon));
  
  Y ~ poisson_log(log(E) + eta);
}
