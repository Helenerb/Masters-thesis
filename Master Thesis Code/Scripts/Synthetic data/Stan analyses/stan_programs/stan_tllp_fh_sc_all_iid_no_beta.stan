 // This STAN model implements the traditional lee-carter model, 
// with soft constraints and with fixed precisions. 

// tllp: traditional lc log prec model/ gaussian Lee-Carter model, with hyperpriors defined on the log precisions
// fh: fixed hypers - fixed hyperparameters
// iid: all random effects modelled as iid 
// no_constr: no constraints applied to random effects 

data {
  int<lower=0> T;  // number of periods
  int<lower=0> X;  // number of age groups
  
  int x[X*T];   // age group indices in the data
  int t[X*T];   // period indices in the data
  
  //real mr[X*T];   // observed mortality rate (force of mortality = D/E) - response variable
  real exp_mr[X*T];
  
  real nx;  //  number of age steps as a float
  real nt;  // number of period steps as a float
}

parameters {
  //real theta_alpha;  // log precision of alpha
  //real theta_beta;  // log-precision of beta
  //real theta_kappa;  // log-precision of kappa
  //real theta_epsilon;  // log-precision of epsilon
  
  real intercept;  // included intercept
  
  vector[X] alpha;  //  before sum-to-zero is implemented
  //vector[X] beta;  //  before sum-to-unit is implemented
  vector[T] kappa;  //  before sum-to-zero is implemented
  
  //vector[X*T] epsilon;  // unconstrained
}

transformed parameters {
  //vector[X] alpha = append_row(alpha_raw, 0 - sum(alpha_raw));
  //vector[X] beta = append_row(beta_raw, 1 - sum(beta_raw));
  //real tau_alpha = exp(theta_alpha);  // precision of random walk, alpha
  //real tau_beta = exp(theta_beta);  // precision of normal distribution of beta
  //real tau_kappa  = exp(theta_kappa);  // precision of random walk, kappa
  //real tau_epsilon = exp(theta_epsilon);  // precision of normal distribution of epsilon - error term
  real tau_alpha = 1.96;
  //real tau_beta = 100;
  real tau_kappa = 70;
  real tau_epsilon = 400;
  
  // split kappa into driftless rw and linear term
  //vector[X*T] eta = rep_vector(intercept, X*T) + alpha[x] + beta[x].*kappa[t] + epsilon;
  //vector[X*T] eta = rep_vector(intercept, X*T) + alpha[x] + beta[x].*kappa[t];
  vector[X*T] eta = rep_vector(intercept, X*T) + alpha[x] + kappa[t];
}

model {
  intercept ~ normal(0, 1/sqrt(0.001));
  
  //alpha[1] ~ normal(0, 100);  // flat prior
  //alpha[2:X] ~ normal(alpha[1:X-1], 1/sqrt(tau_alpha));
  
  alpha ~ normal(0, 1/sqrt(tau_alpha));  // alpha as iid
  sum(alpha) ~ normal(0, 0.001*nx);  // soft sum-to-zero
  
  //beta ~ normal(0, 1/sqrt(tau_beta));
  //sum(beta) ~ normal(1, 0.001*nx);  // soft sum-to-one
  
  //kappa[1] ~ normal(0, 100);
  //kappa[2:T] ~ normal(kappa[1:T-1], 1/sqrt(tau_kappa));
  
  kappa ~ normal(0, 1/sqrt(tau_kappa));  // kappa as iid
  sum(kappa) ~ normal(0, 0.001*nt);  // soft sum-to-zero
  
  exp_mr ~ normal(eta, 1/sqrt(tau_epsilon));
}

