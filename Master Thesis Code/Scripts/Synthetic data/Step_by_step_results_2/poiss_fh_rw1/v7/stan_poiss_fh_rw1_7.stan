// poiss: poisson lee-carter
// fh: fixed hyperparameters
// rw1: alpha and kappa as random walks

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
  //real intercept;  // included intercept
  
  vector[X] alpha;
  vector[X] beta;  
  vector[T] kappa;
  vector[X*T] epsilon;  // unconstrained
  
}

transformed parameters {
  real tau_alpha = 1.96;
  real tau_beta = 202;
  real tau_kappa = 30;
  //real tau_kappa = 336;
  real tau_epsilon = 420; 
  
  vector[X*T] eta = alpha[x] + beta[x].*kappa[t] + epsilon;
}

model {
  //intercept ~ normal(0, 1/sqrt(0.001));
  
  alpha[1] ~ normal(0, 100);  // flat prior
  alpha[2:X] ~ normal(alpha[1:X-1], 1/sqrt(tau_alpha));
  //sum(alpha) ~ normal(0, 0.001*nx);  // soft sum-to-zero
  
  beta ~ normal(0, 1/sqrt(tau_beta));
  sum(beta) ~ normal(1, 0.001*nx);
  
  // simpler implementation - kappa as rw1
  kappa[1] ~ normal(0, 100);
  kappa[2:T] ~ normal(kappa[1:T-1], 1/sqrt(tau_kappa));
  
  sum(kappa) ~ normal(0, 0.001*nt);
  
  epsilon ~ normal(0, 1/sqrt(tau_epsilon));
  
  Y ~ poisson_log(log(E) + eta);
}

