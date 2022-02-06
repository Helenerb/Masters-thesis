// This STAN model implements the traditional lee-carter model, with the following
// model choices:

// gaus: Gaussian Lee-Carter model
// lin: Linear model (no beta term))
// gp: gamma priors
// rw1: alpha and kappa as rw1
// sc: soft constraints

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
  
  real log_mr[X*T];
  
  real nx;  //  number of age steps as a float
  real nt;  // number of period steps as a float
}

parameters {
  real intercept;  // included intercept
  
  vector[X] alpha_us;  //  before sum-to-zero is implemented
  vector[T] kappa_us;  //  before sum-to-zero is implemented
  
  // log-precisions: 
  real theta_alpha;
  real theta_kappa;
  real theta_epsilon;
  
}

transformed parameters {
  real tau_alpha = exp(theta_alpha);
  real tau_kappa = exp(theta_kappa);
  real tau_epsilon = exp(theta_epsilon);
  
  vector[X] alpha = (1/sqrt(tau_alpha))*alpha_us;
  vector[X] kappa = (1/sqrt(tau_kappa))*kappa_us;
  
  vector[X*T] eta = rep_vector(intercept, X*T) + alpha[x] + kappa[t];
}

model {
  // hyperparameters:
  theta_alpha ~ log_gamma(1, 0.00005);
  theta_kappa ~ log_gamma(1, 0.005);  // prior for higher variance
  theta_epsilon ~ log_gamma(1, 0.00005);
  
  intercept ~ normal(0, 1/sqrt(0.001));
  
  //alpha[1] ~ normal(0, 100);  // flat prior
  //alpha[2:X] ~ normal(alpha[1:X-1], 1/sqrt(tau_alpha));
  
  target += -0.5*dot_self(alpha_us[2:X] - alpha_us[1:(X - 1)]);
  sum(alpha_us) ~ normal(0, 0.001*nx);  // soft sum-to-zero
  
  // simpler implementation - kappa as rw1
  // kappa[1] ~ normal(0, 100);
  // kappa[2:T] ~ normal(kappa[1:T-1], 1/sqrt(tau_kappa));
  
  target += -0.5*dot_self(kappa_us[2:T] - kappa_us[1:(T-1)]);
  sum(kappa_us) ~ normal(0, 0.001*nt);

  // send exp_mr as response in the first place - to get identical model to inlabru
  
  log_mr ~ normal(eta, 1/sqrt(tau_epsilon));
}

