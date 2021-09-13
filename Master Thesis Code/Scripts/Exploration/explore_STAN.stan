// Information about the program
data {
  int<lower=0> T;
  int<lower=0> X;
  real Y[(X+1)*(T+1)];
  real E[(X+1)*(T+1)];
}

// The parameters accepted by the model. Our model
// accepts the parameters tau_alpha, tau_beta, tau_kappa, tau_epsilon, phi, 
// alpha, beta_raw, kappa_raw, epsilon. 
parameters {
  real tau_alpha;  // precision of random walk, alpha
  real tau_beta;   // precision of normal distribution of beta
  real tau_kappa;  // precision of random walk, kappa
  // real tau_gamma;
  real tau_epsilon;
  real phi;        // drift of random walk, kappa
  
  vector[X] alpha;
  vector[X-1] beta_raw;  // one of the methods for sum-to-zero/unit
  vector[T-1] kappa_raw;  // will sum-to-zero
  // vector[T*X -1 ] gamma_raw;  // will sum-to-zero
  // vector[T*X] epsilon;
  matrix[X,T] epsilon;
}

transformed parameters {
  vector[X] beta = append_row(beta_raw, 1 - sum(beta_raw));
  row_vector[T] kappa = append_row(kappa_raw, -sum(kappa_raw));
  // vector[X*T] gamma = append_row(gamma_raw, -sum(gamma_raw));
  //vector[T*X] eta = alpha + beta*kappa + epsilon; // unsure if this is correct implementation
  matrix[X,T] eta = rep_array(alpha, T) + beta*kappa + epsilon;
}

model {
  // hyperpriors: use default distributions
  
  // prior distributions
  alpha[2:X] ~ normal(alpha[1:X-1], 1/sqrt(tau_alpha));
  beta ~ normal(0, 1/sqrt(tau_beta));
  kappa[2:T] ~ normal(phi + kappa[1:T-1], 1/sqrt(tau_kappa)); 
  epsilon ~ normal(0, 1/sqrt(tau_epsilon))
  
  Y ~ poisson_log(E .* eta);  // elementwise multiplocation
}

