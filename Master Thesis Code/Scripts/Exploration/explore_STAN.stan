//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> T;
  int<lower=0> X;
  real Y[(X+1)*(T+1)];
  real E[(X+1)*(T+1)];
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real tau_alpha;
  real tau_beta;
  real tau_kappa;
  real tau_gamma;
  // real tau_epsilon;
  
  vector[X] alpha;
  vector[X-1] beta_raw;  // one of the methods for sum-to-zero/unit
  vector[T-1] kappa_raw;  // will sum-to-zero
  vector[(T+1)*(X+1) -1 ] gamma_raw;  // will sum-to-zero
  vector[(T+1)*(X+1)] epsilon;
}

transformed parameters {
  vector[X] beta = append_row(beta_raw, 1 - sum(beta_raw));
  vector[T] kappa = append_row(kappa_raw, -sum(kappa_raw));
  vector[(X+1)*(T+1)] gamma = append_row(gamma_raw, -sum(gamma_raw));
  vector[(T+1)*(X+1)] eta;
  // TODO: give eta the correct value --> on vector form. 
  // eta = alpha + beta*kappa + gamma + epsilon
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  //TODO: Y ~ poisson(eta);
}

