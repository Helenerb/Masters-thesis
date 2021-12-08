# script running only stan - analysis of poiss_gp_rw1_4.R

#   ----   Load libraries and set workspace   ----
set_workspace <- function(markov=TRUE){
  if(markov){
    .libPaths("~/Documents/R_libraries")
    setwd("~/Documents/GitHub/Masteroppgave/Masters-thesis/Master Thesis Code")
  } else {
    setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code")
  }
}

#   ----   TODO: Change the following lines to change to and from Markov  ----
set_workspace(markov=TRUE)
#set_workspace(markov=FALSE)

library("tidyverse")
library("inlabru")
library("ggplot2")
library("INLA")
library("patchwork")
library("rstan")

investigation.name <- "poiss_gp_rw1"

# Path to where files and results are stored
# If you run locally - change this to where you have stored the code
investigation.path <- file.path(investigation.name, "v4")

#   ----    Retrieve the data   ----

synthetic.male.lung.v4 <- function(){
  #TODO: If you run this locally - change to where you have stored the data
  obs <- read.csv("Data/synthetic_male_lung_4.csv")
  
  obs.trad <- obs %>% 
    select(c(x, t, xt, age.int, year, x.c, alpha, beta, kappa, intercept, epsilon,
             eta, tau.alpha, tau.beta, tau.kappa, tau.epsilon, E)) %>%
    mutate(eta.no.error = intercept + alpha + beta*kappa) %>%
    mutate(mr_gaussian = exp(eta)) %>%
    mutate(Y_gaussian  = mr_gaussian * E)
  
  underlying.effects <- list(obs = obs.trad, nx = 18, nt = 18,
                             alpha.true = {obs %>% filter(t == 0)}$alpha,
                             beta.true = {obs %>% filter(t == 0)}$beta,
                             kappa.true = {obs %>% filter(x == 0)}$kappa,
                             intercept = unique(obs$intercept),
                             age.intercept.true = unique(obs$intercept),
                             tau.alpha.true = unique(obs$tau.alpha),
                             tau.beta.true = unique(obs$tau.beta),
                             tau.kappa.true = unique(obs$tau.kappa),
                             tau.epsilon.true = unique(obs$tau.epsilon))
  
  return(list(obs = obs, underlying.effects = underlying.effects))
}

# We use this data for both inlabru and stan
config.data <- synthetic.male.lung.v4()
obs <- config.data$obs
underlying.effects <- config.data$underlying.effects

#   ----   Run STAN analysis   ----
# Running traditional lc version of stan, 
# implemented with log-precisions, random effects as iid, and no constraints

# TODO: If you run this locally - change it to where you have your folder containing these files
stan.output  <- file.path("Scripts/Synthetic data/Step_by_step_results", investigation.path)
source("Scripts/Synthetic\ data/run_stan_functions.R")

run_stan <- function(stan_program, obs, chains, warmup, iter, output.path, config.name, markov=TRUE){
  
  stan_fit <- run_stan_program_lc(
    list(obs = obs), chains=chains,warmup=warmup,
    iter=iter, stan_program=stan_program)
  
  store_stan_results(
    fit=stan_fit, output.path=output.path, config=config.name,
    chains=chains, warmup=warmup, iter=iter, stan_program=stan_program)
}

#TODO: if you run locally, change path to where you have stored stan program
run_stan(
  stan_program="Scripts/Synthetic data/Stan analyses/stan_programs/step_by_step_results/stan_pois_gp_rw1_sc.stan",
  obs = obs, chains=4, warmup = 30, iter = 300, output.path = stan.output,
  config.name = investigation.name, markov=F)