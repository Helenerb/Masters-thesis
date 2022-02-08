# Masters-thesis
This repository contains the code used to produce the results of my master's thesis, "Bayesian Mortality Modeling with Linearized Integrated Nested Laplace Approximations"

## Programming Language 
All scripts are written in the R programming language

## Working Directory
Unless otherwise specified, the scripts assume that the working directory is set to "~/Master\ Thesis\ Code"

## Installation
The script "Scripts/installation.R"

## Data
The data is stored in the folder "Data".
The script "Functions/formatters.R" contains the code for formatting and storing the data in the correct format

## Comparison of Random Walk models
The folder "Synthetic\ data/Random\ walk\ implementations/Rw1" contains the code for comparison of random walk models in Stan and inlabru

## Comparison of Constraints
The folder "Synthetic\ data/Random\ walk\ implementations/Rw1_constraints" contains the code for comparison of constraints in Stan and inlabru

## Simulation Study
The folder "Scripts/Synthetic\ data/Step_by_step_resutls_2" contains the code used in the Simulation study. The code for producing the results for each of the models are ordered in subfolders as:

* G1lin: gauss_lin_fh_rw1
* G2lin: gauss_lin_gp_rw1
* G1: gauss_fh_rw1
* G2: gauss_gp_rw1
* P1lin: poiss_lin_fh_rw1
* P2lin poiss_lin_gp_rw1
* P1: poiss_fh_rw1
* P2: poiss_gp_rw1

## Application to German Cancer Data
The folder "Scripts/Real\ data/Analyses" contains the code used in the chapter "Application to German Cancer Data". The code for producing the results for each of the sections in the chapter are ordered in subfolders as:

### Application to the full data
* Female lung cancer: lung_female_full
* Male lung cancer: lung_male_full
* Female stomach cancer: stomach_female_full
* Male stomach cancer: stomach_male_full

### PLCC model without an error term
* Male lung cancer: lung_male_full_no_error

### PLCC model for forecasting
* Female lung cancer: lung_female_predict
* Male lung cancer: lung_male_predict
* Female stomach cancer: stomach_female_predict
* Male stomach cancer: stomach_male_predict

