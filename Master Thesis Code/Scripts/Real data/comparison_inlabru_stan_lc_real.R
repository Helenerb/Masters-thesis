# running inlabru and stan with real data! 

setwd("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Synthetic data")

# where should figures produced in this script be stored
path.to.folder = '/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Output/Figures/Stomach cancer'

source("~/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Scripts/Functions/formatters.R")

population <- format_population_data(path_to_data = "/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data/population-germany.xlsx",
                       path_to_output = "/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data/population-germany.Rda")

cancer.data <- format_cancer_data(path_to_data = "/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data/stomachCancer-germany.xls",
                                     path_to_output = "/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data/stomachCancer-germany.Rda")

