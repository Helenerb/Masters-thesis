# script for running inlabru analyses with the real data. 

load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data/population-germany.Rda")

load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data/lungCancer-germany.Rda")
lung.cancer <- cancer.data

load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data/stomachCancer-germany.Rda")
stomach.cancer <- cancer.data

# working directory at <Master\ Thesis\ Code>
source("Scripts/Synthetic data/Inlabru analyses/inlabru_analyses.R")

res.stomach <- inlabru.rw2.cohort.2(stomach.cancer)

source("Scripts/Real data/plot_real_data.R")

plots.stomach <- plot.inlabru.real(
  res.stomach, stomach.cancer, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2")

plot.hypers.inlabru.real(
  res.stomach, stomach.cancer, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2")

# run for lung cancer
res.lung <- inlabru.rw2.cohort.2(lung.cancer, max_iter = 100)

plots.lung <- plot.inlabru.real(
  res.lung, lung.cancer, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_rw2")

# run with extraconstraints on gamma, stomach
res.stomach.extraconstr <- inlabru.rw2.cohort.2.gamma.extraconstr(stomach.cancer, real_data = TRUE)

plot.hypers.inlabru.real(
  res.stomach.extraconstr, stomach.cancer, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2_extraconstr")

plots.stomach.extraconstr <- plot.inlabru.real(
  res.stomach.extraconstr, stomach.cancer, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2_extraconstr")

# run with extraconstraints on gamma, lung
res.stomach.extraconstr <- inlabru.rw2.cohort.2.gamma.extraconstr(stomach.cancer, real_data = TRUE)

plot.hypers.inlabru.real(
  res.stomach.extraconstr, stomach.cancer, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2_extraconstr")




