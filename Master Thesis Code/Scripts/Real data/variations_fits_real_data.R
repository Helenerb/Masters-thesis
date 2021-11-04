# test different variations of fits to the real data

load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data/population-germany.Rda")

load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data/lungCancer-germany.Rda")
lung.cancer <- cancer.data

load("/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data/stomachCancer-germany.Rda")
stomach.cancer <- cancer.data

# working directory at <Master\ Thesis\ Code>
source("Scripts/Synthetic data/Inlabru analyses/inlabru_analyses.R")

# run analyisis for male and female cancer mortality separately:

#   ----   reformat data   ----

# calculate male and female mortality:

lung.cancer <- lung.cancer %>%
  mutate(female.mr = female/female.t, male.mr = male/male.t)

stomach.cancer <- stomach.cancer %>%
  mutate(female.mr = female/female.t, male.mr = male/male.t) 
# separate dfs for male and female mortality:

lung.cancer.female <- lung.cancer %>%
  select(x, x.c, t, xt, cohort, c, age, age.int, year, birth.year, female.t, female, female.mr) %>%
  mutate(E = female.t, Y = female, `mortality rate` = female.mr)

stomach.cancer.female <- stomach.cancer %>%
  select(x, x.c, t, xt, cohort, c, age, age.int, year, birth.year, female.t, female, female.mr) %>%
  mutate(E = female.t, Y = female, `mortality rate` = female.mr)

lung.cancer.male <- lung.cancer %>%
  select(x, x.c, t, xt, cohort, c, age, age.int, year, birth.year, male.t, male, male.mr) %>%
  mutate(E = male.t, Y = male, `mortality rate` = male.mr)

stomach.cancer.male <- stomach.cancer %>%
  select(x, x.c, t, xt, cohort, c, age, age.int, year, birth.year, male.t, male, male.mr) %>%
  mutate(E = male.t, Y = male, `mortality rate` = male.mr)

#   ----   perform inlabru analyses   ----
source("Scripts/Synthetic data/Inlabru analyses/inlabru_analyses.R")

res.lung.lc.f <- inlabru.rw2.lc.2(lung.cancer.female)

res.stomach.lc.f <- inlabru.rw2.lc.2(stomach.cancer.female)

res.lung.lc.m <- inlabru.rw2.lc.2(lung.cancer.male)

res.stomach.lc.m <- inlabru.rw2.lc.2(stomach.cancer.male)

#   ----   plot and save results   ----

source("Scripts/Real data/plot_real_data.R")

# female lung
plots.lung.female <- plot.inlabru.real(
  res.lung.lc.f, lung.cancer.female, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_rw2_lc/female",
  cohort=FALSE)

plot.hypers.inlabru.real(
  res.lung.lc.f, lung.cancer.female, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_rw2_lc/female",
  cohort=FALSE)

# female stomach
plots.stomach.female <- plot.inlabru.real(
  res.stomach.lc.f, stomach.cancer.female, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2_lc/female",
  cohort=FALSE)

plot.hypers.inlabru.real(
  res.stomach.lc.f, stomach.cancer.female, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2_lc/female")

# male lung
plots.lung.male <- plot.inlabru.real(
  res.lung.lc.m, lung.cancer.male, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_rw2_lc/male",
  cohort=FALSE)

plot.hypers.inlabru.real(
  res.lung.lc.m, lung.cancer.male, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/lung_rw2_lc/male")

# male stomach
plots.stomach.male <- plot.inlabru.real(
  res.stomach.lc.m, stomach.cancer.male, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2_lc/male",
  cohort=FALSE)

plot.hypers.inlabru.real(
  res.stomach.lc.m, stomach.cancer.male, save=TRUE, 
  path.to.storage = "Scripts/Real data/Output/Figures/stomach_rw2_lc/male")
