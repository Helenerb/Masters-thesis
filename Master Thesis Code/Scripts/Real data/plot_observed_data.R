# plots observed cancer and population data
library("ggplot2")
library("tidyverse")
library("patchwork")

source("Scripts/Misc/palette.R")
source("Scripts/Functions/plotters.R")

load("Data/population-germany.Rda")

load("Data/lungCancer-germany.Rda")
lung.cancer <- cancer.data %>% mutate(female.mr = female/female.t, male.mr = male/male.t) %>%
  mutate(year = parse_integer(year))

load("Data/stomachCancer-germany.Rda")
stomach.cancer <- cancer.data %>% mutate(female.mr = female/female.t, male.mr = male/male.t) %>%
  mutate(year = parse_integer(year))

output.path <- "Scripts/Real data/Output/Figures/observed_data"

#   ----   create summaries of dataframes   ----:


stomach.cancer.by.age <- stomach.cancer %>%
  group_by(age.int, age) %>%
  summarise(Y.mean = mean(Y), Y.all = sum(Y),
            male.mean = mean(male), male.all = sum(male),
            female.mean = mean(female), female.all = sum(female),
            total.t.mean = mean(total.t), total.t.all = sum(total.t),
            male.t.mean = mean(male.t), male.t.all = sum(male.t),
            female.t.mean = mean(female.t), female.t.all = sum(female.t),
            total.mr.mean = mean(`mortality rate`), 
            male.mr.mean = mean(male.mr),
            female.mr.mean = mean(female.mr)
            )

stomach.cancer.by.period <- stomach.cancer %>%
  group_by(year) %>%
  summarise(Y.mean = mean(Y), Y.all = sum(Y),
            male.mean = mean(male), male.all = sum(male),
            female.mean = mean(female), female.all = sum(female),
            total.t.mean = mean(total.t), total.t.all = sum(total.t),
            male.t.mean = mean(male.t), male.t.all = sum(male.t),
            female.t.mean = mean(female.t), female.t.all = sum(female.t),
            total.mr.mean = mean(`mortality rate`), 
            male.mr.mean = mean(male.mr),
            female.mr.mean = mean(female.mr)
  )

stomach.cancer.by.cohort <- stomach.cancer %>%
  group_by(cohort, birth.year) %>%
  summarise(Y.mean = mean(Y), Y.all = sum(Y),
            male.mean = mean(male), male.all = sum(male),
            female.mean = mean(female), female.all = sum(female),
            total.t.mean = mean(total.t), total.t.all = sum(total.t),
            male.t.mean = mean(male.t), male.t.all = sum(male.t),
            female.t.mean = mean(female.t), female.t.all = sum(female.t),
            total.mr.mean = mean(`mortality rate`), 
            male.mr.mean = mean(male.mr),
            female.mr.mean = mean(female.mr)
  )

lung.cancer.by.age <- lung.cancer %>%
  group_by(age.int, age) %>%
  summarise(Y.mean = mean(Y), Y.all = sum(Y),
            male.mean = mean(male), male.all = sum(male),
            female.mean = mean(female), female.all = sum(female),
            total.t.mean = mean(total.t), total.t.all = sum(total.t),
            male.t.mean = mean(male.t), male.t.all = sum(male.t),
            female.t.mean = mean(female.t), female.t.all = sum(female.t),
            total.mr.mean = mean(`mortality rate`), 
            male.mr.mean = mean(male.mr),
            female.mr.mean = mean(female.mr)
  )

lung.cancer.by.period <- lung.cancer %>%
  group_by(year) %>%
  summarise(Y.mean = mean(Y), Y.all = sum(Y),
            male.mean = mean(male), male.all = sum(male),
            female.mean = mean(female), female.all = sum(female),
            total.t.mean = mean(total.t), total.t.all = sum(total.t),
            male.t.mean = mean(male.t), male.t.all = sum(male.t),
            female.t.mean = mean(female.t), female.t.all = sum(female.t),
            total.mr.mean = mean(`mortality rate`), 
            male.mr.mean = mean(male.mr),
            female.mr.mean = mean(female.mr)
  )

lung.cancer.by.cohort <- lung.cancer %>%
  group_by(cohort, birth.year) %>%
  summarise(Y.mean = mean(Y), Y.all = sum(Y),
            male.mean = mean(male), male.all = sum(male),
            female.mean = mean(female), female.all = sum(female),
            total.t.mean = mean(total.t), total.t.all = sum(total.t),
            male.t.mean = mean(male.t), male.t.all = sum(male.t),
            female.t.mean = mean(female.t), female.t.all = sum(female.t),
            total.mr.mean = mean(`mortality rate`), 
            male.mr.mean = mean(male.mr),
            female.mr.mean = mean(female.mr)
  )


#   ----   save plots   ----

p.mr.lung.by.age <- ggplot(data = lung.cancer.by.age) +
  geom_line(aes(x = age.int, y = male.mr.mean, color = "Male")) + 
  geom_line(aes(x = age.int, y = female.mr.mean, color = "Female")) + 
  geom_line(aes(x = age.int, y = total.mr.mean, color = "All sexes")) + 
  scale_color_manual(name = "", values = palette ) +
  theme_classic() + 
  labs(title = "Average mortality rate for lung cancer 1999 - 2016", x= "age", y = "")

save.figure(p.mr.lung.by.age, name="mr_lung_by_age", path=output.path)


p.mr.lung.by.period <- ggplot(data = lung.cancer.by.period) +
  geom_line(aes(x = year, y = male.mr.mean, color = "Male")) + 
  geom_line(aes(x = year, y = female.mr.mean, color = "Female")) + 
  geom_line(aes(x = year, y = total.mr.mean, color = "All sexes")) + 
  scale_color_manual(name = "", values = palette ) +
  theme_classic() + 
  labs(title = "Average mortality rate for lung cancer, ages 0 - 85+", x= "year", y = "")
p.mr.lung.by.period

save.figure(p.mr.lung.by.period, name="mr_lung_by_period", path=output.path)


p.mr.lung.by.cohort <- ggplot(data = lung.cancer.by.cohort) +
  geom_histogram(data = lung.cancer, aes(x = cohort, y = after_stat(density), color = "Cohort observations", fill = "Cohort observations"), alpha = 0.5) + 
  geom_line(aes(x = cohort, y = male.mr.mean, color = "Male", fill = "Male")) + 
  geom_line(aes(x = cohort, y = female.mr.mean, color = "Female", fill="Female")) + 
  geom_line(aes(x = cohort, y = total.mr.mean, color = "All sexes", fill = "All sexes")) + 
  scale_color_manual(name = "", values = palette ) +
  scale_fill_manual(name = "", values = palette ) +
  theme_classic() + 
  labs(title = "Average mortality rate for lung cancer, ages 0 - 85+", x= "cohort", y = "")

save.figure(p.mr.lung.by.cohort, name="mr_lung_by_cohort", path=output.path)

# stomach cancer

p.mr.stomach.by.age <- ggplot(data = stomach.cancer.by.age) +
  geom_line(aes(x = age.int, y = male.mr.mean, color = "Male")) + 
  geom_line(aes(x = age.int, y = female.mr.mean, color = "Female")) + 
  geom_line(aes(x = age.int, y = total.mr.mean, color = "All sexes")) + 
  scale_color_manual(name = "", values = palette ) +
  theme_classic() + 
  labs(title = "Average mortality rate for stomach cancer 1999 - 2016", x= "age", y = "")

save.figure(p.mr.stomach.by.age, name="mr_stomach_by_age", path=output.path)

p.mr.stomach.by.period <- ggplot(data = stomach.cancer.by.period) +
  geom_line(aes(x = year, y = male.mr.mean, color = "Male")) + 
  geom_line(aes(x = year, y = female.mr.mean, color = "Female")) + 
  geom_line(aes(x = year, y = total.mr.mean, color = "All sexes")) + 
  scale_color_manual(name = "", values = palette ) +
  theme_classic() + 
  labs(title = "Average mortality rate for stomach cancer, ages 0 - 85+", x= "year", y = "")

save.figure(p.mr.stomach.by.period, name="mr_stomach_by_period", path=output.path)


p.mr.stomach.by.cohort <- ggplot(data = stomach.cancer.by.cohort) +
  geom_histogram(data = lung.cancer, aes(x = cohort, y = after_stat(density), color = "Cohort observations", fill = "Cohort observations"), alpha = 0.5) + 
  geom_line(aes(x = cohort, y = male.mr.mean, color = "Male", fill = "Male")) + 
  geom_line(aes(x = cohort, y = female.mr.mean, color = "Female", fill="Female")) + 
  geom_line(aes(x = cohort, y = total.mr.mean, color = "All sexes", fill = "All sexes")) + 
  scale_color_manual(name = "", values = palette ) +
  scale_fill_manual(name = "", values = palette ) +
  theme_classic() + 
  labs(title = "Average mortality rate for stomach cancer, ages 0 - 85+", x= "cohort")

save.figure(p.mr.stomach.by.cohort, name="mr_stomach_by_cohort", path=output.path)
