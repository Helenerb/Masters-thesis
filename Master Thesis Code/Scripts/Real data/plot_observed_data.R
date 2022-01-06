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
  geom_line(aes(x = age.int, y = male.mr.mean, color = "Male", linetype = "Male")) + 
  geom_line(aes(x = age.int, y = female.mr.mean, color = "Female", linetype = "Female")) + 
  geom_line(aes(x = age.int, y = total.mr.mean, color = "All sexes", linetype = "All sexes")) + 
  scale_color_manual(name = "", values = palette ) +
  scale_linetype_manual(name = "", values = c(1,2,3)) + 
  theme_classic() + 
  labs(title = "Average mortality rate for lung cancer 1999 - 2016", x= "Age", y = "")

ggsave("mr_lung_by_age.pdf", p.mr.lung.by.age, path = output.path, height = 4, width = 6.5, dpi = "retina")
#save.figure(p.mr.lung.by.age, name="mr_lung_by_age", path=output.path)

p.mr.lung.by.period <- ggplot(data = lung.cancer.by.period) +
  geom_line(aes(x = year, y = male.mr.mean, color = "Male", linetype = "Male")) + 
  geom_line(aes(x = year, y = female.mr.mean, color = "Female", linetype = "Female")) + 
  geom_line(aes(x = year, y = total.mr.mean, color = "All sexes", linetype = "All sexes")) + 
  scale_color_manual(name = "", values = palette ) +
  scale_linetype_manual(name = "", values = c(1,2,3)) + 
  theme_classic() + 
  labs(title = "Average mortality rate for lung cancer, ages 0 - 85+", x= "Year", y = "")

ggsave("mr_lung_by_period.pdf", p.mr.lung.by.period, path = output.path, height = 4, width = 6.5, dpi = "retina")
#save.figure(p.mr.lung.by.period, name="mr_lung_by_period", path=output.path)

p.mr.lung.by.cohort.w.hist <- ggplot(data = lung.cancer.by.cohort) +
  geom_histogram(data = lung.cancer, aes(x = cohort, y = after_stat(density), color = "Cohort observations", fill = "Cohort observations", linetype = "Cohort observations"), alpha = 0.5) + 
  geom_line(aes(x = cohort, y = male.mr.mean, color = "Male", fill = "Male", linetype = "Male")) + 
  geom_line(aes(x = cohort, y = female.mr.mean, color = "Female", fill="Female", linetype = "Female")) + 
  geom_line(aes(x = cohort, y = total.mr.mean, color = "All sexes", fill = "All sexes", linetype = "All sexes")) + 
  scale_color_manual(name = "", values = palette ) +
  scale_fill_manual(name = "", values = palette ) +
  scale_linetype_manual(name = "", values = c(2,1,3,4)) + 
  theme_classic() + 
  labs(title = "Average mortality rate for lung cancer, ages 0 - 85+", x= "Cohort", y = "")

ggsave("mr_lung_by_cohort_w_hist.pdf", p.mr.lung.by.cohort.w.hist, path = output.path, height = 4, width = 6.5, dpi = "retina")

p.mr.lung.by.cohort <- ggplot(data = lung.cancer.by.cohort) +
  geom_line(aes(x = cohort, y = male.mr.mean, color = "Male", fill = "Male", linetype = "Male")) + 
  geom_line(aes(x = cohort, y = female.mr.mean, color = "Female", fill="Female", linetype = "Female")) + 
  geom_line(aes(x = cohort, y = total.mr.mean, color = "All sexes", fill = "All sexes", linetype = "All sexes")) + 
  scale_color_manual(name = "", values = palette ) +
  scale_fill_manual(name = "", values = palette ) +
  scale_linetype_manual(name = "", values = c(1,2,3)) + 
  theme_classic() + 
  labs(title = "Average mortality rate for lung cancer, ages 0 - 85+", x= "Cohort", y = "")

ggsave("mr_lung_by_cohort.pdf", p.mr.lung.by.cohort, path = output.path, height = 4, width = 6.5, dpi = "retina")

# stomach cancer

p.mr.stomach.by.age <- ggplot(data = stomach.cancer.by.age) +
  geom_line(aes(x = age.int, y = male.mr.mean, color = "Male", linetype = "Male")) + 
  geom_line(aes(x = age.int, y = female.mr.mean, color = "Female", linetype = "Female")) + 
  geom_line(aes(x = age.int, y = total.mr.mean, color = "All sexes", linetype = "All sexes")) + 
  scale_color_manual(name = "", values = palette ) +
  scale_linetype_manual(name = "", values = c(1,2,3)) + 
  theme_classic() + 
  labs(title = "Average mortality rate for stomach cancer 1999 - 2016", x= "Age", y = "")

ggsave("mr_stomach_by_age.pdf", p.mr.stomach.by.age, path = output.path, height = 4, width = 6.5, dpi = "retina")

p.mr.stomach.by.period <- ggplot(data = stomach.cancer.by.period) +
  geom_line(aes(x = year, y = male.mr.mean, color = "Male", linetype = "Male")) + 
  geom_line(aes(x = year, y = female.mr.mean, color = "Female", linetype = "Female")) + 
  geom_line(aes(x = year, y = total.mr.mean, color = "All sexes", linetype = "All sexes")) + 
  scale_color_manual(name = "", values = palette ) +
  scale_linetype_manual(name = "", values = c(1,2,3)) + 
  theme_classic() + 
  labs(title = "Average mortality rate for stomach cancer, ages 0 - 85+", x= "Year", y = "")

ggsave("mr_stomach_by_period.pdf", p.mr.stomach.by.period, path = output.path, height = 4, width = 6.5, dpi = "retina")


p.mr.stomach.by.cohort.w.hist <- ggplot(data = stomach.cancer.by.cohort) +
  geom_histogram(data = lung.cancer, aes(x = cohort, y = after_stat(density), color = "Cohort observations", fill = "Cohort observations", linetype = "Cohort observations"), alpha = 0.5) + 
  geom_line(aes(x = cohort, y = male.mr.mean, color = "Male", fill = "Male", linetype = "Male")) + 
  geom_line(aes(x = cohort, y = female.mr.mean, color = "Female", fill="Female", linetype = "Female")) + 
  geom_line(aes(x = cohort, y = total.mr.mean, color = "All sexes", fill = "All sexes", linetype = "All sexes")) + 
  scale_color_manual(name = "", values = palette ) +
  scale_fill_manual(name = "", values = palette ) +
  scale_linetype_manual(name = "", values = c(2,1,3,4)) + 
  theme_classic() + 
  labs(title = "Average mortality rate for stomach cancer, ages 0 - 85+", x= "Cohort", y = "")

ggsave("mr_stomach_by_cohort_w_hist.pdf", p.mr.stomach.by.cohort.w.hist, path = output.path, height = 4, width = 6.5, dpi = "retina")

p.mr.stomach.by.cohort <- ggplot(data = stomach.cancer.by.cohort) +
  geom_line(aes(x = cohort, y = male.mr.mean, color = "Male", fill = "Male", linetype = "Male")) + 
  geom_line(aes(x = cohort, y = female.mr.mean, color = "Female", fill="Female", linetype = "Female")) + 
  geom_line(aes(x = cohort, y = total.mr.mean, color = "All sexes", fill = "All sexes", linetype = "All sexes")) + 
  scale_color_manual(name = "", values = palette ) +
  scale_fill_manual(name = "", values = palette ) +
  scale_linetype_manual(name = "", values = c(1,2,3)) + 
  theme_classic() + 
  labs(title = "Average mortality rate for stomach cancer, ages 0 - 85+", x= "Cohort", y = "")

ggsave("mr_stomach_by_cohort.pdf", p.mr.stomach.by.cohort, path = output.path, height = 4, width = 6.5, dpi = "retina")

###   ----   Plot absoulute observed data   ----

p.all.total.by.year <- ggplot(data = lung.cancer.by.period) + 
  geom_point(aes(x = year, y = male.t.all, color = "Male", shape = "Male")) + 
  geom_line(aes(x = year, y = male.t.all, color = "Male")) + 
  geom_point(aes(x = year, y = female.t.all, color = "Female", shape = "Female")) + 
  geom_line(aes(x = year, y = female.t.all, color = "Female")) + 
  scale_color_manual(name = "", values = palette) +
  scale_shape_manual(name = "", values = c(16, 17)) + 
  theme_classic() + 
  labs(title = "Total male and female population, years 1999 - 2016", x= "Year", y = "")

ggsave("all_total_by_year.pdf", p.all.total.by.year, path = output.path, height = 4, width = 6.5, dpi = "retina")

p.all.lung.by.year <- ggplot(data = lung.cancer.by.period) + 
  geom_point(aes(x = year, y = male.all, color = "Male", shape = "Male")) + 
  geom_line(aes(x = year, y = male.all, color = "Male")) + 
  geom_point(aes(x = year, y = female.all, color = "Female", shape = "Female")) + 
  geom_line(aes(x = year, y = female.all, color = "Female")) + 
  scale_color_manual(name = "", values = palette) +
  scale_shape_manual(name = "", values = c(16, 17)) + 
  theme_classic() + 
  labs(title = "Male and female lung cancer deaths, years 1999 - 2016", x= "Year", y = "")

ggsave("all_lung_by_year.pdf", p.all.lung.by.year, path = output.path, height = 4, width = 6.5, dpi = "retina")

p.all.stomach.by.year <- ggplot(data = stomach.cancer.by.period) + 
  geom_point(aes(x = year, y = male.all, color = "Male", shape = "Male")) + 
  geom_line(aes(x = year, y = male.all, color = "Male")) + 
  geom_point(aes(x = year, y = female.all, color = "Female", shape = "Female")) + 
  geom_line(aes(x = year, y = female.all, color = "Female")) + 
  scale_color_manual(name = "", values = palette) +
  scale_shape_manual(name = "", values = c(16, 17)) + 
  theme_classic() + 
  labs(title = "Male and female stomach cancer deaths, years 1999 - 2016", x= "Year", y = "")

ggsave("all_stomach_by_year.pdf", p.all.stomach.by.year, path = output.path, height = 4, width = 6.5, dpi = "retina")


###    ----    Plot total observed population and deaths for each age, for years 1999, 2007 and 2016    -----

stomach.cancer.3.years <- stomach.cancer %>% filter(year %in% c(1999, 2007, 2016))
lung.cancer.3.years <- lung.cancer %>% filter(year %in% c(1999, 2007, 2016))

p.3.year.stomach.by.age <- ggplot(data = stomach.cancer.3.years) + 
  geom_point(aes(x = age.int, y = male, color = as.factor(year), shape = "Male")) + 
  geom_point(aes(x = age.int, y = female, color = as.factor(year), shape = "Female")) + 
  geom_line(aes(x = age.int, y = male, color = as.factor(year), linetype = as.factor(year))) + 
  geom_line(aes(x = age.int, y = female, color = as.factor(year), linetype = as.factor(year))) + 
  scale_color_manual(name = "", values = palette) +
  scale_shape_manual(name = "", values = c(1, 2)) + 
  scale_linetype_manual(name = "", values= c(1,5,3)) + 
  theme_classic() + 
  labs(title = "Male and female stomach cancer deaths, ages 0 - 85+", x= "Year", y = "")

ggsave("all_stomach_by_age.pdf", p.3.year.stomach.by.age, path = output.path, height = 4, width = 6.5, dpi = "retina")

p.3.year.lung.by.age <- ggplot(data = lung.cancer.3.years) + 
  geom_point(aes(x = age.int, y = male, color = as.factor(year), shape = "Male")) + 
  geom_point(aes(x = age.int, y = female, color = as.factor(year), shape = "Female")) + 
  geom_line(aes(x = age.int, y = male, color = as.factor(year), linetype = as.factor(year))) + 
  geom_line(aes(x = age.int, y = female, color = as.factor(year), linetype = as.factor(year))) + 
  scale_color_manual(name = "", values = palette) +
  scale_shape_manual(name = "", values = c(1, 2)) + 
  scale_linetype_manual(name = "", values= c(1,5,3)) + 
  theme_classic() + 
  labs(title = "Male and female lung cancer deaths, ages 0 - 85+", x= "Year", y = "")

ggsave("all_lung_by_age.pdf", p.3.year.lung.by.age, path = output.path, height = 4, width = 6.5, dpi = "retina")

p.3.year.total.by.age <- ggplot(data = lung.cancer.3.years) + 
  geom_point(aes(x = age.int, y = male.t, color = as.factor(year), shape = "Male")) + 
  geom_point(aes(x = age.int, y = female.t, color = as.factor(year), shape = "Female")) + 
  geom_line(aes(x = age.int, y = male.t, color = as.factor(year), linetype = as.factor(year))) + 
  geom_line(aes(x = age.int, y = female.t, color = as.factor(year), linetype = as.factor(year))) + 
  scale_color_manual(name = "", values = palette) +
  scale_shape_manual(name = "", values = c(1, 2)) + 
  scale_linetype_manual(name = "", values= c(1,5,3)) + 
  theme_classic() + 
  labs(title = "Male and female German population, ages 0 - 85+", x= "Year", y = "")

ggsave("all_total_by_age.pdf", p.3.year.total.by.age, path = output.path, height = 4, width = 6.5, dpi = "retina")

