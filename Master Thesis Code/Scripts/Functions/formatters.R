# script containing functions for formatting data
library(readxl)
library(tidyverse)
library(lubridate)

#' Formats raw population data in excel file to a dataframe on the correct 
#' format. Saves the dataframe to a specified location
#'
#'@param path_to_data ("<filename>.xlsx") location of raw data in xlsx file
#'@param path_to_output ("path/to/output/location") location of output data
#'@param filename ("<filename>") the name of the outout file, without file type suffix
format_population_data <- function(path_to_data, path_to_output="", filename="", save=TRUE){
  population <- read_excel(path_to_data, sheet=1,
                           col_names = c("age", "male", "female", "total"),
                           skip = 6)
  
  years <- population %>% select(age) %>% 
    filter(!is.na(as.Date(age, format="%d.%m.%Y", optional = TRUE))) 
  
  population <- population %>% slice(1:1848) %>% 
    mutate(year = rep(as.vector(years$age), each=88)) %>%
    filter(year != age) %>% filter(age != "Insgesamt") %>%
    mutate(age.number =  parse_number(age)) %>%
    mutate(age.number = replace(age.number, age == "unter 1 Jahr", 0)) %>%
    mutate(age.int = 5*(age.number%/%5)) %>%
    mutate(age.int = paste(age.int, "-", age.int + 4)) %>%
    mutate(age.int = replace(age.int, age.int == "85 - 89", "85 +")) %>%
    group_by(age.int, year) %>%
    summarize(total = sum(total), male = sum(male), female = sum(female)) %>%
    mutate(year = format(as.POSIXct(year, format="%d.%m.%Y"), format="%Y")) %>%
    filter(year < 2017)
  
  if(save){
    save(population, file=file.path(path_to_output, paste(filename, ".Rda", sep = "")))
    write.csv(population, file = file.path(path_to_output, paste(filename, ".csv", sep="")), row.names = FALSE)
  }
  return(population)
}

# run format_population_data

# population <- format_population_data(
#   "/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data/population-germany.xlsx",
#   "/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data",
#   "population-germany")

#' Formats raw cancer mortality data (made for lung or stomach cancer) 
#' 
#'@param path_to_data ("<filename>.xlsx") location of raw data in xlsx file
#'@param path_to_output ("path/to/output/location") location of output dataframe
#'@param filename ("<filename>") the name of the outout file, without file type suffix
#'@param population (data.frame) formatted population data
format_cancer_data <- function(path_to_data, population, path_to_output="", filename="", save=TRUE){
  cancer.data  <- read_excel(path_to_data) %>%
    rename(sex = "...1") %>% rename(age = "...2") %>%
    mutate(age = replace(age, age == "85", "85 +")) %>%
    pivot_longer(!c(sex,age), names_to="year", values_to="deaths") %>%
    mutate(sex = replace(sex, sex == "männlich", "male")) %>%
    mutate(sex = replace(sex, sex == "weiblich", "female")) %>%
    pivot_wider(names_from = sex, values_from = deaths) %>% 
    mutate(total = male + female) %>% mutate(Y = total) %>%
    mutate(t = as.integer(year)-1999) %>% 
    mutate(age.int = parse_number(age)) %>%
    mutate(x = parse_number(age)%/%5) %>% mutate(x.1 = x) %>% mutate(x.c = x) %>%
    mutate(xt = (x*(2016-1998) +t)) %>% mutate(t.1 = t) %>%
    mutate(cohort = 5 * (max(x) - x) + t) %>% mutate(c = cohort) %>%
    mutate(birth.year = as.integer(year) - as.integer(age.int)) %>%
    mutate(birth.year = str_c(birth.year - 5, birth.year, sep = " - ")) %>%
    left_join(population, by = c("year" = "year", "age" = "age.int"), suffix = c("", ".t")) %>%
    mutate("mortality rate" = total/total.t) %>% 
    mutate(E = total.t)
  
  if(save){
    save(cancer.data, file=file.path(path_to_output, paste(filename, ".Rda", sep = "")))
    write.csv(cancer.data, file = file.path(path_to_output, paste(filename, ".csv", sep="")), row.names = FALSE)
  }
  return(cancer.data)
}

#   ----   Run format cancer data for stomach cancer data   ----
# stomach.cancer <- format_cancer_data(
#   "/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data/stomachCancer-germany.xls",
#   "/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data",
#   "stomachCancer-germany",
#   population=population)


#   ----   Run format cancer data for lung cancer data   ----

# lung.cancer <- format_cancer_data(
#   "/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data/lungCancer-germany.xls",
#   "/Users/helen/Desktop/Masteroppgave/Masters-thesis/Master Thesis Code/Data",
#   "lungCancer-germany",
#   population=population)


