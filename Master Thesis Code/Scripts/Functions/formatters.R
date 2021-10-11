# script containing functions for formatting data

#' Formats raw population data in excel file to a dataframe on the correct 
#' format. Saves the dataframe to a specified location
#'
#'@param path_to_data ("<filename>.xlsx") location of raw data in xlsx file
#'@param path_to_output ("<filename>.Rda") location of output dataframe
format_population_data <- function(path_to_data, path_to_output){
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
    filter(year < 2017) %>%
    pivot_longer(c(total, male, female), names_to = "sex", values_to="population")
  
  save(population, file=path_to_output)
  return(population)
}

#' Formats raw cancer mortality data (made for lung or stomach cancer) 
#' 
#'@param path_to_data ("<filename>.xlsx") location of raw data in xlsx file
#'@param path_to_output ("<filename>.Rda") location of output dataframe
format_cancer_data <- function(path_to_data, path_to_output){
  cancer.data  <- read_excel(path_to_data) %>%
    rename(sex = "...1") %>% rename(age = "...2") %>%
    mutate(age = replace(age, age == "85", "85 +")) %>%
    pivot_longer(!c(sex,age), names_to="year", values_to="deaths") %>%
    mutate(sex = replace(sex, sex == "männlich", "male")) %>%
    mutate(sex = replace(sex, sex == "weiblich", "female")) %>%
    pivot_wider(names_from = sex, values_from = deaths) %>% 
    mutate(total = male + female) %>%
    pivot_longer(c(total, male, female), names_to = "sex", values_to = "cases") %>%
    mutate(t = as.integer(year)-1999) %>% 
    mutate(age.int = parse_number(age)) %>%
    mutate(x = parse_number(age)%/%5) %>% 
    mutate(k = 5 * (max(x) - x) + t) %>%
    mutate(birth.year = as.integer(year) - as.integer(age.int)) %>%
    mutate(birth.year = str_c(birth.year - 5, birth.year, sep = " - ")) %>%
    left_join(population, by = c("year" = "year", "age" = "age.int", "sex" = "sex"), suffix = c("", ".t")) %>%
    filter(sex != "total") %>%
    mutate(s = recode(sex, male = 0, female = 1))  %>%
    mutate(xts = (x*(2016-1998) + t + s*(max(x)*(2016-1998) + max(t)))) %>% 
    mutate("mortality rate" = cases/population)
  
  save(cancer, file=path_to_output)
  return(cancer)
}

