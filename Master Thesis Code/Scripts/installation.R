# script for installation of libraries
.libPaths("~/Documents/R_libraries")

#install.packages("ggplot2")
#install.packages("tidyverse")

#  STAN
#Sys.setenv(DOWNLOAD_STATIC_LIBV8 = 1) # only necessary for Linux without the nodejs library / headers
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)

# inla
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

# inlabru
#install.packages("inlabru")

# patchwork - for plotting
#install.packages("patchwork")

#    ----   Second version   ----

# re-install inlabru:
remove.packages("inlabru")
install.packages("remotes")
remotes::install_github("inlabru-org/inlabru", ref="stable")
