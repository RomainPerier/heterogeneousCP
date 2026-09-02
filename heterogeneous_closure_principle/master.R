#1 --- Generates data for figure 3, 5 and 7 (the seed is fixed within the file) ------------------------------------------
source("simu_evalues_shortcuts.R") 

#2 --- Save the figures inside the /figures folder ------------------------------------------
setwd("../figures")
source("script_figure.R") 
setwd("..")
rm(list = ls())
