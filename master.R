#1 --- Generates data for figures 1, 2, 3 and 4 (the seed is fixed within the file) ------------------------------------------
source("simu_evalues_shortcuts.R") 

#2 --- Generates and saves the figures inside the /figures folder ------------------------------------------
setwd("../figures")
# --- Figure 1 and 2 (q=0.25) -------------------------------------------------------
source("script_figure.R") 
# --- Figure 3 and 4 (q=0.4) -------------------------------------------------------
source("script_figure_2.R") 
setwd("..")
rm(list = ls())
