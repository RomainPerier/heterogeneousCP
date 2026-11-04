setwd("../scripts")
source("../scripts/discrete_functions.R")
source("../scripts/calibration_functions.R")
setwd("../config_files")
source("../config_files/config_full_closure.R")
setwd("../simus")

library(DiscreteFDR)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(forcats)

# --- Set seed to 1 ------------------------------------------

set.seed(1)

# --- Generate the pvalues and the list of cdfs ----------------------------------------------

dat <- simu.one.array.fisher(N,m0_prime,m0_second,m1,p1,p2,v)
test.result <- generate.pvalues(dat, "fisher")
list_cdf <-test.result$get_pvalue_supports()
pvalues <-test.result$get_pvalues()


# --- Get the rejected sets for each of the closure sets ------------------------------------------

local_evalues <- get.local.evalues(pvalues, list_cdf, alpha, m)
rejected_sets_HBY <- list()
rejected_sets_CBY <- list()
rejected_sets_DBY <- list()
cpt <- 1
for (i in 1:m){
  rejected_sets_HBY[[i]] <- list()
  rejected_sets_CBY[[i]] <- list()
  rejected_sets_DBY[[i]] <- list()
  a <- sapply(local_evalues[[i]], function(x){
    exit_loop_HBY <- FALSE
    exit_loop_CBY <- FALSE
    exit_loop_DBY <- FALSE
    for (j in 1:m){
      for (y in local_evalues[[j]]){
        nb_elements_in_intersection <- length(intersect(x$subset_of_hyp,y$subset_of_hyp))
        if (y$evalue_HBY < (nb_elements_in_intersection/(alpha*j))){
          exit_loop_HBY <- TRUE
        }
        if (y$evalue_CBY < (nb_elements_in_intersection/(alpha*j))){
          exit_loop_CBY <- TRUE
        }
        if (y$evalue_DBY < (nb_elements_in_intersection/(alpha*j))){
          exit_loop_DBY <- TRUE
        }
        if (exit_loop_HBY & exit_loop_CBY & exit_loop_DBY){
          break
        }
        if (exit_loop_HBY & exit_loop_CBY & exit_loop_DBY){
        break
        }
      }
      if (exit_loop_HBY & exit_loop_CBY & exit_loop_DBY){
        break
      }
    }
    if (!exit_loop_HBY){
      rejected_sets_HBY[[i]][[length(rejected_sets_HBY[[i]])+1]] <<- x$subset_of_hyp
    }
    if (!exit_loop_CBY){
      rejected_sets_CBY[[i]][[length(rejected_sets_CBY[[i]])+1]] <<- x$subset_of_hyp
    }
    if (!exit_loop_DBY){
      rejected_sets_DBY[[i]][[length(rejected_sets_DBY[[i]])+1]] <<- x$subset_of_hyp
    }
  }
  )
}


result_closedBY <- closedBY(pvalues,alpha = alpha)

# --- Final data frame for rejection sets ------------------------------------------

results_df <- data.frame(
  pi0 = pi0,
  q = q,
  m0_prime = m0_prime,
  m0_second = m0_second,
  m1 = m1,
  rejected_sets_HBY = rejected_sets_HBY,
  rejected_sets_CBY = rejected_sets_CBY,
  rejected_sets_DBY = rejected_sets_DBY
)

save(results_df,file = file.path("resData/PBD/rejection_sets.RData"))

rm(list = ls())
