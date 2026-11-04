setwd("../scripts")
source("../scripts/discrete_functions.R")
source("../scripts/calibration_functions.R")
setwd("../config_files")
source("../config_files/config_shortcuts.R")
setwd("../simus")

library(DiscreteFDR)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(forcats)
library(eClosure)

# --- Set seed to 1 ------------------------------------------

set.seed(1)

q <- 0.4
v <- c(p3,q)

# --- Apply the simulations from the paper with nb_repl replications with the the parameter q = 0.4 ----------------------------------------------
results_list <- lapply(1:nb_repl, function(i){
  print(i)
  dat <- simu.one.array.fisher(N,m0_prime,m0_second,m1,p1,p2,v)
  H0 <- 1:(m0_prime+m0_second)
  test.result <- generate.pvalues(dat, "fisher")
  list_cdf <-test.result$get_pvalue_supports()
  pvals <-test.result$get_pvalues()
  order_pvals <- order(pvals)
  compound_pvals <- compute.heterogeneous.compound.pvalues(pvals,list_cdf)
  order_compound_pvals <- order(compound_pvals)
  res <- nb_rejected_hyp(pvals, alpha, list_cdf, H0 = H0)
  DBY_vanilla_set <- BY_vanilla(compound_pvals,alpha)
  length_set_DBY <- length(DBY_vanilla_set)
  length_set_cBY <- closedBY(pvals, alpha = alpha)
  FDP_DBY <- 0
  FDP_cBY <- 0
  if (length_set_cBY>0){
    nb_fp <- length(intersect(order_pvals[1:length_set_cBY],H0))
    if(nb_fp>0){
      FDP_cBY <- nb_fp/length_set_cBY
    }
  }
  if (length_set_DBY>0){
    nb_fp <- length(intersect(order_pvals[1:length_set_DBY],H0))
    if(nb_fp>0){
      FDP_DBY <- nb_fp/length_set_DBY
    }
  }
  #  if (length(DBY_vanilla_set)!=0){
  #    print(is_in_closure_set(pvals, DBY_vanilla_set, alpha, list_cdf)[[2]])
  #  }
  #  if (i==2){
  #    local_evalues <- get.local.evalues.DBY(pvals, list_cdf, alpha, m)
  #    exit_loop_DBY <- FALSE
  #    for (j in 1:m){
  #      for (y in local_evalues[[j]]){
  #        nb_elements_in_intersection <- length(intersect(DBY_vanilla_set,y$subset_of_hyp))
  #        if (y$evalue_DBY < (nb_elements_in_intersection/(alpha*j)) && !isTRUE(all.equal(y$evalue_DBY, nb_elements_in_intersection/(alpha*j)))){
  #          exit_loop_DBY <- TRUE
  #          print(pvals)
  #          print(list_cdf)
  #          print(y)
  #          print(nb_elements_in_intersection/(alpha*j))
  #          print(all.equal(y$evalue_DBY, nb_elements_in_intersection/(alpha*j)))
  #          print(DBY_vanilla_set)
  #        }
  #        if (exit_loop_DBY){
  #          break
  #        }
  #      }
  #      if (exit_loop_DBY){
  #        break
  #      }
  #    }
  #    if (!exit_loop_DBY){
  #      print(TRUE)
  #    }
  #    else{
  #      print(FALSE)
  #    }
  #  }
  data.frame(
    replicate = i,
    length_set_SC_cHBY = res$length_set_SC_cHBY,
    length_set_SC_cDBY = res$length_set_SC_cDBY,
    length_set_cBY  = length_set_cBY,
    length_set_DBY  = length_set_DBY,
    FDP_SC_cHBY = res$FDP_SC_cHBY,
    FDP_SC_cDBY = res$FDP_SC_cDBY,
    FDP_SC_cBY = FDP_cBY,
    FDP_SC_DBY = FDP_DBY
  )
})
final_results_df <- do.call(rbind, results_list)
# --- Final data frame for rejection sets ------------------------------------------

save(final_results_df,file = file.path("resData/PBD/rejection_sets_length_SC2.RData"))


