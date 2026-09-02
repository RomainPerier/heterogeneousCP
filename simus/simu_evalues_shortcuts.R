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

# --- Generate the pvalues and the list of cdfs ----------------------------------------------
results_list <- lapply(1:nb_repl, function(i){
  print(i)
  dat <- simu.one.array.fisher(N,m0_prime,m0_second,m1,p1,p2,v)
  test.result <- generate.pvalues(dat, "fisher")
  list_cdf <-test.result$get_pvalue_supports()
  pvals <-test.result$get_pvalues()
  compound_pvals <- compute.heterogeneous.compound.pvalues(pvals,list_cdf)
  res <- nb_rejected_hyp(pvals, alpha, list_cdf)
  data.frame(
    replicate = i,
    length_set_SC_cHBY = res$length_set_SC_cHBY,
    length_set_SC_cDBY = res$length_set_SC_cDBY,
    length_set_cBY  = closedBY(pvals, alpha = alpha),
    length_set_DBY  = length(BY_vanilla(compound_pvals,alpha))
  )
})
final_results_df <- do.call(rbind, results_list)
# --- Final data frame for rejection sets ------------------------------------------

save(final_results_df,file = file.path("resData/PBD/rejection_sets_length_SC.RData"))


