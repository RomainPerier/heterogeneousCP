# Exact comparison between raw DBY and the exact closed DBY procedure.
# This script is separate from the shortcut-based approximation and is meant to
# directly check the statement from the theorem: raw DBY is contained in the
# exact closed DBY procedure.

project_root <- normalizePath(getwd())
if (!file.exists(file.path(project_root, "scripts", "calibration_functions.R"))) {
  project_root <- normalizePath(file.path(project_root, ".."))
}

library(DiscreteFDR)
library(eClosure)

source(file.path(project_root, "scripts", "discrete_functions.R"))
source(file.path(project_root, "scripts", "calibration_functions.R"))
source(file.path(project_root, "config_files", "config_shortcuts.R"))

# -------------------------------
# Exact DBY comparison utilities
# -------------------------------

BY_vanilla <- function(pvals, alpha){
  m <- length(pvals)
  hm <- sum(1 / (1:m))
  thresholds <- (1:m) * alpha / (m * hm)
  sorted_pvals <- sort(pvals)
  aux <- which(sorted_pvals <= thresholds)
  if (length(aux) > 0) {
    k_star <- aux[length(aux)]
    threshold_final <- alpha * k_star / (m * hm)
  } else {
    threshold_final <- 0
  }
  which(pvals <= threshold_final)
}

compute_heterogeneous_compound_pvalues <- function(pvals, list_cdf){
  sapply(pvals, function(p){
    sum(sapply(list_cdf, function(cdf){ ev.cdf(p, cdf) })) / length(pvals)
  })
}

# Exact local DBY e-value for subset S, matching the paper's definition.
local_evalue_dby_subset <- function(pvals, list_cdf, subset_idx, alpha){
  S <- subset_idx
  m <- length(pvals)
  s <- length(S)
  hs_m <- sum(1 / (1:m))
  cpvals <- compute_heterogeneous_compound_pvalues(pvals[S], list_cdf[S])
  e_vals <- sapply(cpvals, function(x){
    e_calibrator(x, m, alpha * m / s, hs_m)
  })
  mean(e_vals)
}

# Check whether a set R belongs to the exact closure family under the DBY local e-values.
# This is the exact closure membership condition from the paper.
member_of_exact_closed_dby <- function(R, pvals, list_cdf, alpha){
  m <- length(pvals)
  if (length(R) == 0) return(TRUE)

  all_subsets <- unlist(lapply(1:m, function(s) {
    combn(1:m, s, simplify = FALSE)
  }), recursive = FALSE)

  for (S in all_subsets) {
    eS <- local_evalue_dby_subset(pvals, list_cdf, S, alpha)
    lhs <- alpha * eS
    rhs <- length(intersect(R, S)) / max(1, length(R))
    if (lhs < rhs) {
      return(FALSE)
    }
  }
  TRUE
}

# Direct check of the theorem for a single simulation realization.
check_exact_dby_inclusion <- function(pvals, list_cdf, alpha){
  compound_pvals <- compute_heterogeneous_compound_pvalues(pvals, list_cdf)
  R_dby <- BY_vanilla(compound_pvals, alpha)
  in_exact_closed <- member_of_exact_closed_dby(R_dby, pvals, list_cdf, alpha)
  list(
    raw_dby_set = R_dby,
    raw_dby_length = length(R_dby),
    exact_closed_membership = in_exact_closed
  )
}

# -------------------------------
# Small simulation harness
# -------------------------------

run_exact_dby_check <- function(N = 20, m = 8, alpha = 0.05, nrep = 10,
                                pi0 = 0.5, q = 0.25, p1 = 0.01, p2 = 0.10){
  set.seed(1)
  out <- list()

  for (rep in 1:nrep) {
    m0_prime <- round(pi0 * m * 0.5)
    m0_second <- round((pi0 - pi0 * 0.5) * m)
    m1 <- m - m0_prime - m0_second
    v <- c(0.1, q)

    dat <- simu.one.array.fisher(N, m0_prime, m0_second, m1, p1, p2, v)
    test.result <- generate.pvalues(dat, "fisher")
    list_cdf <- test.result$get_pvalue_supports()
    pvals <- test.result$get_pvalues()

    check <- check_exact_dby_inclusion(pvals, list_cdf, alpha)
    out[[rep]] <- data.frame(
      rep = rep,
      raw_dby_length = check$raw_dby_length,
      exact_closed_membership = check$exact_closed_membership
    )
  }

  do.call(rbind, out)
}

if (sys.nframe() == 0) {
  print(run_exact_dby_check())
}
