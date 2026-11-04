# Code the function g from the paper, for some alpha and m ------------------

e_calibrator <- function(x,s,alpha,hs){
  if (hs*x<=alpha){
    return(s/(alpha*(ceiling(s*hs*x/alpha))))
  }
 else{
   return(0)
 }
}

e_calibrator_DBY <- function(x,s,alpha,hm,m){
  if (hm*x*s/m<=alpha){
    return(s/(alpha*(ceiling(s*hm*x/alpha))))
  }
  else{
    return(0)
  }
}


# BY procedure ----------------

BY_vanilla <- function(pvals,alpha){
  m <- length(pvals)  
  hm <- sum(1/1:m)
  sorted_pvals <- sort(pvals)
  tresholds <- 1:m*alpha/(m*hm)
  aux <- which(sorted_pvals <= tresholds)
  if (length(aux)>0){
    k_star <- aux[length(aux)]
    threshold_final <- alpha*k_star/(m*hm)
  }
  else{
    k_star <- 0
    threshold_final <- 0
  }
  return(which(pvals <= threshold_final))
}

# Code the function to get the vector of compound evalues from compound pvalues

cevalues_from_cpvalues <- function(cpvalues,s,alpha,hs){
  return(sapply(cpvalues,function(x){e_calibrator(x,s,alpha,hs)}))
}


get.local.evalues.subset <- function(pvalues, list_cdf, alpha, i, m, subset_of_hyp,vect_local_evalue,hs_vect){
      cpvalues <- compute.heterogeneous.compound.pvalues(pvalues[subset_of_hyp],list_cdf[subset_of_hyp])
      evalue_HBY <- 1/i * sum(cevalues_from_cpvalues(cpvalues,i,alpha,hs_vect[i]))
      evalue_CBY <- 1/i * sum(vect_local_evalue[subset_of_hyp])
      evalue_DBY <- 1/i * sum(cevalues_from_cpvalues(cpvalues,m,alpha*m/i,hs_vect[m]))
      return(list(subset_of_hyp = subset_of_hyp, evalue_HBY = evalue_HBY, evalue_CBY = evalue_CBY, evalue_DBY = evalue_DBY))
      }

get.local.evalues.subset.DBY <- function(pvalues, list_cdf, alpha, i, m, subset_of_hyp,vect_local_evalue,hs_vect){
  cpvalues <- compute.heterogeneous.compound.pvalues(pvalues[subset_of_hyp],list_cdf[subset_of_hyp])
  evalue_DBY <- 1/i * sum(cevalues_from_cpvalues(cpvalues,m,alpha*m/i,hs_vect[m]))
  return(list(subset_of_hyp = subset_of_hyp, evalue_DBY = evalue_DBY))
}


# --- generate the list of compound evalues for all subsets of hypotheses. It gives a list of two elements: 
# the first codes the subsets, the second codes the local evalues 
# each is a list of list such that the i-th element corresponds to the subsets of size i. --------------------------
get.local.evalues <- function(pvalues, list_cdf, alpha, m){
  hs_vect <- cumsum(1/(1:m))
  list_of_evalues <- list()
  for (i in 1:m){
    vect_local_evalue <- sapply(pvalues, function(x){e_calibrator(x,i,alpha,hs_vect[i])})
    list_of_subsets <- combn(1:m,i,simplify = FALSE)
    list_of_evalues[[i]] <- lapply(list_of_subsets,function(subset_of_hyp){
      get.local.evalues.subset(pvalues, list_cdf, alpha, i, m, subset_of_hyp,vect_local_evalue,hs_vect)
      }
    )
  }
  return(list_of_evalues)
}

get.local.evalues.DBY <- function(pvalues, list_cdf, alpha, m){
  hs_vect <- cumsum(1/(1:m))
  list_of_evalues <- list()
  for (i in 1:m){
    vect_local_evalue <- sapply(pvalues, function(x){e_calibrator(x,i,alpha,hs_vect[i])})
    list_of_subsets <- combn(1:m,i,simplify = FALSE)
    list_of_evalues[[i]] <- lapply(list_of_subsets,function(subset_of_hyp){
    get.local.evalues.subset.DBY(pvalues, list_cdf, alpha, i, m, subset_of_hyp,vect_local_evalue,hs_vect)
    })
  }
  return(list_of_evalues)
}


# Auxiliary function that compute worst case compound pvalues for each value of a and b --------------------
compute_f_wc_aux_matrix <-  function(p,R, Rc, list_cdf){
  if (is.null(Rc) || length(Rc) == 0) {
    sum_Rc_vect <- c(0)
  } 
  else {
    ev_cdf_vect_Rc <- sort(unlist(sapply(list_cdf[Rc], function(A) ev.cdf(p, A))), decreasing = TRUE)
    sum_Rc_vect <- c(0, cumsum(ev_cdf_vect_Rc))
  }
  ev_cdf_vect_R <- sort(unlist(sapply((list_cdf[R]),function(A){ev.cdf(p,A)})),decreasing = TRUE)
  sum_R_vect <- cumsum(ev_cdf_vect_R)
  return(outer(sum_R_vect, sum_Rc_vect, '+'))
}

# Gives if R belongs to cDBY or to cDBH, returns a list ----------------------------
is_in_closure_set <-function(pvals,R, alpha,list_cdf){
  m <- length(pvals)
  hs_vect <- cumsum(1/(1:m))
  k <- length(R)
  l <- m-k
  order_R <- order(pvals[R], decreasing = TRUE)
  pvals_sorted_R <- pvals[R[order_R]] 
  if (l>0){
    Rc <- setdiff(1:m, R)
    order_Rc <- order(pvals[Rc], decreasing = TRUE)
    pvals_sorted_Rc <- pvals[Rc[order_Rc]]
  }
  else{
    Rc <- NULL
  }
  list_sums <- list()
  for (i in 1:k){
    list_sums[[i]] <- compute_f_wc_aux_matrix(pvals_sorted_R[i], R, Rc, list_cdf)
  }
  if (l>0){
    for (i in 1:l){
      list_sums[[k+i]] <- compute_f_wc_aux_matrix(pvals_sorted_Rc[i], R, Rc, list_cdf)
    }
  }
  is_in_cHBY <- TRUE
  is_in_cDBY <- TRUE
  a <- 1
  while (a <= k & (is_in_cHBY||is_in_cDBY)){
    b <- 0
    while (b <= l & (is_in_cHBY||is_in_cDBY)) {
      if (b==0){
        if (is_in_cHBY){
           wc_evalue_HBY <-(1/(a+b))* (
          sum(sapply(1:a,function(i){
            e_calibrator((1/(a+b))*(list_sums[[i]][a,(b+1)]),(a+b),alpha,hs_vect[(a+b)])
          })
          )
          )
        }
        if (is_in_cDBY){
          wc_evalue_DBY <-(1/(a+b))* (
            sum(sapply(1:a,function(i){
              e_calibrator((1/(a+b))*(list_sums[[i]][a,(b+1)]),m,((alpha*m)/(a+b)),hs_vect[m]) 
            })
            )
          )
        }
      }
      else{
        if (is_in_cHBY){
          wc_evalue_HBY <-(1/(a+b))* (
            sum(sapply(1:a,function(i){
              e_calibrator((1/(a+b))*(list_sums[[i]][a,(b+1)]),(a+b),alpha,hs_vect[(a+b)])
            }
            )) 
            + sum(sapply(1:b,function(i){
              e_calibrator((1/(a+b))*(list_sums[[(k+i)]][a,(b+1)]),(a+b),alpha,hs_vect[(a+b)])
            }
            ))
          )
        }
        if (is_in_cDBY){
          wc_evalue_DBY <-(1/(a+b))* (
            sum(sapply(1:a,function(i){
              e_calibrator((1/(a+b))*(list_sums[[i]][a,(b+1)]),m,((alpha*m)/(a+b)),hs_vect[m])
            }
            )) 
            + sum(sapply(1:b,function(i){
              e_calibrator((1/(a+b))*(list_sums[[(k+i)]][a,(b+1)]),m,((alpha*m)/(a+b)),hs_vect[m])
            }
            ))
          )
        }
      }
      if (wc_evalue_HBY*alpha < (a/k) && !isTRUE(all.equal(wc_evalue_HBY*alpha,a/k))){
        is_in_cHBY <- FALSE
      }
      if (wc_evalue_DBY*alpha < (a/k) && !isTRUE(all.equal(wc_evalue_DBY*alpha,a/k))){
        is_in_cDBY <- FALSE
      }
      b <- b+1
    }
    a <- a+1
  }
  return(list(is_in_cHBY = is_in_cHBY, is_in_cDBY = is_in_cDBY))
}

# Function that compute the bigger top-k set (for the superset relation) that belong to the heterogeneous closed sets -----------

nb_rejected_hyp <- function(pvals, alpha, list_cdf,H0 = NULL){
  m <- length(pvals)
  order_pvals <- order(pvals)
  nb_rejected_hyp_cHBY <- 0
  nb_rejected_hyp_cDBY <- 0
  FDP_SC_cHBY <- 0
  FDP_SC_cDBY <- 0
  for (i in 1:m){
    result <- is_in_closure_set(pvals,order_pvals[1:i], alpha,list_cdf)
    if (result$is_in_cHBY){
      nb_rejected_hyp_cHBY <- i 
    }
    if (result$is_in_cDBY){
      nb_rejected_hyp_cDBY <- i 
    }
  }
  if (!is.null(H0)){
    if (nb_rejected_hyp_cHBY>0){
      nb_fp <- length(intersect(order_pvals[1:nb_rejected_hyp_cHBY],H0))
      if(nb_fp>0){
        FDP_SC_cHBY <- nb_fp/nb_rejected_hyp_cHBY
      }
    }
    if (nb_rejected_hyp_cDBY>0){
      nb_fp <- length(intersect(order_pvals[1:nb_rejected_hyp_cDBY],H0))
      if(nb_fp>0){
        FDP_SC_cDBY <- nb_fp/nb_rejected_hyp_cDBY
      }
    }
    return(list(length_set_SC_cHBY = nb_rejected_hyp_cHBY, length_set_SC_cDBY = nb_rejected_hyp_cDBY, FDP_SC_cHBY = FDP_SC_cHBY, FDP_SC_cDBY = FDP_SC_cDBY))
  }
  else{
    return(list(length_set_SC_cHBY = nb_rejected_hyp_cHBY, length_set_SC_cDBY = nb_rejected_hyp_cDBY))
  }
}

