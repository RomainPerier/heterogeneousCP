# Generates a realization of m tables 2x2
simu.one.array.fisher <- function(N,m1,m2,m3,p1,p2,v){
  X1 <- c(rbinom(m1,N,p1),rbinom(m2,N,p2),rbinom(m3,N,v[1]))
  X2 <- c(rbinom(m1,N,p1),rbinom(m2,N,p2),rbinom(m3,N,v[2]))
  vec_N <- rep(N, m1+m2+m3)
  Y1 <- vec_N - X1
  Y2 <- vec_N - X2
  df <- data.frame(X1, Y1, X2, Y2)
}

# Transforms continuous p-values to discrete p-values on A
continuous.to.discrete <- function(p,A)  {
  return(A[which(A-p>=0)[1]])
}

# Evaluate the cdf

ev.cdf <- function(p,A)  {
  if (is.na(p)) stop("The p-value is not defined")
  else if (p==1) {return(1)}
  else if (is.na(A[which(A-p<=0)[1]])) {return(0)}
  else {B <-A[which(A-p<=0)]
  return(B[length(B)])}
}

# Compute the vector of locally compound $p$-values

compute.heterogeneous.compound.pvalues <- function(pvalues,list_cdf){
  sapply(pvalues, function(p){sum(sapply(list_cdf,function(cdf){ev.cdf(p,cdf)}))/length(pvalues)})
}


