alpha <- c(0.05,0.1)[1] # level alpha
N <- 50 # total number of "patient"
m <- 10 # total number of hypotheses
pi0_vect <- c(0.1,0.5,0.8) # proportion of null hypotheses
pi0_prime <- c(0.2,0.5,0.8)[2] # proportion of the null hypotheses such that beta = 0.01
q_vect <- c(0.15,0.25,0.4) # signal strength
p1 <- 0.01 # probability of null 
p2 <- 0.1
p3 <- 0.1
q <- q_vect[2]
pi0 <- pi0_vect[1]
m0_prime <- round(pi0*pi0_prime*m) # number of null hypotheses such that beta = 0.01
m0_second <- round((pi0-pi0*pi0_prime)*m) # number of null hypotheses such that beta = 0.1
m1 <- m-m0_prime-m0_second
v <- c(p3,q) # vector of the probabilities of the alternative hypotheses

#configs <- expand.grid(
#  N = N,
#  pi0 = pi0_vect,
#  pi0_prime = pi0_prime,
#  q = q_vect,
#  p1 = p1,
#  p2 =p2 ,
#  p3 =p3)


