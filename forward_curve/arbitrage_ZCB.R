############################################################
# Example of arbitrage with forward-implied ZCB
############################################################

T1 <- 1
T2 <- 2
T3 <- 3
t0 <- 0

f0 <- 0.04
p <- c(exp(-f0*(T1-t0)),exp(-f0*(T2-t0)),exp(-f0*(T3-t0)))

t1 <- 1
f_w1 <- 0.06
f_w2 <- 0.02

p_w1 <- c(1,exp(-f_w1*(T2-t1)),exp(-f_w1*(T3-t1))) # this is becuase at t=1 the value of the ZC1 is one by definition
p_w2 <- c(1,exp(-f_w2*(T2-t1)), exp(-f_w2*(T3-t1)))

A <- matrix(c(p[1],p_w1[1],p_w2[1],
              p[2],p_w1[2],p_w2[2],
              p[3],p_w1[3],p_w2[3]),3,3)
det(A)
b <- c(0,1,1)

bond_units <- round(solve(A,b),0) # unique solution as the determinants is not EXACTLY zero but still very small
bond_units



f <- function(m) class(try(solve(m),silent=T))=="matrix"
# x <- matrix(rep(1,25),nc=5)          # singular
# y <- matrix(1+1e-10*rnorm(25),nc=5)  # very nearly singular matrix
# z <- 0.001*diag(1,5) 
f(A)
