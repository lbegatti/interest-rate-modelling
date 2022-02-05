####################################################################
# Pricing of call and put using the Vasicek short-rate model.
####################################################################

T_1 <- 0.5    # given
T_0 <- 0.25   # given
R_0 <- 0.06   # given
k <- 0.86     # given
theta <- 0.08 # given
sigma <- 0.01 # given 
gamma <- sqrt(k^2+2*(sigma^2))



A_T1 <- -(2*k*theta/(sigma^2))*log((2*gamma*exp((gamma+k)*(T_1/2)))/((gamma+k)*(exp(gamma*T_1)-1)+ 2*gamma))


#AT1 <- B_T1*k*theta - (B_T1^2)*((sigma^2)/2)
B_T1 <- (1 - exp(-k*T_1))/k

A_T0 <- -(2*k*theta/(sigma^2))*log((2*gamma*exp((gamma+k)*(T_0/2)))/((gamma+k)*(exp(gamma*T_0)-1)+ 2*gamma))


#AT0 <- B_T0*k*theta - (B_T0^2)*((sigma^2)/2)
B_T0 <- (1 - exp(-k*T_0))/k

P_1 <- exp(-A_T1 - (B_T1*R_0))
P_0 <- exp(-A_T0 - (B_T0*R_0))
strike <- P_1/P_0           ## ATM strike price

d1 <- log(P_1/(strike*P_0)) + 0.5*((sigma^2)/(k^2))*(exp(-k*T_0)-exp(-k*T_1))^2*(exp(2*k*T_0)-1)/(2*k)/
  sqrt((sigma^2)/(k^2)*(exp(-k*T_0)-exp(-k*T_1))^2*(exp(2*k*T_0)-1)/(2*k))

d2 <- log(P_1/(strike*P_0)) - 0.5*((sigma^2)/(k^2))*(exp(-k*T_0)-exp(-k*T_1))^2*(exp(2*k*T_0)-1)/(2*k)/
  sqrt((sigma^2)/(k^2)*(exp(-k*T_0)-exp(-k*T_1))^2*(exp(2*k*T_0)-1)/(2*k))

price_put <- P_1*pnorm(-d2) - P_1*pnorm(-d1)
round(price_put*10000,2)

price_call <- P_1*pnorm(d1) - P_1*pnorm(d2)
