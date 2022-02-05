########################################################
# ZCB and FRA basic valuation
########################################################

rm(list=ls())

fwd_rates <- c(0.06,0.07,0.05,0.08)
CF <- 100
ZCB_PV <- round(CF/(prod(1+(fwd_rates/4))),2)
#ZCB_PV_continuos <- CF*exp(-prod(1+fwd_rates/4))
print(ZCB_PV)



P0_T0 <- 1/(1+0.06/4)
P0_T1 <- 1/(prod(1+fwd_rates/4))
N <- 100
k <- 0.04
delta <- 1-1/4 # (T1 - T0)
V_FRA <- N*(P0_T0 - P0_T1 - P0_T1*delta*k)
print(V_FRA)
