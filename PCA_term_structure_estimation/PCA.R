####################################################################################
# Application of PCA to term structure estimation
####################################################################################
#install.packages("magrittr")
#install.packages("dplyr")
library(dplyr)
library(magrittr)
library(muStat)
setwd("~/") # set this for the right directory
swissgov_yields <- read.csv(file="swissgov_yields.csv", header = TRUE) # read the file with the info.
colnames(swissgov_yields) <- c("Month","2Y","3Y","4Y","5Y","7Y","10Y","20Y","30Y")
monthly_changes <- sapply(swissgov_yields[2:9]/100, function(x) (x-lag(x))) %>% na.omit() %>% data.frame()

#pca <- prcomp(monthly_changes, center = TRUE, scale. = TRUE)
#summary(pca)
#loadings(pca)

# variance explanation by the first two PCs.
pca1 <- princomp(monthly_changes)
summary(pca1)
loadings(pca1)


# Find the SD of monthly change in portolio value
partial_yi <- function(i, cfi, yi){
  -i*cfi*exp(-i*yi) # partial derivative of portfolio V with respect to y_i
}

partial_t <- function(y2,y3,y4,y5,CF2,CF3,CF4,CF5){
  (-y2*CF2*exp(-2*y2) - y3*CF3*exp(-3*y3) - y4*CF4*exp(-4*y4) - y5*CF5*exp(-5*y5))
}

partial_t1 <- function(y2,y3,y4,y5,CF2,CF3,CF4,CF5){
  (+y2*CF2*exp(-2*y2) + y3*CF3*exp(-3*y3) + y4*CF4*exp(-4*y4) + y5*CF5*exp(-5*y5))
}


yields10jul <- swissgov_yields[120,2:9]/100

loading <- loadings(pca1)[1:4,1:2]
loading_full <- loadings(pca1)[,1:8]

mean_yield_change <- apply(monthly_changes,2,FUN = mean, na.rm = TRUE)
pc_values <- pca1$scores[119,1:2]
pc_values_full <- pca1$scores[119,]

delta_Y = list()

for (i in 1:nrow(loading)){
  delta_Y[i] = mean_yield_change[i] + (pc_values[1]*loading[i,1] + pc_values[2]*loading[i,2])
  # delta_Y[i] <- mean_yield_change[i] + (pc_values_full[1]*loading_full[i,1] + pc_values_full[2]*loading_full[i,2]
  #                                           + pc_values_full[3]*loading_full[i,3] + pc_values_full[4]*loading_full[i,4]
  #                                           + pc_values_full[5]*loading_full[i,5] + pc_values_full[6]*loading_full[i,6]
  #                                           + pc_values_full[7]*loading_full[i,7] + pc_values_full[8]*loading_full[i,8])
}

d_Y <- data.frame(matrix(unlist(delta_Y), nrow = length(delta_Y), byrow = T))
colnames(d_Y) <- c("DV") # 1 is 2Y

delta_V <- partial_t(yields10jul[1],yields10jul[2],yields10jul[3], yields10jul[4], 80,70,150,40) * (1/12) +
  (partial_yi(2,80,yields10jul[1]) * d_Y[1,1] + partial_yi(3,70,yields10jul[2])*d_Y[2,1] + partial_yi(4,150,yields10jul[3])*d_Y[3,1] + partial_yi(5,40,yields10jul[4])*d_Y[4,1])


delta_V1 <- partial_yi(2,80,yields10jul[1]) * delta_Y[1] + partial_yi(3,70,yields10jul[2]) * delta_Y[2] + partial_yi(4,150,yields10jul[3]) * delta_Y[3] + partial_yi(5,40,yields10jul[4]) * delta_Y[4]

delta_V2 <- partial_t1(yields10jul[1],yields10jul[2],yields10jul[3], yields10jul[4], 80,70,150,40) * (1/12) +
  (partial_yi(2,80,yields10jul[1]) * delta_Y[1] + partial_yi(3,70,yields10jul[2])*delta_Y[2] + partial_yi(4,150,yields10jul[3])*delta_Y[3] + partial_yi(5,40,yields10jul[4])*delta_Y[4])



colnames(delta_V) <- c("Portfolio Value 1m")
colnames(delta_V1) <- c("Portfolio Value 1m")
colnames(delta_V2) <- c("Portfolio Value 1m")



PVs_orginal <- list()
swiss_yields <- subset(swissgov_yields, select = -c(1))
swiss_yields <- swiss_yields/100
swiss_yields_only_rel <- subset(swiss_yields,select = -c(5,6,7,8))
for (i in 1:nrow(swiss_yields_only_rel)){
  PVs_orginal[i] <- 80*exp(-2*swiss_yields_only_rel[i,1])+70*exp(-3*swiss_yields_only_rel[i,2])+150*exp(-4*swiss_yields_only_rel[i,3])+40*exp(-5*swiss_yields_only_rel[i,4])
}
PV_orig <- data.frame(matrix(unlist(PVs_orginal), nrow = length(PVs_orginal), byrow = T))
colnames(PV_orig) <- c("Portfolio Value 1m")

PV_monthly_changes_orig <- sapply(PV_orig, function(x) (x-lag(x))) %>% na.omit() # this are the monthly changes on the orginal portfolio (no delta_V)
PV_monthly_changes_origV <- rbind(PV_monthly_changes_orig,delta_V)  # this is with negative derivative wrt time

PV_monthly_changes_origV1 <- rbind(PV_monthly_changes_orig, delta_V1) # this is WITHOUT derivative wrt time
PV_monthly_changes_origV2 <- rbind(PV_monthly_changes_orig, delta_V2) # this is with positive derivative wrt time

portfolio_monthly_chages_sd <-sapply(PV_monthly_changes_origV,sd)
pd <- data.frame(PV_monthly_changes_origV)
sd(pd$Portfolio.Value.1m)
stdev(pd$Portfolio.Value.1m,unbiased = FALSE)

portfolio_monthly_chages_sd1 <- sapply(PV_monthly_changes_origV1,sd)
pd_1 <- data.frame(PV_monthly_changes_origV1)
sd(pd_1$Portfolio.Value.1m)
stdev(pd_1$Portfolio.Value.1m,unbiased = FALSE)

portfolio_monthly_chages_sd2 <- sapply(PV_monthly_changes_origV2,sd)
pd_2 <- data.frame(PV_monthly_changes_origV2)
sd(pd_2$Portfolio.Value.1m)
stdev(pd_2$Portfolio.Value.1m,unbiased = FALSE)
