########################################################################################
# Code for duration/convexity hedging and respective price change of the bond portfolio.
########################################################################################

rm(list=ls())

labels = c("Portfolio", "Bond1", "Bond2")

yields = c(6, 5.8, 5.62, 5.46, 5.33)/100
times = 1:5
disc = exp(-yields*times)

CF = matrix(c(6,8,106,7,9, 4,104,0,0,0, 10,10,10,110,0), 5, 3)

prices = apply(CF * disc, 2, sum)
duration = apply(CF * disc * times, 2, sum)/prices

convexity = apply(CF * disc * (times^2), 2, sum)/prices


#duration hedge
q = (- duration[1]*prices[1])/(duration[3]*prices[3])
print(q) ### This answer is accepted in the grading

#rel price change duration hedge portfolio
yieldsShock = yields + 0.02
discShock = exp(-yieldsShock * times)
pricesShock = apply(CF * discShock, 2, sum)
originalPrice = prices[1] + q * prices[3]
shockPrice = pricesShock[1] + q * pricesShock[3]
print(10^4 * (shockPrice-originalPrice)/originalPrice)


#convexity hedge
LHS <- matrix(c(-duration[2]*prices[2],convexity[2]*prices[2], -duration[3]*prices[3], convexity[3]*prices[3]),2,2)
RHS <- c(duration[1]*prices[1],-convexity[1]*prices[1])
q1_q2 <- solve(LHS,RHS)
print(q1_q2[2])
# A <- matrix(c(1, 2, -1, 2), 2, 2)
# b <- c(2,1)
# solve(A,b, fractions = FALSE)


# rel price change convexity hedge portfolio
# print(yieldsShock)
# print(discShock)
# print(pricesShock)

init_price <- prices[1] + q1_q2[1]*prices[2] + q1_q2[2]*prices[3]
final_price <- pricesShock[1]+ q1_q2[1]*pricesShock[2] + q1_q2[2]*pricesShock[3]
relative_chg <- 10^4*(final_price-init_price)/init_price
print(round(relative_chg,2))
