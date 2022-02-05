##############################################################################################################
# Find the IVol under Black/Bachelier of the ATM swaption
## Find the price when IVol = 50%
##############################################################################################################
import pandas as pd
import numpy as np
from scipy import optimize
from scipy.stats import norm
from math import sqrt,erf

## BLACK VOLATILITY
delta = 0.25
periods = [0.25, 0.5, 0.75, 1., 1.25, 1.5, 1.75, 2.]
forwardRates = [0.06, 0.08, 0.09, 0.10, 0.10, 0.10, 0.09, 0.09]

swaptionprice = [0.01]

# find the simple rate
simpleRates = [forwardRates[0]]
for i in range(1, len(periods)):
    # we derive the simple spot rate from the Law of OnePrice
    simpleRates.append(((1 + periods[i-1] * simpleRates[i - 1]) * (1 + delta * forwardRates[i]) - 1) / periods[i])

# discount factors
discountFactors = []
for i in range(0, len(periods)):
    # this is just calculating the discount factors : 1/(1+delta*L(t))
    discountFactors.append(1. / (1. + periods[i] * simpleRates[i])) # because t=0 so (T-t) becomes T

# strike price ATM
def atmStrikeRate(years):
    return (discountFactors[3] - discountFactors[(4 * years)-1]) / \
           (delta * sum([discountFactors[i] for i in range(4, 4 * years)]))


atmStrikeRate(2)

def d1_black(strike, rswap, vol):
    return (np.log(rswap/strike) + 0.5*((vol**2) * periods[3]))/(vol * np.sqrt(periods[3]))

def d2_black(strike, rswap, vol):
    return d1_black(strike, rswap, vol) - vol * np.sqrt(periods[3])

def payerSwaptionBlackPrice(period_index,rswap, strike,vol):
    return delta * discountFactors[period_index] * (rswap * norm.cdf(d1_black(strike, rswap, vol),0,1) -
                                                    (strike * norm.cdf(d2_black(strike, rswap, vol),0,1)))

def sumswaptionBlack(year, strike, rswap, vol):
    return sum([payerSwaptionBlackPrice(period_index,rswap, strike,vol) for period_index in range(4,4 * year)])

def impliedBlackVol(year):
    strike = atmStrikeRate(year) #ATM_R_swap[year]
    rswap = atmStrikeRate(year)
    f = lambda vol: sumswaptionBlack(year,strike,rswap, vol) - swaptionprice[0]
    return optimize.newton(f,0.05)

impliedVolsBlack = impliedBlackVol(2)
round(impliedVolsBlack*100,2)



# BACHELIER 
def d_bachelier(strike, rswap, vol):
    return (rswap - strike)/(vol*np.sqrt(periods[3]))

def payerSwaptionBachelierPrice(period_index,strike, rswap, vol):
    return delta * (discountFactors[period_index] * vol * np.sqrt(periods[3]) *
                    ((d_bachelier(strike, rswap, vol) * norm.cdf(d_bachelier(strike, rswap, vol),0,1))
                     + norm.pdf(d_bachelier(strike, rswap, vol))))
def sumswaptionBachelier(year, strike, rswap, vol):
    return sum([payerSwaptionBachelierPrice(period_index, strike, rswap, vol) for period_index in range(4,4 * year)])

def impliedBachelierVol(year):
    strike = atmStrikeRate(year) #ATM_R_swap[year]
    rswap = atmStrikeRate(year)
    f = lambda vol: sumswaptionBachelier(year,strike,rswap, vol) - swaptionprice[0]
    #return optimize.brentq(f,-0.5,0.5)
    return optimize.newton(f,0.05)

impliedBach = impliedBachelierVol(2)
round(impliedBach*10000,2)



# Price when IVol = 50%
round(sumswaptionBlack(2,atmStrikeRate(2),0.5)*100,2)
