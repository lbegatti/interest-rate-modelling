#################################################################################################################################################
# Credits to Jason Schultz for the help.

# Calibration of the volatility parameters a,b in the Gaussian HJM framework to Cap market data (Cap prices are given)
# Find parameters such that the estimated cap prices are the close as possible to the market data. 
# In practice we minimize the mean squared error of implied volatilities instead of the direct cap prices
# using vega weights: 

# 1) assumption => C_model = C_market + vega*(vol_model - vol_market)
## 2) C_model = C_market + vega*(vol_model - vol_market) --> (vol_model - vol_market)^2 = 1/(vega^2)*(C_model - C_market)^2
### 3) find the value that satifie WLS{Weighted least-squars} optimization: min[(vol_model - vol_market)^2 = 1/(vega^2)*(C_model - C_market)^2] 
#################################################################################################################################################

import numpy as np
from scipy import optimize
from scipy.stats import norm

delta = 0.5
periods = [0, 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4.]
forwardRates = [0.06, 0.08, 0.09, 0.10, 0.10, 0.10, 0.09, 0.09]

# maturities: 1Y, 2Y, 3Y, 4Y;
## semi-annual cashflows 
### first reset date at 6m
capMarketPrices = [0.002, 0.008, 0.012, 0.016]

discountFactors = [1.0]*len(periods)
for i in range(1, len(discountFactors)):
    discountFactors[i] = discountFactors[i-1] / (1. + delta * forwardRates[i-1])


def atmStrikeRate(years):
    # from cap/floor parity we infer that ATM caps strike rates are equal to swap rate with same structure
    return (discountFactors[1] - discountFactors[2 * years]) / \
           (delta * sum([discountFactors[i] for i in range(2, 2 * years + 1)]))


# BlackIvol Cap

def d1_black(period_index, strike, vol):
    return (np.log(forwardRates[period_index] / strike) + 0.5 * (vol ** 2 * periods[period_index])) / (
            vol * np.sqrt(periods[period_index]))


def d2_black(period_index, strike, vol):
    return d1_black(period_index, strike, vol) - vol * np.sqrt(periods[period_index])


def capletBlackPrice(period_index, strike, vol):
    return delta * discountFactors[period_index + 1] * (
            forwardRates[period_index] * norm.cdf(d1_black(period_index, strike, vol), 0, 1)
            - strike * norm.cdf(d2_black(period_index, strike, vol), 0, 1))


def capBlackPrice(year, strike, vol):
    return sum([capletBlackPrice(period_index, strike, vol) for period_index in range(1, 2 * year)])

def impliedBlackVol(year):
    strike = atmStrikeRate(year)
    f = lambda vol: capBlackPrice(year, strike, vol) - capMarketPrices[year-1]
    return optimize.newton(f, 0.5)

impliedVols = [impliedBlackVol(year) for year in range(1, 5)]


# BLACK CAP VEGA --> vega weights
def capletVega(period_index, strike, vol):
    return delta * discountFactors[period_index + 1] * forwardRates[period_index] * \
           np.sqrt(periods[period_index]) * norm.pdf(d1_black(period_index, strike, vol), 0, 1)


def vega(year):
    vol = impliedBlackVol(year)
    strike = atmStrikeRate(year)
    return sum([capletVega(period_index, strike, vol) for period_index in range(1, 2 * year)])

vegas = [vega(year) for year in range(1, 5)]

# Gaussian HJM model --> volatility specification

def HJMvol(beta, v, period_index):
    return (v ** 2 / beta ** 2) * (np.exp(-beta * (periods[period_index])) -
                                   np.exp(- beta * periods[period_index + 1])) ** 2 * \
           (np.exp(2 * beta * (periods[period_index])) - 1) / (2 * beta)

def HJM_d1(period_index, strike, beta, v):
    return (np.log(
        discountFactors[period_index + 1] / discountFactors[period_index] * (1 + delta * strike)) +
            0.5 * HJMvol(beta,v,period_index)) / (np.sqrt(HJMvol(beta, v, period_index)))


def HJM_d2(period_index, strike, beta, v):
    return (np.log(
        discountFactors[period_index + 1] / discountFactors[period_index] * (1 + delta * strike)) -
            0.5 * HJMvol(beta,v,period_index)) / (np.sqrt(HJMvol(beta, v, period_index)))


def HJM_Caplet(period_index, strike, beta, v):
    return (discountFactors[period_index] * norm.cdf(-HJM_d2(period_index, strike, beta, v), 0, 1)) - \
           (1 + delta * strike) * (discountFactors[period_index + 1] * norm.cdf(-HJM_d1(period_index, strike, beta, v), 0, 1))


def HJM_Cap(year, beta, v):
    strike = atmStrikeRate(year)
    return sum([HJM_Caplet(period_index, strike, beta, v) for period_index in range(1, 2 * year)])

def calibrateModel():
    f = lambda p: sum(
        [((1 / vegas[year - 1] ** 2) * (HJM_Cap(year, p[0], p[1]) - capMarketPrices[year - 1]) ** 2) for year in
         range(1, 5)])
    return optimize.minimize(f, [0.01, 0.03])

# get vol parameter b and v, respectively.
calibrateModel()
