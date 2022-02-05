################################################################################################
# Pricing cap and floors under Black/Normal volatility
# Implying volatilities given a price
################################################################################################

## BLACK VOLATILITY

delta = 0.25
periods = [0.25, 0.5, 0.75, 1., 1.25, 1.5, 1.75, 2.]
forwardRates = [0.06, 0.08, 0.09, 0.10, 0.10, 0.10, 0.09, 0.09]

# first reset date at 1/4y
capMarketPrices = [0.01]

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
    # from cap/floor parity we infer that ATM caps strike rates are equal to swap rate with same structure
    return (discountFactors[0] - discountFactors[4 * years-1]) / (delta * sum([discountFactors[i] for i in range(1, 4 * years)]))

atmStrikeRate(2) # correct

def d1_black(period_index, strike, vol):
    return (np.log(forwardRates[period_index]/strike) + 0.5*(vol**2 * periods[period_index-1]))/(vol * np.sqrt(periods[period_index - 1]))


def d2_black(period_index, strike, vol):
    return d1_black(period_index, strike, vol) - vol * np.sqrt(periods[period_index-1])

def capletBlackVega(period_index,strike,vol):
    return delta * discountFactors[period_index] * forwardRates[period_index] * np.sqrt(periods[period_index-1]) * norm.pdf(d1_black(period_index, strike, vol),0,1)

def capletBlackPrice(period_index,strike,vol):
    return delta * discountFactors[period_index] * (forwardRates[period_index] * norm.cdf(d1_black(period_index, strike, vol),0,1) - strike * norm.cdf(d2_black(period_index, strike, vol),0,1))

# the capPrice is the sum of the capletPrices
def capBlackPrice(year,strike,vol):
    return sum([capletBlackPrice(period_index, strike, vol) for period_index in range(1,4 * year)])

def impliedBlackVol(year):
    strike = atmStrikeRate(year) #ATM_R_swap[year]
    f = lambda vol: capBlackPrice(year, strike, vol) - capMarketPrices[0]
    return optimize.brentq(f,0.005,0.25)

impliedVols = impliedBlackVol(2)
round(impliedBlackVol(2)*100,2)


# BACHELIER VOLATILITY

def d1_bachelier(period_index, strike, vol_norm):
    return (forwardRates[period_index] - strike)/ (vol_norm * (np.sqrt(periods[period_index-1])))


def bacheliercapletvega(period_index, strike, vol_norm):
    return delta * discountFactors[period_index] * np.sqrt(periods[period_index - 1]) * \
           norm.pdf(d1_bachelier(period_index, strike, vol_norm),0,1)


def bacheliercapletprice(period_index, strike, vol_norm):
    return delta * discountFactors[period_index] * vol_norm * np.sqrt(periods[period_index - 1]) * \
           (d1_bachelier(period_index, strike, vol_norm) * norm.cdf(d1_bachelier(period_index, strike, vol_norm),0,1) +
            norm.pdf(d1_bachelier(period_index, strike, vol_norm),0,1))


def bacheliercapprice(year, strike, vol_norm):
    return sum([bacheliercapletprice(period_index, strike, vol_norm) for period_index in range(1,4 * year)])


def impliedBachVol(year):
    strike = atmStrikeRate(year) #ATM_R_swap[year]
    f = lambda vol: bacheliercapprice(year, strike, vol) - capMarketPrices[0]
    return optimize.brentq(f,0.005,0.25)

round(impliedBachVol(2)*10000,2)


## Practical example to find the price given the volatility.
### suppose the black vol for the cap is 14.1%, what is the price?
atmStrikeRate(2)
round(capBlackPrice(2,atmStrikeRate(2),0.141)*100,2)
