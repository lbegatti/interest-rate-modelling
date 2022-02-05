####################################################################################################
# Pseudoinverse method to find the discount curve with swap rates
## price of an ATM cap at specific maturity under HJM framework
### Black IVol of the cap
#### Normal IVol 
####################################################################################################

import pandas as pd
import numpy as np
from datetime import datetime
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt

spot_date = '24/01/2022'

info = pd.DataFrame(
    {
        'MATURITY': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30],
        'SWAP_RATE': [0.36, 0.52, 0.93, 1.21, 1.46, 1.66, 1.84, 1.99, 2.13, 2.21, 2.63, 2.73, 2.71],
        'SWAP_PRICES': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    }
)
date = datetime.strptime(spot_date, '%d/%m/%Y')
start_payments = date + relativedelta(days=+2)


def date_count(start, end, day_of_month=1):
    dates = [start]
    next_date = start.replace(day=day_of_month)
    if day_of_month > start.day:
        dates.append(next_date)
    while next_date < end.replace(day=day_of_month):
        next_date += relativedelta(next_date, months=+6)
        dates.append(next_date)
    return dates


col_names = ['SWAP1', 'SWAP2', 'SWAP3', 'SWAP4', 'SWAP5', 'SWAP6',
             'SWAP7', 'SWAP8', 'SWAP9', 'SWAP10', 'SWAP11', 'SWAP12', 'SWAP13']

flows_dates = pd.DataFrame(columns=col_names, index=range(1, 61))
maturity_dates = pd.DataFrame(columns=col_names, index=[0])
for i in range(0, len(info.MATURITY)):
    cashflow_swap = date_count(start_payments, start_payments + relativedelta(years=+info.MATURITY[i]), day_of_month=26)
    date_conversion = [datetime.strftime(x, '%d/%m/%Y') for x in cashflow_swap]
    maturity_dates.iloc[:, i] = date_conversion[-1]
    flows_dates.iloc[:, i] = pd.Series(date_conversion)

# PSEUDOINVERSE
N = len(flows_dates.SWAP13)
W = np.zeros([N, N])  # diagonal matrix

W[0, 0] = 1 / np.sqrt(0.5)
for i in range(1, N):
    W[i, i] = 1 / np.sqrt(0.5)

# the M matrix
M = np.zeros([N, N])
for i in range(N):
    M[i, i] = 1
    if i < N - 1:
        M[i + 1, i] = -1
foo = np.zeros(N)  # vector
foo[0] = 1

# cashflow matrix
cashflow_dates = pd.DataFrame({"date": flows_dates.SWAP13})
dates = [datetime.strptime(date, '%d/%m/%Y') for date in cashflow_dates.date]
info['NEXT_SWAP_PAYMENT'] = start_payments
info['MATURITY_DATE'] = [datetime.strptime(x, '%d/%m/%Y') for x in maturity_dates.iloc[0]]

sparse_matrix = np.zeros((len(info.SWAP_RATE), len(dates)))

for i in range(0, len(info.SWAP_RATE)):  # rows
    for ii in range(0, len(dates)):  # column
        if (info.NEXT_SWAP_PAYMENT[i] != dates[ii]) & (info.NEXT_SWAP_PAYMENT[i] != info.MATURITY_DATE[i]):
            new_date = info.NEXT_SWAP_PAYMENT[i] + relativedelta(months=+6)
            data = dates[ii]
            for g in range(0, len(dates)):
                if (new_date != data) & (new_date != info.MATURITY_DATE[i]):
                    new_date += relativedelta(months=+6)
                elif (new_date == data) & (new_date == info.MATURITY_DATE[i]):
                    sparse_matrix[i, ii] = 1 + info.SWAP_RATE[i] / 100 * 0.5
                elif (new_date == data) & (new_date != info.MATURITY_DATE[i]):
                    sparse_matrix[i, ii] = info.SWAP_RATE[i] / 100 * 0.5
                else:
                    pass
        elif (info.NEXT_SWAP_PAYMENT[i] == dates[ii]) & (info.NEXT_SWAP_PAYMENT[i] == info.MATURITY_DATE[i]):
            sparse_matrix[i, ii] = 1 + info.SWAP_RATE[i] / 100 * 0.5
        elif (info.NEXT_SWAP_PAYMENT[i] == dates[ii]) & (info.NEXT_SWAP_PAYMENT[i] != info.MATURITY_DATE[i]):
            sparse_matrix[i, ii] = info.SWAP_RATE[i] / 100 * 0.5
        else:
            pass

Prices = list(info.SWAP_PRICES)
##
A = sparse_matrix.dot(np.linalg.inv(M)).dot(np.linalg.inv(W))
term1 = np.transpose(A).dot(np.linalg.inv(A.dot(np.transpose(A))))
term2 = Prices - sparse_matrix.dot(np.linalg.inv(M).dot(np.transpose(foo)))
delta = term1.dot(term2)

d = np.linalg.inv(M).dot(np.linalg.inv(W).dot(delta) + np.transpose(foo))  # discount factor curve
plt.plot(d)



## Pricing Gaussian HJM model

from scipy.stats import norm
from scipy import optimize

periods = [i for i in np.arange(0.5, 30.5, 0.5)]
v1 = 0.01
v2 = 0.02
beta1 = 0.3
beta2 = 0.5

discountFactors = list(d)
time_increment = 0.5

sr = []  # simple rates
for i in range(0, len(periods)):
    sr.append(((1 / discountFactors[i]) - 1) / time_increment)

fwdrates = [sr[0]]
for i in range(1, len(periods)):
    fwdrates.append((1 / time_increment) * ((discountFactors[i - 1] / discountFactors[i]) - 1))


def atmStrikeRate(years):
    # from cap/floor parity we infer that ATM caps strike rates are equal to swap rate with same structure
    return (discountFactors[0] - discountFactors[2 * years - 1]) / (time_increment * sum([discountFactors[i] for i in
                                                                                          range(1,
                                                                                                2 * years)]))  # this is to skip the unnecessary swap rates, i.e. odd index


def HJMvol(v1, v2, beta1, beta2, maturity):  # maturity will be periods[period_index]
    return (((v1 ** 2 / beta1 ** 2) * (np.exp(-beta1 * (maturity - 0.5)) - np.exp(- beta1 * maturity)) ** 2 * (
            np.exp(2 * beta1 * (maturity - 0.5)) - 1) / (2 * beta1)) +
            ((v2 ** 2 / beta2 ** 2) * (np.exp(-beta2 * (maturity - 0.5)) - np.exp(- beta2 * maturity)) ** 2 * (
                    np.exp(2 * beta2 * (maturity - 0.5)) - 1) / (2 * beta2)))


def HJM_d1(period_index, strike, v1, v2, beta1, beta2):
    return (np.log(discountFactors[period_index] / discountFactors[period_index - 1] * (1 + time_increment * strike)) +
            0.5 * HJMvol(v1, v2, beta1, beta2, periods[period_index])) / (
               np.sqrt(HJMvol(v1, v2, beta1, beta2, periods[period_index])))


def HJM_d2(period_index, strike, v1, v2, beta1, beta2):
    return (np.log(discountFactors[period_index] / discountFactors[period_index - 1] * (1 + time_increment * strike)) -
            0.5 * HJMvol(v1, v2, beta1, beta2, periods[period_index])) / (
               np.sqrt(HJMvol(v1, v2, beta1, beta2, periods[period_index])))


def HJM_Caplet(period_index, strike, v1, v2, beta1, beta2):
    return (discountFactors[period_index - 1] * norm.cdf(-HJM_d2(period_index, strike, v1, v2, beta1, beta2), 0, 1)) - \
           (1 + time_increment * strike) * (
                   discountFactors[period_index] * norm.cdf(-HJM_d1(period_index, strike, v1, v2, beta1, beta2), 0,
                                                            1))


def HJM_Cap(year, v1, v2, beta1, beta2):
    price = []
    strike = atmStrikeRate(year)
    impl_vol = sum([HJM_Caplet(period_index, strike, v1, v2, beta1, beta2) for period_index in range(1, 2 * year)])
    price.append(impl_vol)
    return price

cap_price = HJM_Cap(30, v1=v1, v2=v2, beta1=beta1, beta2=beta2)


### Black IVol in pct

def d1_black(period_index, strike, vol):
    return (np.log(fwdrates[period_index] / strike) + 0.5 * (vol ** 2 * periods[period_index - 1])) / (
            vol * np.sqrt(periods[period_index - 1]))


def d2_black(period_index, strike, vol):
    return d1_black(period_index, strike, vol) - vol * np.sqrt(periods[period_index - 1])


def capletBlackPrice(period_index, strike, vol):
    return time_increment * discountFactors[period_index] * (
            fwdrates[period_index] * norm.cdf(d1_black(period_index, strike, vol), 0, 1)
            - strike * norm.cdf(d2_black(period_index, strike, vol), 0, 1))


def capBlackPrice(year, strike, vol):
    return sum([capletBlackPrice(period_index, strike, vol) for period_index in range(1, 2 * year)])


def impliedBlackVol(year):
    strike = atmStrikeRate(year)
    f = lambda vol: capBlackPrice(year, strike, vol) - cap_price[0]
    return optimize.newton(f, 0.5, )


impliedBlackVol(30)
impliedVolsBlack = impliedBlackVol(30)
round(impliedVolsBlack * 100, 2)


#### BACHELIER IVol in bps

def d1_bachelier(period_index, strike, vol_norm):
    return (fwdrates[period_index] - strike) / (vol_norm * (np.sqrt(periods[period_index - 1])))


def bacheliercapletprice(period_index, strike, vol_norm):
    return time_increment * discountFactors[period_index] * vol_norm * np.sqrt(periods[period_index - 1]) * \
           (d1_bachelier(period_index, strike, vol_norm) * norm.cdf(d1_bachelier(period_index, strike, vol_norm), 0,
                                                                    1) +
            norm.pdf(d1_bachelier(period_index, strike, vol_norm), 0, 1))


def bacheliercapprice(year, strike, vol_norm):
    return sum([bacheliercapletprice(period_index, strike, vol_norm) for period_index in range(1, 2 * year)])


def impliedBachVol(year):
    strike = atmStrikeRate(year)  # ATM_R_swap[year]
    f = lambda vol: bacheliercapprice(year, strike, vol) - cap_price[0]
    return optimize.newton(f, 0.005)


round(impliedBachVol(30) * 10000, 2)
