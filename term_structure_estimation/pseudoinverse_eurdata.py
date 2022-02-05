#############################################################################
# Pseudoinverse estimation. Credits to Nicolas Leveroni for the help.
#############################################################################

import pandas as pd
import numpy as np
import datetime as dt
import os

os.chdir("/Users...") # load the file by pointing this directory
import matplotlib.pyplot as plt

# import calendar and obtain a vector "days" that counts the day to each cashflow.
Calendar = pd.read_csv("Calendar.csv")
t_0 = dt.date(2012, 10, 3) # spot date
t_i = [0] * 39
days = [0] * 39
for i in range(0, 39):
    days[i] = dt.datetime.strptime(Calendar.Tenor[i], '%d/%m/%Y').date() - t_0
    days[i] = days[i].days
days = pd.DataFrame({"Days": days})
Calendar = Calendar.join(days)

# import rates
Rates = pd.read_csv("Rates.csv")
# define global variables
N = Calendar.shape[0]
n = Rates.shape[0]
# initializing vectors
P = [0] * N  # discount Factors
F = [0] * N  # Forward rates
L = [0] * N  # Simple Rates
Prices = [0] * n # 17 prices

# Build MATRIX C (Cash flow matrix) in three steps, MM, Futures and Swaps

## Build LIBOR rows for C matrix under C_L
Rates_L = Rates[Rates['Source'] == 'LIBOR']
Rates_L = Rates_L.reset_index(drop=True)
Calendar_L = Calendar[Calendar['TenorCode'].str.startswith('S')]
Calendar_L = Calendar_L.reset_index(drop=True)
n_L = Rates_L.shape[0]
C_L = np.zeros([n_L, N]) ## 4 x 39
for i in range(n_L):
    for t in range(N):
        if Calendar['Tenor'][t] == Rates_L['Maturity Dates'][i]:
            C_L[i, t] = 1 + (Rates_L['Market Quotes'][i] / 100) * (Calendar_L['Days'][i] / 360)

## Build Futures rows for C matrix under C_F
Rates_F = Rates[Rates['Source'] == 'Futures']
Rates_F = Rates_F.reset_index(drop=True)
Calendar_F = Calendar[Calendar['TenorCode'].str.startswith('T')]
Calendar_F = Calendar_F.reset_index(drop=True)
n_F = Rates_F.shape[0]
C_F = np.zeros([n_F, N])
for i in range(n_F):
    start = Calendar_F['Tenor'][i] # start with the first date
    end = Calendar_F['Tenor'][i + 1] # second date, then it goes into a loop
    for t in range(N):
        if Calendar['Tenor'][t] == start:
            C_F[i, t] = -1 # because the FIRST cash flow for the futures is -1 (when you buy it)
        elif Calendar['Tenor'][t] == end:
            C_F[i, t] = 1 + (100 - Rates_F['Market Quotes'][i]) / 100 * (
                        Calendar_F['Days'][i + 1] - Calendar_F['Days'][i]) / 360

## Build Swap rows for C matrix under C_S
Rates_S = Rates[Rates['Source'] == 'Swap']
Rates_S = Rates_S.reset_index(drop=True)
Calendar_S = Calendar[Calendar['TenorCode'].str.startswith('U')]
Calendar_S = Calendar_S.reset_index(drop=True)
n_S = Rates_S.shape[0]
C_S = np.zeros([n_S, N])

for i in range(n_S):
    Maturity = Rates_S['Maturity Dates'][i]
    DaysToMaturity = int(Calendar_S[Calendar_S['Tenor'] == Rates_S['Maturity Dates'][i]]['Days']) # takes the number from the days column coinciding with the dates in both tables
    for t in range(N):
        if Calendar['Days'][t] < DaysToMaturity:
            if Calendar['TenorCode'][t].startswith('U'):
                if Calendar['TenorCode'][t] == 'U1':
                    PrevDate = 0
                else:
                    PrevDate = int(Calendar_S[Calendar_S['Days'] < Calendar['Days'][t]]['Days'].tail(1)) # take the previous to last number of days
                C_S[i, t] = Rates_S['Market Quotes'][i] / 100 * (Calendar['Days'][t] - PrevDate) / 360
        elif Calendar['Days'][t] == DaysToMaturity:
            PrevDate = int(Calendar_S[Calendar_S['Days'] < DaysToMaturity]['Days'].tail(1))
            C_S[i, t] = 1 + Rates_S['Market Quotes'][i] / 100 * (DaysToMaturity - PrevDate + 1) / 360
        else:
            C_S[i, t] = 0


# Put 3 Matrices Together and Print
C = np.vstack((C_L, C_F))
C = np.vstack((C, C_S))
C = pd.DataFrame(C)
C.columns = Calendar['TenorCode']
pd.set_option('display.max_columns', None)
C
# Build Prices Vector - using the instructions in the slides
i = 0
for r in Rates['Source']:
    if r == 'LIBOR':
        Prices[i] = 1
    elif r == 'Futures':
        Prices[i] = 0
    elif r == 'Swap':
        Prices[i] = 1
    i += 1
    

# PSEUDOINVERSE CALCULATION

# Create W Matrix
W = np.zeros([N, N])  # W needs to be NxN and it is a diagonal matrix
W[0, 0] = 1 / np.sqrt(Calendar['Days'][0] / 360) # first entry
for i in range(1, N): # note that this starts from 1
    W[i, i] = 1 / np.sqrt((Calendar['Days'][i] - Calendar['Days'][i - 1]) / 360)

# Create M Matrix
M = np.zeros([N, N])
for i in range(N): # this starts from 0
    M[i, i] = 1
    if i < N - 1:
        M[i + 1, i] = -1
foo = np.zeros(N)  # this is the vector (1,0,0,0,0,0,0...,0)
foo[0] = 1
A = C.dot(np.linalg.inv(M)).dot(np.linalg.inv(W))
term1 = np.transpose(A).dot(np.linalg.inv(A.dot(np.transpose(A))))
term2 = Prices - C.dot(np.linalg.inv(M).dot(np.transpose(foo)))
delta = term1.dot(term2)



# using this to estimate the discount bond prices (so d has to be compared to P)
d = np.linalg.inv(M).dot(np.linalg.inv(W).dot(delta) + np.transpose(foo))

# here we estimate the fwd rates using d instead of P
fwd = [0] * N
fwd[0]=0.095 # this is in percentage points
for f in range(1,N):
    fwd[f] = ((d[f-1]/d[f])-1)/((Calendar['Days'][f] - Calendar['Days'][f-1])/360)*100
plt.plot(Calendar['Days']/365,fwd)
plt.show()

# find the simple forward rate derived using the discounts from the pseudoinverse method.
round(fwd[38],2)
