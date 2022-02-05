######################################################################################################################################################
# Pseudoinverse application to UK gilt data to find discount curve and then price of a portfolio of cashflows
## sensitivity to the first bond of the portfolio --> units of BOND_1 to buy to hedge the portfolio against fluctuation in price of BOND_1
######################################################################################################################################################

import pandas as pd
import numpy as np
from datetime import datetime
import datetime as dt
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt

t0 = dt.date(1996, 9, 4)  # spot date
spot_date = '04/09/1996' # as string

cashflow_dates = ['26/09/1996', '13/10/1996', '06/11/1996', '15/11/1996', '07/12/1996', '19/01/1997', '27/02/1997',
                  '03/03/1997', '08/03/1997', '26/03/1997', '13/04/1997', '06/05/1997', '07/06/1997', '19/07/1997',
                  '27/08/1997', '03/09/1997', '08/09/1997', '26/09/1997', '13/10/1997', '06/11/1997', '07/12/1997',
                  '19/01/1998', '27/02/1998', '03/03/1998', '08/03/1998', '26/03/1998', '13/04/1998', '06/05/1998',
                  '07/06/1998', '27/08/1998', '03/09/1998', '08/09/1998', '26/09/1998', '13/10/1998', '06/11/1998',
                  '07/12/1998', '27/02/1999', '03/03/1999', '08/03/1999', '26/03/1999', '13/04/1999', '06/05/1999',
                  '07/06/1999', '27/08/1999', '03/09/1999', '08/09/1999', '13/10/1999', '06/11/1999', '07/12/1999',
                  '27/02/2000', '03/03/2000', '08/03/2000', '13/04/2000', '06/05/2000', '07/06/2000', '27/08/2000',
                  '08/09/2000', '13/10/2000', '06/11/2000', '07/12/2000', '27/02/2001', '08/03/2001', '13/04/2001',
                  '06/05/2001', '07/06/2001', '27/08/2001', '08/09/2001', '13/10/2001', '06/11/2001', '07/12/2001',
                  '27/02/2002', '08/03/2002', '13/04/2002', '07/06/2002', '27/08/2002', '08/09/2002', '13/10/2002',
                  '07/12/2002', '08/03/2003', '13/04/2003', '07/06/2003', '08/09/2003', '13/10/2003', '07/12/2003',
                  '08/03/2004', '13/04/2004', '07/06/2004', '08/09/2004', '13/10/2004', '07/12/2004', '08/03/2005',
                  '13/04/2005', '07/06/2005', '08/09/2005', '13/10/2005', '07/12/2005', '08/03/2006', '13/04/2006',
                  '08/09/2006', '13/10/2006', '13/04/2007', '13/10/2007', '13/04/2008', '13/10/2008']

info = pd.DataFrame(
    {
        'BOND': ['BOND1', 'BOND2', 'BOND3', 'BOND4', 'BOND5', 'BOND6', 'BOND7', 'BOND8', 'BOND9'],
        'ANNUAL_COUPON_PCT': [10, 9.75, 12.25, 9, 7, 9.75, 8.5, 7.75, 9],
        'NEXT_COUPON': ['15/11/1996', '19/01/1997', '26/09/1996', '03/03/1997', '06/11/1996', '27/02/1997',
                        '07/12/1996', '08/03/1997', '13/10/1996'],
        'MATURITY_DATE': ['15/11/1996', '19/01/1998', '26/03/1999', '03/03/2000', '06/11/2001', '27/08/2002',
                          '07/12/2005', '08/09/2006', '13/10/2008'],
        'PRICE': [103.82, 106.04, 118.44, 106.28, 101.15, 111.06, 106.24, 98.49, 110.87]

    }
)

# pseudo-inverse
N = len(cashflow_dates)
W = np.zeros([N, N])  # diagonal matrix

W[0,0] = 1/np.sqrt((min(dt.datetime.strptime(cashflow_dates[0],'%d/%m/%Y').day,30)+max(30 - dt.datetime.strptime(spot_date,'%d/%m/%Y').day,0))/360 +
                   (dt.datetime.strptime(cashflow_dates[0],'%d/%m/%Y').month - dt.datetime.strptime(spot_date,'%d/%m/%Y').month - 1)/12 +
                   (dt.datetime.strptime(cashflow_dates[0],'%d/%m/%Y').year - dt.datetime.strptime(spot_date,'%d/%m/%Y').year))
for i in range(1, N):
    W[i,i] = 1/np.sqrt((min(dt.datetime.strptime(cashflow_dates[i],'%d/%m/%Y').day,30)+max(30 - dt.datetime.strptime(cashflow_dates[i-1],'%d/%m/%Y').day,0))/360 +
                       (dt.datetime.strptime(cashflow_dates[i],'%d/%m/%Y').month - dt.datetime.strptime(cashflow_dates[i-1],'%d/%m/%Y').month - 1)/12 +
                       (dt.datetime.strptime(cashflow_dates[i],'%d/%m/%Y').year - dt.datetime.strptime(cashflow_dates[i-1],'%d/%m/%Y').year))

# the M matrix
M = np.zeros([N, N])
for i in range(N):
    M[i, i] = 1
    if i < N - 1:
        M[i + 1, i] = -1
foo = np.zeros(N)  # vector
foo[0] = 1

# cashflow matrix
cashflow_dates = pd.DataFrame({"date": cashflow_dates})
dates = [datetime.strptime(date, '%d/%m/%Y') for date in cashflow_dates.date]
info['NEXT_COUPON'] = [datetime.strptime(x, '%d/%m/%Y') for x in info.NEXT_COUPON]
info['MATURITY_DATE'] = [datetime.strptime(x, '%d/%m/%Y') for x in info.MATURITY_DATE]


sparse_matrix = np.zeros((len(info.BOND), len(dates)))

for i in range(0, len(info.BOND)):  # rows
    for ii in range(0, len(dates)):  # column
        if (info.NEXT_COUPON[i] != dates[ii]) & (info.NEXT_COUPON[i] != info.MATURITY_DATE[i]):
            new_date = info.NEXT_COUPON[i] + relativedelta(months=+6)
            data = dates[ii]
            for g in range(0, len(dates)):
                if (new_date != data) & (new_date != info.MATURITY_DATE[i]):
                    new_date += relativedelta(months=+6)
                elif (new_date == data) & (new_date == info.MATURITY_DATE[i]):
                    sparse_matrix[i,ii] = (info.ANNUAL_COUPON_PCT[i]/2) + 100
                elif (new_date == data) & (new_date != info.MATURITY_DATE[i]):
                    sparse_matrix[i, ii] = (info.ANNUAL_COUPON_PCT[i]/2)
                else:
                    pass
        elif (info.NEXT_COUPON[i] == dates[ii]) & (info.NEXT_COUPON[i] == info.MATURITY_DATE[i]):
            sparse_matrix[i, ii] = (info.ANNUAL_COUPON_PCT[i]/2) + 100
        elif (info.NEXT_COUPON[i] == dates[ii]) & (info.NEXT_COUPON[i] != info.MATURITY_DATE[i]):
            sparse_matrix[i,ii] = (info.ANNUAL_COUPON_PCT[i]/2)
        else:
            pass

Prices = list(info.PRICE)


A = sparse_matrix.dot(np.linalg.inv(M)).dot(np.linalg.inv(W))
term1 = np.transpose(A).dot(np.linalg.inv(A.dot(np.transpose(A))))
term2 = Prices - sparse_matrix.dot(np.linalg.inv(M).dot(np.transpose(foo)))
delta = term1.dot(term2)

d = np.linalg.inv(M).dot(np.linalg.inv(W).dot(delta) + np.transpose(foo)) # discount factor curve

PV = 80*d[74] + 100*d[95] + 60*d[98] + 250*d[103]
portfolio_cashflow = [80.0,100.0,60.0,250.0]
M_inv = np.linalg.inv(M)
Mshape=M_inv.shape
W_inv = np.linalg.inv(W)
Wshape=W_inv.shape
term1shape=term1.shape

b = (np.linalg.inv(M).dot(np.linalg.inv(W).dot(term1)))

sens_to_bond1=b[74,0]*portfolio_cashflow[0]+b[95,0]*portfolio_cashflow[1]+b[98,0]*portfolio_cashflow[2]+b[103,0]*portfolio_cashflow[3]
round(sens_to_bond1,2)
