##################################################################################################
# In a CIR framework find the mean reversion level.
##################################################################################################

import numpy as np
from scipy.optimize import fsolve

# CIR model
## parameters
a = 0.2
sigma = 0.1
r_0 = 0.05
L = 0.05
T = 1.0
t = 0.0
h = np.sqrt(a**2 + 2*(sigma**2))



# Bond pricing
## using the fact that L(0,1) = 0.05

#ZCB = 1/(1+(T-t)*L)

def A(a,b,h,T,t):
    return (((2*h)*np.exp((a+h)*(T-t)/2))/(2*h + (a + h)*(np.exp((T-t)*h)-1)))**((2*a*b)/(sigma**2))
def B(a,h,T,t):
    return (2*(np.exp((T-t)*h)-1))/((2*h) + (a+h)*(np.exp((T-t)*h)-1))

def P(a,b,h,T,t):
    return A(a,b,h,T,t)*np.exp(-B(a,h,T,t)*r_0)


f = lambda b: (((1/P(a,b,h,T,t))-1)-0.05)
fsolve(f,0.05)
