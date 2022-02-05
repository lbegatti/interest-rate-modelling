####################################################################################################
# Compute the convexity adjustment under Vasicek framework, which is a Gaussian HJM model
## with deterministic volatility sigma(t,T) {"f" in the code}
####################################################################################################

import numpy as np
from sympy import integrate,Symbol,exp

k = 0.86
theta = 0.08
sigma = 0.01
r_0 = 0.06
T1 = 1.25
T0 = 1
t = 0
gamma = np.sqrt(k**2 + 2*(sigma**2))

B_T1 = (2*(np.exp(gamma*T1)-1))/((gamma+k)*(np.exp(gamma*T1)-1)+2*gamma)
B_T0 = (2*(np.exp(gamma*T0)-1))/((gamma+k)*(np.exp(gamma*T0)-1)+2*gamma)

A_T1 = (-2*k*theta/(sigma**2))*np.log((2*gamma*np.exp((gamma+k)*T1/2))/((gamma+k)*(np.exp(gamma*T1)-1)+2*gamma))
A_T0 = (-2*k*theta/(sigma**2))*np.log((2*gamma*np.exp((gamma+k)*T0/2))/((gamma+k)*(np.exp(gamma*T0)-1)+2*gamma))

P_T1 = np.exp(-A_T1-(B_T1*r_0))     # P(0,1.25)
P_T0 = np.exp(-A_T0-(B_T0*r_0))     # P(0,1)

delta_ratio = (1/(T1-T0))*(P_T0/P_T1)

### we have deterministic forward rate volatility sigma(t,T)
## here we solve the integral with Python.

s = Symbol('s')
u = Symbol('u')
f = sigma*(exp(-k*(u-s)))
first = integrate(f,(u,s,T1))

v = Symbol('v')
f2 = sigma*(exp(-k*(v-s)))
second = integrate(f2,(v,T0,T1))

together = first*second

final_integration = integrate(together,(s,t,T0))

convexity = delta_ratio * (exp(final_integration)-1)

round(convexity*10000,2)
