########################################################################
# Instantaneous convexity adjustment under Vasicek framework.
########################################################################
import numpy as np


k = 0.86
theta = 0.08
sigma = 0.01
r_0 = 0.06
T1 = 1
T0 = 0
t = 0

E_Q_rT = (theta-((sigma**2)/2*(k**2)))*(1-np.exp(-k*(T1-t))) + (sigma**2)/(2*(k**2))*(np.exp(-k*(T1-t))-np.exp(-2*k*(T1-t)))

f_t_T = (sigma**2)/(4*k)*2*((1-np.exp(-k*(T1-t)))/k)*np.exp(-k*(T1-t))-((((k**2)*theta)-(sigma**2)/2)/(k**2))*(np.exp(-k*(T1-t))-1)

a = theta*(1-np.exp(-k*(T1-t)))
convexity_0_1 = a - f_t_T

round(convexity_0_1*10000,2)
