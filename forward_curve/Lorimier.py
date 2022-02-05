######################################################################################
# Estimate the forward curve with Lorimier's method and find a specific yield on ZCB.
######################################################################################

import pandas as pd
import numpy as np

alpha = 0.1
Times = np.array([0,2,3,4,5,7,10,20,30])
Yields = np.array([0,-0.79, -0.73, -0.65, -0.55, -0.33, -0.04, 0.54, 0.73])

# the form should be something like Ax=b

# A = HH
HH = np.array([(Times[i]*Times[j])*(1+min(Times[i],Times[j])/2) - ((min(Times[i],Times[j])**3)/6) for i in range(len(Times)) for j in range(len(Times))]).reshape(9,9)
np.fill_diagonal(HH,HH.diagonal()+(1/alpha))
HH[0,] = [0,2,3,4,5,7,10,20,30]
HH[0:9,0] = [0,2,3,4,5,7,10,20,30]

# b = B
B = np.array([Times[i]*Yields[i] for i in range(len(Times))]).reshape(-1,1)

# x = betas
betas = np.linalg.solve(HH,B)


# find the yield at T = 6
T_6 = 6

iteration_part = [(betas[i]*(Times[i]*min(T_6, Times[i])*(1+((1/2)*min(T_6, Times[i])))-(1/T_6)*min(T_6, Times[i])**3+max((T_6-Times[i]),0)*(Times[i]*(1+Times[i])-(1/2)*(Times[i]**2)))) for i in range(1,len(Times))]

yield_6 = betas[0] + (1/T_6)*sum(iteration_part)
print(np.round(yield_6,2))
