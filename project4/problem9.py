"""
FYS4150
Project 4

Plotting mean energy, mean magnetisation, heat capacity and susceptibility
for lattice size L = 40, 60, 80, 100
as function of temperature
"""

import numpy as np
from scipy.stats import linregress

TC = 2.269 # analytical value for infinite lattice

L_inv = np.array([1/(20*(2+i)) for i in range(4)])
TC_CV = np.array([2.305,2.293,2.290,2.287])
TC_chi = np.array([2.334,2.321,2.313,2.300])

res_CV = linregress(L_inv, TC_CV)
res_chi = linregress(L_inv, TC_chi)
print("From CV data:")
print(res_CV.intercept, "+/-", res_CV.intercept_stderr)
print("relative error", np.abs((res_CV.intercept-TC)/TC)*100, "%")
print()
print("From chi data:")
print(res_chi.intercept, "+/-", res_chi.intercept_stderr)
print("relative error", np.abs((res_chi.intercept-TC)/TC)*100, "%")

"""# doing linear regression
# setting up design matrix
X = np.stack((np.ones(len(L_inv)), L_inv), axis=1)

coeffs_CV = np.linalg.pinv(X.T@X)@X.T@TC_CV
coeffs_chi = np.linalg.pinv(X.T@X)@X.T@TC_chi
print("CV:", coeffs_CV)
print("chi:", coeffs_chi)"""