"""
FYS4150
project 1
problem 5a)

plotting matrix size N vs. number of iterations needed for Jacobi to converge

"""

import numpy as np
import matplotlib.pyplot as plt

N = np.array([4, 5, 6, 7, 8, 9, 10, 15, 20]) # matrix size
i = np.array([20,35,52,70,101,111,154,313,595]) # iterations needed to converge

plt.plot(N,i, label="i(N)") # plotting datapoints
plt.plot(N, 1.5*N**2, linestyle="--", label="3/2 N^2") # approximate behaviour
plt.xlabel("N")
plt.ylabel("Iterations")
plt.legend()
plt.suptitle("Convergence rate of Jacobi")
plt.show()