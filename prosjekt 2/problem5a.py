"""
FYS4150
project 2
problem 5a)

plotting matrix size N vs. number of iterations needed for Jacobi to converge

"""

import numpy as np
import matplotlib.pyplot as plt

N = np.array([4, 5, 6, 7, 8, 9, 10, 15, 20,40,60,80,100]) # matrix size
#i = np.array([20,35,52,70,101,111,154,313,595]) # iterations needed to converge
i = np.array([6,30,35,65,80,121,129,342,621,2640,5911,10685,16316]) # iterations needed to converge

plt.plot(N,i, label="$i(N)$") # plotting datapoints
plt.plot(N, 0.5*np.pi*N**2, linestyle="--", label="$\pi/2 \ N^2$") # approximate behaviour
plt.xlabel("N")
plt.ylabel("Iterations")
plt.legend()
plt.suptitle("Convergence rate of Jacobi")
plt.savefig("problem5a.pdf")
plt.show()