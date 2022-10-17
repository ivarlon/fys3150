"""
FYS4150
project 2
problem 6

plotting eigenvectors

"""

import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

# boundary values
x0 = 0.; xn = 1.
u0 = 0.; un = 0.

# load eigenvector data from Armadillo
eigenvecs = pa.mat()
eigenvecs.load("eigenvecs.bin")

N = len(eigenvecs)**0.5
N = int(N)
nsteps = N+1

x = np.linspace(x0, xn, nsteps+1)

eigenvecs = np.asarray(eigenvecs)
eigenvecs = eigenvecs.reshape(N,N)

# plot
for i in range(3):
    u_i = np.zeros(nsteps+1)
    u_i[0] = u0
    u_i[-1] = un
    u_i[1:-1] = eigenvecs[i]
    plt.plot(x, u_i, label=f"u_{i+1}")


plt.legend()
plt.xlabel("$\hat{x}$"); plt.ylabel("$u$")
plt.savefig(f"problem6 N = {N}.pdf")
plt.show()