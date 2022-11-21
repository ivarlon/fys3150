"""
FYS4150
Project 4
calculating some analytical values for L = 2, T = 1
"""



import numpy as np
from scipy.special import binom

L = 2
N = L**2
def lattice(N_plus):
    s = np.ones(N)
    for i in range(N_plus):
        s[i] = -1
    s = s.reshape(L,L)
    return s

def energy(s):
    E = 0
    for i in range(L):
        for j in range(L):
            E -= s[i,j]* (s[(i+1)%L,j] + s[i, (j+1)%L])
    return E

def magnetisation(s):
    M = 0
    for i in range(L):
        for j in range(L):
            M += s[i,j]
    return M

print("N+, E, M, d")

Z = 0
E = np.zeros(5)
M = np.zeros(5)
d = np.zeros(5)

for N_plus in range(5):
    s = lattice(N_plus)
    E_ = energy(s)
    M_ = magnetisation(s)
    d_ = binom(N, N_plus)
    print(N_plus, E_, M_, d_)
    for array, val in zip([E, M, d], [E_,M_,d_]):
        array[N_plus] = val
    Z += d_*np.exp(-E_)

print("partition function Z =", Z)

mean_E = 1/Z * np.sum(d*np.exp(-E)*E)
print("mean eps", mean_E/N)
mean_E2 = 1/Z * np.sum(d*np.exp(-E)*E**2)
print("mean eps^2", 1/N * mean_E2)
mean_M = 1/Z * np.sum(d*np.exp(-E)*np.abs(M))
print("mean |m|", 1/N * mean_M)
mean_M2 = 1/Z * np.sum(d*np.exp(-E)*M**2)
print("mean m^2", 1/N * mean_M2)
print("CV", 1/N * (mean_E2 - mean_E**2))
print("chi", 1/N * (mean_M2 - mean_M**2))