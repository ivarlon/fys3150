"""
FYS4150
Project 4
calculating some analytical values for L = 2, T = 1
"""



import numpy as np
from scipy.special import binom

L = 2
N = L**2
def lattice(N_plus,permute=False):
    s = np.ones(N)
    if not permute:
        for i in range(N_plus):
            s[i] = -1
    else:
        for i in range(N_plus):
            s[1+i] = -1
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

E = []
M = []
d = []

for N_plus in range(5):
    s = lattice(N_plus)
    E_ = energy(s)
    M_ = magnetisation(s)
    d_ = binom(N, N_plus)
    print(N_plus, E_, M_, d_ - (N_plus==2)*2)
    for array, val in zip([E, M, d], [E_,M_,d_]):
        array.append(val)
    if N_plus==2:
        s = lattice(N_plus, permute=True)
        E_ = energy(s)
        M_ = magnetisation(s)
        d_ = 2
        d[-1] -= d_
        print(N_plus, E_, M_, d_)
        for array, val in zip([E, M, d], [E_,M_,d_]):
            array.append(val)

E = np.array(E)
M = np.array(M)
d = np.array(d)

Z = np.sum(d*np.exp(-E))
print("partition function Z =", Z)
#Z = 2*np.exp(8) + 12 + 2*np.exp(-8)
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