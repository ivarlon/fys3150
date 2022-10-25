import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def get_position(method, n_i):
    # method = "RK4" or "FE"
    # n_i = 4000,8000,16000,32000
    position = pd.read_csv(f'position_{method}_{n_i}.csv', header = None)
    position = position.dropna(axis='columns').values # drop final column of nans
    
    N = position.shape[0]
    r = np.empty(shape=(N, 3))
    v = np.empty(shape=(N, 3))
    for i in range(N):
        for j in range(3):
            r[i, j] = position[i][j]
    return r

def analytic(t):
    B0 = 9.65e1
    V0 = 2.41e6
    d = 500.
    q = 1
    m = 40.08
    omega0 = q*B0/m # cyclotron frequency
    omegaz2 = 2*q*V0/(m*d**2) # axial frequency
    omegaplus = 0.5*( omega0 + np.sqrt(omega0**2 - 2*omegaz2) ) # reduced cyclotron frequency
    omegaminus = 0.5*( omega0 - np.sqrt(omega0**2 - 2*omegaz2) ) # magneton frequency
    r0 = np.array([20, 0, 20])
    v0 = np.array([0, 25, 0])
    z = r0[2]*np.cos(omegaz2**0.5 * t)
    Aplus = (v0[1] + omegaminus*r0[0])/(omegaminus - omegaplus)
    Aminus = -(v0[1] + omegaplus*r0[0])/(omegaminus - omegaplus)
    x = Aplus * np.cos(omegaplus*t) + Aminus * np.cos(omegaminus*t)
    y = - Aplus * np.sin(omegaplus*t) - Aminus*np.sin(omegaminus*t)
    return np.array([x, y, z]).T

n_list = np.array([4000,8000,16000,32000])
T = 50.
fig0, ax0 = plt.subplots(figsize=(5,5))
fig1, ax1 = plt.subplots(figsize=(5,5))
axs = [ax0,ax1]
delta_max = [np.zeros(4), np.zeros(4)] # for calculating error conv. rate
for j, method in enumerate(["RK4", "FE"]):
    ax = axs[j]
    for i, n_i in enumerate(n_list):
        t = np.linspace(0.,T, n_i+1)
        r = get_position(method, n_i)
        r_anal = analytic(t)
        rel_error = np.linalg.norm(r - r_anal, axis = 1) / np.linalg.norm(r_anal, axis = 1)
        ax.plot(t, rel_error, label=f"h=50/{n_i}")
        delta_max[j][i] = np.max( np.abs(r_anal - r) )
        print("delta_max", delta_max[j][i])
    ax.legend()
    ax.set_xlabel(u"Time $t_i$ [\u03bcs]")
    ax.set_ylabel("Relative error")
    if method=="RK4":
        ax.set_title("4th order Runge-Kutta")
    else:
        ax.set_title("Forward Euler")
fig0.tight_layout(); fig1.tight_layout()
h = T/n_list
r_err_RK4 = 1/3 * np.sum( np.log( delta_max[0][1:] / delta_max[0][:-1] ) / np.log( h[1:]/h[:-1] ) )
r_err_FE = 1/3 * np.sum( np.log( delta_max[1][1:] / delta_max[1][:-1] ) / np.log( h[1:]/h[:-1] ) )

print("Error convergence rates:")
print("RK4: ", r_err_RK4)
print("FE: ", r_err_FE)

plt.show()