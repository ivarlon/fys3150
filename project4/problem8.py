"""
FYS4150
Project 4

Plotting mean energy, mean magnetisation, heat capacity and susceptibility
for lattice size L = 40, 60, 80, 100
as function of temperature
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def heat_capacity(E_vals, L, T):
    return 1/L**2 * (1/T**2) * np.var(E_vals)

def susceptibility(M_vals, L, T):
    return 1/L**2 /T * np.var(np.abs(M_vals))

L_list = [40,60,80,100]

# broad range of temperatures
T_min = 2.1
T_max = 2.4
n_T = 10
T1 = np.linspace(T_min, T_max, n_T)

# narrow range of temperatures
T_min = 2.27
T_max = 2.35
T2 = np.linspace(T_min, T_max, n_T)

# add all points together
T_arr = np.concatenate((T1, T2))
T_arr = np.sort(T_arr)

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8,6))

for L in L_list:
    epsilon = []
    m = []
    CV = []
    chi = []
    for T in T_arr:
        """
        in order to get correct filename, must use 6 decimals for T label
        """
        T_label = (str(np.round(T, 6)) + "0"*10)[:8]
        filename = f"output_L_{L}_T_{T_label}_disordered.csv"
        df = pd.read_csv(filename)
        burn_in = 500000
        E_vals = df['E'].values[burn_in:]
        M_vals = df['M'].values[burn_in:]

        eps_ = np.mean(E_vals)/L**2
        m_ = np.mean(np.abs(M_vals))/L**2
        CV_ = heat_capacity(E_vals, L, T)
        chi_ = susceptibility(M_vals, L, T)
        for arr, val in zip([epsilon, m, CV, chi], [eps_, m_, CV_, chi_]):
            arr.append(val)
    
    # now use data for fine-grained scan
    for idx,arr in enumerate([epsilon, m, CV, chi]):
        axs[idx//2, idx%2].plot(T_arr, arr, label=f"L={L}")
for ax in axs.ravel():
    ax.set_xlabel("T [J/kB]")
    ax.legend()
axs[0,0].set_ylabel("$\\langle \\epsilon \\rangle$")
axs[0,1].set_ylabel("$\\langle |m| \\rangle$")
axs[1,0].set_ylabel("$C_V$")
axs[1,1].set_ylabel("$\\chi$")
fig.tight_layout()
#fig.savefig("problem8.pdf")
plt.show()
