import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

L = 20

fig1, ax1 = plt.subplots(ncols=2, figsize=(9,4))
fig2, ax2 = plt.subplots(ncols=2, figsize=(9,4))
for j,T in enumerate(["1.000000", "2.400000"]):
    for ordered in ["ordered", "disordered"]:
        # loading data
        filename = f"output_L_{L}_T_{T}_{ordered}.csv"
        df = pd.read_csv(filename)
        E_vals = df['E'].values#[300000:]
        M_vals = df['M'].values#[300000:]
        
        # creating arrays to plot
        n_points = 20
        epsilon = np.zeros(n_points)
        m = np.zeros(n_points)
        n_cycles = len(E_vals) - 1
        cycles = n_cycles//n_points*np.arange(1,n_points+1)
        for i,cycle in enumerate(cycles):
            epsilon[i] = np.mean(E_vals[:cycle])/L**2
            m[i] = np.mean(np.abs(M_vals[:cycle]))/L**2
        
        # plotting
        ax1[j].plot(cycles, epsilon, label=ordered)
        ax1[j].set_title("T = " + T[:3])
        ax2[j].plot(cycles, m, label=ordered)
        ax2[j].set_title("T = " + T[:3])

for ax in ax1:
    ax.set_xlabel("No. of cycles"); ax.set_ylabel("$\\langle \\epsilon \\rangle$")
    ax.legend()
fig1.tight_layout()
#fig1.savefig("problem5_eps.pdf")
for ax in ax2:
    ax.set_xlabel("No. of cycles"); ax.set_ylabel("$\\langle |m| \\rangle$")
    ax.legend()
fig2.tight_layout()
#fig2.savefig("problem5_m.pdf")
plt.show()