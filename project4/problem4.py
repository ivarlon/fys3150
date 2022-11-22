"""
FYS4150
Project 4

Ising spin model + Markov Chain Monte Carlo

"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

L_list = [2]
T_list = [1.]
T = T_list[0]

CV_anal = 0.03208
chi_anal = 0.004011

def heat_capacity(E_vals, L, T):
    return 1/L**2 * (1/T**2) * np.var(E_vals)

def susceptibility(M_vals, L, T):
    return 1/L**2 /T * np.var(np.abs(M_vals))

for L in L_list:
    # loading data
    filename = f"output_L_{L}_T_1.000000_disordered.csv"
    df = pd.read_csv(filename)
    # calculating mean energy per spin
    E_vals = df['E'].values
    n_cycles = len(E_vals) - 1
    mean_epsilon = np.mean(E_vals)/L**2
    print("epsilon_mean =", mean_epsilon)
    # calculating mean (abs.) magnetisation per spin
    M_vals = df['M'].values
    mean_m = np.mean(np.abs(M_vals))/L**2
    print("m_mean =", mean_m)
    # calculating heat cap. and susceptibility
    n_points = 20
    C_V = np.zeros(n_points)
    chi = np.zeros(n_points)
    cycles = n_cycles//n_points*np.arange(n_points)
    for i,cycle in enumerate(cycles):
        C_V[i] = heat_capacity(E_vals[:cycle+1], L, T)
        chi[i] = susceptibility(M_vals[:cycle+1], L, T)
    fig1, ax1 = plt.subplots(figsize=(4,3))
    ax1.plot(cycles, C_V, label="numerical")
    ax1.set_xlabel("No. of cycles"); ax1.set_ylabel("$C_V$")
    ax1.plot([cycles[0], cycles[-1]], CV_anal*np.ones(2), 'b--', label="analytic")
    ax1.legend()
    ax1.set_xlim([-n_cycles*0.05,n_cycles*1.05])
    fig1.tight_layout()
    #fig1.savefig("problem4_CV.pdf")

    fig2, ax2 = plt.subplots(figsize=(4,3))
    ax2.plot(cycles, chi, label="numerical")
    ax2.set_xlabel("No. of cycles"); ax2.set_ylabel("$\\chi$")
    ax2.plot([cycles[0], cycles[-1]], chi_anal*np.ones(2), 'b--', label="analytic")
    ax2.legend()
    ax2.set_xlim([-n_cycles*0.05,n_cycles*1.05])
    fig2.tight_layout()
    #fig2.savefig("problem4_chi.pdf")
    plt.show()