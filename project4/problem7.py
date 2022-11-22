"""
FYS4150
Project 4

Results from some timing tests of parallelisation
using L = 10, n_cycles = 10^6, T = 1.5-2.5

"""

import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(5,4))

data_paral = {
    5: [10.73, 10.01, 10.46, 10.49, 10.06],
    10: [20.60, 21.07, 20.93, 20.67, 21.24],
    15: [28.38, 29.56, 30.66, 28.26, 28.23]
}

data_no_paral = {
    5: [18.44, 18.46, 18.63, 18.55, 18.52],
    10: [36.95, 36.92, 36.98, 36.93, 36.92],
    15: [55.39, 55.25,55.33, 55.59, 55.39]
}


t_paral = [np.mean(t_vals) for key, t_vals in data_paral.items()]
dt_paral = [np.std(t_vals, ddof=1) for key, t_vals in data_paral.items()]

t_no_paral = [np.mean(t_vals) for key, t_vals in data_no_paral.items()]
dt_no_paral = [np.std(t_vals, ddof=1) for key, t_vals in data_no_paral.items()]

print("Estimated speed-up:")
speed_up = np.array(t_no_paral)/np.array(t_paral)
print(speed_up)
print(speed_up/np.array(t_paral) * np.array(dt_paral)) # estimated uncertainty

n_T = [key for key, vals in data_paral.items()]

plt.errorbar(n_T, t_paral, dt_paral, label="OMP", marker="s")
plt.errorbar(n_T, t_no_paral, dt_no_paral, label="no OMP", marker="s")

plt.xlabel("No. of temperature points"); plt.ylabel("t [s]")
plt.legend()
fig.tight_layout()
plt.savefig("problem7.pdf")
plt.show()