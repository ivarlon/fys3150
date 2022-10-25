import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

n_omegas = 151
omega_array = np.linspace(1.3,1.45,n_omegas)

fig, ax = plt.subplots(figsize=(5,4))

for i, partInt in enumerate(["off", "on"]):
    n_remaining = pd.read_csv(f'particles_remaining_finestruct_partint_{partInt}.csv', header = None).values.ravel()
    ax.plot(omega_array, n_remaining, label = f"Int. {partInt}")
    ax.set_xlabel("$\omega$ [MHz]")
    ax.set_ylabel("% remaining")
ax.legend(loc="best")
fig.tight_layout()
plt.show()