import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

amplitudes = ["0.1", "0.4", "0.7"]
n_omegas = 116
omega_array = np.linspace(0.2,2.5,n_omegas)

fig, ax = plt.subplots(figsize=(5,3.5))

for i, f in enumerate(amplitudes):
    n_remaining = pd.read_csv(f'particles_remaining_{f}.csv', header = None).values.ravel()
    ax.plot(omega_array, n_remaining, label = f"f = {f}")
    ax.set_xlabel("$\omega$ [MHz]")
    ax.set_ylabel("% remaining")
ax.legend(loc="lower left")
fig.tight_layout()
plt.show()