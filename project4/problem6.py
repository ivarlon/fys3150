import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

L = 20

fig, ax = plt.subplots(nrows=2, figsize=(5,6))
for j,T in enumerate(["1.000000", "2.400000"]):
    # loading data
    filename = f"output_L_{L}_T_{T}_ordered.csv"
    df = pd.read_csv(filename)
    burn_in = 500000
    E_vals = df['E'].values[burn_in:]
    epsilon = E_vals/L**2
    
    print("Mean =", np.mean(epsilon))
    variance = np.var(epsilon,ddof=1)
    print("Var =", variance, "sigma =", variance**0.5)
    print()

    # making histogram
    density, bins = np.histogram(epsilon,bins=300)
    # normalising
    density = density/len(epsilon)
    ax[j].stairs(density, bins, fill=False)
    ax[j].set_title("T = " + T[:3])
    ax[j].set_xlabel("$\\epsilon \ [J]$")
    ax[j].set_ylabel("$p(\\epsilon)$")

fig.tight_layout()
fig.savefig("problem6.pdf")

plt.show()