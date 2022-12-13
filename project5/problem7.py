import numpy as np
import pyarma as pa
import matplotlib.pyplot as plt

T = 0.008
for potential in ["no_potential"]:#, "double_slit"]:
    U = pa.cube()
    U.load(potential + "_p.bin")
    U = np.array(U)
    p = np.sum(U, axis=1)
    p = np.sum(p, axis=1)
    t = np.linspace(0.,T, len(p))
    fig, ax = plt.subplots(figsize=(4,4))
    ax.plot(t, p-1)
    if potential=="no_potential":
        title = "No barrier"
    elif potential=="double_slit":
        title = "Double slit"
    ax.set_title(title)
    ax.set_xlabel("t")
    ax.set_ylabel("$p-1$")
    fig.tight_layout()
    plt.show()
    #fig.savefig(potential+".pdf")
