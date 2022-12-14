"""
FYS4150

Project 5

Solving Schr. eq.

Plotting p(y|x=0.8, t=0.002)
"""

import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa

T = 0.002

for potential in ["single_slit", "double_slit", "triple_slit"]:
    fig, ax = plt.subplots(figsize=(5,4))
    p = pa.cube()
    p.load(potential + "_T_" + (str(T))[:5] + "_p.bin")
    p = np.array(p)[-1] # choose probability for final timestep
    M = p.shape[0] + 2 # no. of points in grid
    x_idx = int((M-1)*0.8) - 1 # idx. of x = 0.8
    p_y = p[x_idx,:]
    p_y = p_y/np.sum(p_y)
    y = np.linspace(0.,1.,M-2)
    if potential=="single_slit":
        label = "enkeltspalte"
    elif potential=="double_slit":
        label = "dobbelspalte"
    elif potential=="triple_slit":
        label = "trippelspalte"
    ax.plot(y, p_y, label=label)
    ax.set_xlabel("$y$")
    ax.set_ylabel("$p(y|x=0.8, t=0.002)$")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"problem9_{potential}.pdf")
    plt.show()
