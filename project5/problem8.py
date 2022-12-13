import numpy as np
import pyarma as pa
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl

T = 0.002 # integration time

potential = "double_slit"

# load potential matrix in order to visualise potential
V = pd.read_csv(potential + ".dat", sep=" ", header=None)

# load probability function
p = pa.cube()
p.load(potential + "_T_" + (str(T))[:5] + "_p.bin")
p = np.array(p)

n_timepoints = p.shape[0]
print(n_timepoints)

# load real part of wave function
u_real = pa.cube()
u_real.load(potential + "_T_" + (str(T))[:5] + "_real.bin")
u_real = np.array(u_real)

# load imaginary part of wave function
u_imag = pa.cube()
u_imag.load(potential + "_T_" + (str(T))[:5] + "_imag.bin")
u_imag = np.array(u_imag)

# create plots
figsize=(8,6)
fig_p, axs_p = plt.subplots(ncols=3, figsize=figsize)
fig_re, axs_re = plt.subplots(ncols=3, figsize=figsize)
fig_im, axs_im = plt.subplots(ncols=3, figsize=figsize)

for i in range(3):
    ax_p = axs_p[i]
    ax_re = axs_re[i]
    ax_im = axs_im[i] 
    if i==2:
        idx = n_timepoints -1 # final timestep
    else:
        idx = n_timepoints*i//2
    
    p_i = p[idx]
    u_real_i = u_real[idx]
    u_imag_i = u_imag[idx]
    
    for ax, arr in zip([ax_p, ax_re, ax_im], [p_i, u_real_i, u_imag_i]):
        arr[np.where(V!=0)] = None # set region of barrier to None in order to highlight it
        arr = np.rot90(arr) # rotate to get x horizontal and y vertical
        if ax == ax_p:
            arr = arr
            img = ax.imshow(arr)#,norm=norm)
        else:
            norm = mpl.cm.colors.Normalize(vmin=0.0, vmax=np.max(arr))
            img = ax.imshow(arr, norm=norm)
        plt.get_cmap().set_bad(color="black") # highlights barrier
        cbar = plt.colorbar(img, ax=ax, shrink=0.3)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xlabel("$x$")
        ax.set_ylabel("$y$")
        ax.set_title(f"$t = {T*i/2}$")
        if i==2:
            if ax==ax_p:
                label="$p$"
            elif ax==ax_re:
                label="Re($u$)"
            elif ax==ax_im:
                label="Im($u$)"
            cbar.set_label(label)

for fig in [fig_p, fig_re, fig_im]:
    fig.tight_layout()
# #fig.savefig(potential+".pdf")
plt.show()
