import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

def get_position_and_velocity(partint):
    # partint = "on" or "off"
    position = pd.read_csv(f'position_partints_{partint}.csv', header = None)
    position = position.dropna(axis='columns').values # drop final column of nans
    velocity = pd.read_csv(f'velocity_partints_{partint}.csv', header = None)
    velocity = velocity.dropna(axis='columns').values # drop final column of nans

    N = position.shape[0] # no. of points = nSteps + 1
    nParticles = position.shape[1]//3 # each particle has three r components
    r = np.empty(shape=(N, nParticles, 3))
    v = np.empty(shape=(N, nParticles, 3))
    for i in range(N):
        for j in range(nParticles):
            for k in range(3):
                r[i, j, k] = position[i][j*3 + k]
                v[i, j, k] = velocity[i][j*3 + k]
    return r, v

nParticles = 2
fig_x, axs_x = plt.subplots(figsize=(8,4), ncols=2)
fig_z, axs_z = plt.subplots(figsize=(8,4),ncols=2)
fig_off = plt.figure(figsize=(4,4))
ax_off = fig_off.add_subplot(1,1,1,projection="3d")
fig_on = plt.figure(figsize=(4,4))
ax_on = fig_on.add_subplot(1,1,1,projection="3d")
ax_3d = [ax_off, ax_on]
for i, partint in enumerate(["off", "on"]):
    
    r, v = get_position_and_velocity(partint)
    N = r.shape[0]
    T = 50.
    t = np.linspace(0., T, N)

    for j in range(nParticles):
        #plotting phase spaces
        x = r[:,j,0]
        vx = v[:,j,0]
        y = r[:,j,1]
        z = r[:,j,2]
        vz = v[:,j,2]
        axs_x[i].plot(x, vx, label=f"particle {j+1}")
        axs_x[i].scatter([x[0], x[-1]], [vx[0], vx[-1]])
        axs_z[i].plot(z, vz, label=f"particle {j+1}")
        axs_z[i].scatter([z[0], z[-1]], [vz[0], vz[-1]])    
        ax_3d[i].plot(x,y,z, label=f"particle {j+1}")
    
    axs_x[i].set_xlabel(u"$x$ [\u03bcm]")
    axs_x[i].set_ylabel(u"$v_x$ [\u03bcm/\u03bcs]")
    axs_x[i].set_title(f"Particle interactions {partint}")
    axs_z[i].set_xlabel(u"$z$ [\u03bcm]")
    axs_z[i].set_ylabel(u"$v_z$ [\u03bcm/\u03bcs]")
    ax_3d[i].set_xlabel(u"$x$ [\u03bcm]")
    ax_3d[i].set_ylabel(u"$y$ [\u03bcm]")
    ax_3d[i].set_zlabel(u"$z$ [\u03bcm]")
    axs_z[i].set_title(f"Particle interactions {partint}")
    ax_3d[i].set_title(f"Particle interactions {partint}")
    axs_x[i].legend(); axs_z[i].legend(); ax_3d[i].legend()
    axs_x[i].axis("equal"); axs_z[i].axis("equal")
for fig in [fig_x, fig_z, fig_off, fig_on]:
    fig.tight_layout()
plt.show()