"""
Plots trajectory of particle in Penning trap
Prints expected frequency of oscillations in z
"""


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

data = pd.read_csv('position.csv', header = None)
data = data.dropna(axis='columns')
data = data.values

N = data.shape[0]
nParticles = data.shape[1]//3
r = np.empty(shape=(N, nParticles, 3))
for i in range(N):
    for j in range(nParticles):
        for k in range(3):
            r[i, j, k] = data[i][j*3 + k]

def analytic(t):
    B0 = 9.65e1
    V0 = 2.41e6
    d = 500.
    q = 1
    m = 40.08
    omega0 = q*B0/m # cyclotron frequency
    omegaz2 = 2*q*V0/(m*d**2) # axial frequency
    omegaplus = 0.5*( omega0 + np.sqrt(omega0**2 - 2*omegaz2) ) # reduced cyclotron frequency
    omegaminus = 0.5*( omega0 - np.sqrt(omega0**2 - 2*omegaz2) ) # magneton frequency
    r0 = np.array([20, 0, 20])
    v0 = np.array([0, 25, 0])
    z = r0[2]*np.cos(omegaz2**0.5 * t)
    Aplus = (v0[1] + omegaminus*r0[0])/(omegaminus - omegaplus)
    Aminus = -(v0[1] + omegaplus*r0[0])/(omegaminus - omegaplus)
    x = Aplus * np.cos(omegaplus*t) + Aminus * np.cos(omegaminus*t)
    y = - Aplus * np.sin(omegaplus*t) - Aminus*np.sin(omegaminus*t)
    return np.array([x, y, z]).T


def get_period(z, t):
    # rough numerical approach to finding maxima and thereby the period
    dz = z[1:] - z[:-1]
    eps = 50.*np.min(np.abs(dz))
    extrema = np.argwhere(np.abs(dz) <= eps)
    z_max = np.max(z)
    z_near_max = np.argwhere(z[:-1] - 0.9*z_max >= 0.)
    maxima = np.intersect1d(extrema, z_near_max)
    # clean up potential indices that are too close
    idx = np.argwhere(t[maxima][1:] - t[maxima][:-1] < 0.9*np.max(t[maxima][1:] - t[maxima][:-1] ) )
    maxima = np.delete(maxima,idx)
    periods = t[maxima][1:] - t[maxima][:-1]
    print("Period: ", np.mean(periods))

T = 50.
t = np.linspace(0., T, N)

fig1, ax1 = plt.subplots(figsize=(4,4))
fig2, (ax_x, ax_y, ax_z) = plt.subplots(figsize=(4,6), nrows=3)
fig3 = plt.figure(figsize=(5,4))
ax3 = fig3.gca(projection="3d")
for i in range(nParticles):
    x = r[:,i,0]
    y = r[:,i,1]
    z = r[:,i,2]
    get_period(z,t) # print period of axial oscillationz
    ax1.plot(x,y)
    ax_x.plot(t,x,label="Simulated")
    ax_y.plot(t,y,label="Simulated")
    ax_z.plot(t,z,label="Simulated")
    ax3.plot(x,y,z,label="Simulated")
    if i==0:
        # plotting analytic solution
        r_anal = analytic(t)
        x_anal = r_anal[:,0]
        y_anal = r_anal[:,1]
        z_anal = r_anal[:,2]
        ax_x.plot(t,x_anal,linestyle="--",label="Analytic")
        ax_y.plot(t,y_anal,linestyle="--",label="Analytic")
        ax_z.plot(t,z_anal,linestyle="--",label="Analytic")
        ax3.plot(x_anal,y_anal,z_anal,linestyle="--",label="Analytic")

ax1.axis("equal")
ax1.set_xlabel(u"$x$ [\u03bcm]")
ax1.set_ylabel(u"$y$ [\u03bcm]")
ax_x.set_xticks([])
ax_x.set_ylabel(u"$x$ [\u03bcm]")
ax_y.set_xticks([])
ax_y.set_ylabel(u"$y$ [\u03bcm]")
ax_z.set_xlabel(u"Time [\u03bcs]")
ax_z.set_ylabel(u"$z$ [\u03bcm]")
ax3.set_xlabel(u"$x$ [\u03bcm]")
ax3.set_ylabel(u"$y$ [\u03bcm]")
ax3.set_zlabel(u"$z$ [\u03bcm]")

for fig in [fig1, fig2, fig3]:
    fig.tight_layout()
for ax in [ax_x, ax_y, ax_z, ax3]:
    ax.legend()

# calculating theoretical values
pi = 3.141593
B0 = 9.65e1
V0 = 2.41e6
d = 500.
q = 1
m = 40.08
omega0 = q*B0/m # cyclotron frequency
omegaz2 = 2*q*V0/(m*d**2) # axial frequency
omegaplus = 0.5*( omega0 + np.sqrt(omega0**2 - 2*omegaz2) ) # reduced cyclotron frequency
omegaminus = 0.5*( omega0 - np.sqrt(omega0**2 - 2*omegaz2) ) # magneton frequency
print("Theoretically:")
print("Axial frequency omegaz =", omegaz2**0.5, "with period =", 2*pi/(omegaz2**0.5))
print("Cyclotron frequency omega0 =", omega0, "with period =", 2*pi/omega0)
print("Reduced cyclotron frequency omega+ =", omegaplus, "with period =", 2*pi/omegaplus)
print("Magneton frequency omega- =", omegaminus, "with period =", 2*pi/omegaminus)

plt.show()