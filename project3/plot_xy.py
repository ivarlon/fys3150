import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

data = pd.read_csv('xy_list.txt', sep = ' ', header = None).values

nSteps = data.shape[0]
nParticles = data.shape[1]//6
print(nParticles)
#print(data)
r = np.empty(shape=(nSteps, nParticles, 3))
v = np.empty(shape=(nSteps, nParticles, 3))
for i in range(nSteps):
    for j in range(nParticles):
        for k in range(3):
            r[i, j, k] = data[i][j*6 + k]
            v[i, j, k] = data[i][j*6 + 3 + k]


T = 50.
t = np.linspace(0., T, nSteps)
#fig, ax = plt.subplots(figsize=(5,5))
#ax.plot(t, z, color="black", label="z")
#ax.legend()
#ax.set_xlabel(u"Time [\u03bcs]"); ax.set_ylabel("z")
fig1, ax1 = plt.subplots()
fig2, (ax_x, ax_y, ax_z) = plt.subplots(ncols=3)
fig3 = plt.figure()
ax3 = fig3.gca(projection="3d")
for i in range(nParticles//10):
    x = r[:,i,0]
    y = r[:,i,1]
    z = r[:,i,2]
    ax1.plot(x,y)
    ax1.scatter([x[0],x[-1]], [y[0],y[-1]])
    ax1.axis("equal")
    ax_x.plot(t,x)
    ax_y.plot(t,y)
    ax_z.plot(t,z)
#surf = ax3.plot(x,y,z)

def analytic(t):
    B0 = 9.65e1
    V0 = 2.41e6
    d = 500.
    q = 10
    m = 20.
    omega0 = q*B0/m
    omegaz2 = 2*q*V0/(m*d**2)
    print(omega0, omegaz2**0.5)
    omegaplus = 0.5*( omega0 + np.sqrt(omega0**2 - 2*omegaz2) )
    omegaminus = 0.5*( omega0 - np.sqrt(omega0**2 - 2*omegaz2) )
    #r0 = np.array([0, 25, 0])
    #v0 = np.array([0, 40, 5])
    print(omegaplus,omegaminus)
    r0 = np.array([20, 0, 20])
    v0 = np.array([0, 25, 0])
    z = r0[2]*np.cos(omegaz2**0.5 * t)
    Aplus = (v0[1] + omegaminus*r0[0])/(omegaminus - omegaplus)
    Aminus = -(v0[1] + omegaplus*r0[0])/(omegaminus - omegaplus)
    print("A", Aplus, Aminus)
    x = Aplus * np.cos(omegaplus*t) + Aminus * np.cos(omegaminus*t)
    y = - Aplus * np.sin(omegaplus*t) - Aminus*np.sin(omegaminus*t)
    return x, y, z

t = np.linspace(0., T, nSteps)
#x, y, z = analytic(t)
#ax1.plot(x,y, label="anal")
#ax1.legend()

#plt.plot(x2,y2)
#plt.scatter([x2[0],x2[-1]], [y2[0],y2[-1]],color="orange")
plt.show()
#plt.plot(t, x1, color="black")
#plt.show()
#idx = np.where(np.abs(z - np.roll(z, shift=1)) < 1e-1)
#print(idx[0])
#ax.plot(t[idx[0]], z[idx[0]], color="red", label="z-z")
#plt.show()