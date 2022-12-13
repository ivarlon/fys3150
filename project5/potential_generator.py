"""
FYS4150

Project 5

Generates potential matrix V for a given setup

"""
import matplotlib.pyplot as plt
import numpy as np

M = 101 # no. of points in x / y direction
h = 1./(M-1) # spatial step length
V0 = 1e10

V = np.zeros((M,M))
# ===================
# no potential
# ===================
# V = V[1:-1, 1:-1]
# np.savetxt("no_potential.dat", V)

# =========================
# double slit
# =========================
# thickness = 0.02
# x_position = 0.5
# aperture = 0.05
# distance_bw_slits = 0.05
# y_position = 0.5


# i_start = int((x_position - 0.5*thickness)/h)
# i_end = int((x_position + 0.5*thickness)/h)
# V[i_start+1:i_end+1,:] = V0
# print(i_end/M,i_start/M,h)
# j_start = int( (y_position - 0.5*distance_bw_slits - aperture)/h )
# j_end = int( (y_position - 0.5*distance_bw_slits)/h )
# V[:,j_start+1:j_end+1] = 0.

# print(j_start, j_end)

# j_start = int( (y_position + 0.5*distance_bw_slits)/h )
# j_end = int( (y_position + 0.5*distance_bw_slits + aperture)/h )
# V[:,j_start+1:j_end+1] = 0.

# print(j_start, j_end)

# V = V[1:-1,1:-1] # truncate boundaries

# plt.imshow(V)
# plt.show()

# np.savetxt("double_slit.dat", V)


# ===============================
# Single slit
# ===============================
thickness = 0.02
x_position = 0.5
aperture = 0.05
y_position = 0.5

i_start = int((x_position - 0.5*thickness)/h)
i_end = int((x_position + 0.5*thickness)/h)
V[i_start+1:i_end+1,:] = V0

j_start = int( (y_position - 0.5*aperture)/h )
j_end = int( (y_position + 0.5*aperture)/h )
V[:,j_start+1:j_end+1] = 0.

V = V[1:-1,1:-1] # truncate boundaries

np.savetxt("single_slit.dat", V)

# ===============================
# Triple slit
# ===============================
# thickness = 0.02
# x_position = 0.5
# aperture = 0.05
# distance_bw_slits = 0.05
# y_position = 0.5


# i_start = int((x_position - 0.5*thickness)/h)
# i_end = int((x_position + 0.5*thickness)/h)
# V[i_start+1:i_end+1,:] = V0

# for slit in [1,2,3]:
#     j_start = int( (y_position - 0.5*aperture + (slit-2)*(distance_bw_slits + aperture) )/h )
#     j_end = int( (y_position + 0.5*aperture + (slit-2)*(distance_bw_slits + aperture) )/h )
#     V[:,j_start+1:j_end+1] = 0.
#     print(j_start, j_end)

# V = V[1:-1,1:-1] # truncate boundaries

# plt.imshow(V)
# plt.show()

# np.savetxt("triple_slit.dat", V)