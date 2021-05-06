"""
Plotting function
"""
from create_geometry import BladeGeometry
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

nspan = 5
ntheta = 5

solver = BladeGeometry(50, 8, 10, 3, nspan, ntheta)
cp = solver.cp
solver.compute_ring()
solver.discretize_blade()
blade = solver.bladepanels
rings = solver.filaments

import random

# fig = plt.figure()
# ax = Axes3D(fig)


# ax.scatter(rings[0,:24,1:26],rings[1,:24,1:26],rings[2,:24,1:26])
# ax.scatter(rings[3,:24,1:26],rings[4,:24,1:26],rings[5,:24,1:26])
# ax.scatter(rings[0,24:48],rings[1,:24],rings[2,:24])
# ax.scatter(rings[0,24:48],rings[1,:24],rings[2,:24])
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlaebl('z')
# plt.show()

############################################### PLOT BLADE PANELS ######################################################
fig = plt.figure(dpi=150)
ax1 = fig.add_subplot(131, projection="3d")
ax2 = fig.add_subplot(132, projection="3d")
ax3 = fig.add_subplot(133, projection="3d")
axes = [ax1, ax2, ax3]

# split blade panels
blades_split = np.split(blade, 3)

for blade_i, ax in zip(blades_split, axes):
    x_le = blade_i[:, 0]
    y_le = blade_i[:, 1]
    z_le = blade_i[:, 2]
    x_te = blade_i[:, 9]
    y_te = blade_i[:, 10]
    z_te = blade_i[:, 11]
    x = np.vstack((x_le, x_te))
    y = np.vstack((y_le, y_te))
    z = np.vstack((z_le, z_te))
    ax.view_init(10, 180)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.plot_wireframe(x, y, z)
ax1.plot(solver.filaments[0,0,:],solver.filaments[1,0,:],solver.filaments[2,0,:])

fig = plt.figure()
ax = Axes3D(fig)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

ax.plot_wireframe(rings[0, :nspan-2, :ntheta-1],rings[1, :nspan-2, :ntheta-1],rings[2, :nspan-2, :ntheta-1], color='green')
ax.plot_wireframe(rings[0, nspan-1:(2*nspan)-3, :ntheta-1],rings[1, nspan-1:(2*nspan)-3, :ntheta-1],rings[2, nspan-1:(2*nspan)-3, :ntheta-1], color='blue')
ax.plot_wireframe(rings[0, (2*nspan)-2:(3*nspan)-4, :ntheta-1],rings[1, (2*nspan)-2:(3*nspan)-4, :ntheta-1],rings[2, (2*nspan)-2:(3*nspan)-4, :ntheta-1], color='red')

plt.show()