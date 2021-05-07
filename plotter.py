"""
Plotting function
"""
from create_geometry import BladeGeometry, doubleRotor
from lifting_line_solver import LiftingLineSolver
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

nspan = 25
ntheta = 25
nblades = 3
nrotor = 2

solver = BladeGeometry(50, 8, 10, 3, nspan, ntheta,spacing = 'equal')
cp = solver.cp
solver.compute_ring()
solver.discretize_blade()
blade = solver.bladepanels
rings = solver.filaments
prop_geo = BladeGeometry(radius=50.0, tsr=8, v_inf=10.0, n_blades=3, n_span=nspan, n_theta=ntheta, spacing='cosine')
solver = LiftingLineSolver(geo=prop_geo, u_rot=10/50, r_rotor=50, tol=1e-3, n_iter=100)
data=solver.run_solver()

# =============================================================================
# Double Rotor Plotting
# =============================================================================

def plotDoubleRotor(xshift,yshift,zshift):
    prop_geo = BladeGeometry(radius=50.0, tsr=8, v_inf=10.0, n_blades=3, n_span=nspan, n_theta=ntheta, spacing='cosine')
    doubleRotor(prop_geo,xshift,yshift,zshift)
    rings = prop_geo.filaments
       
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    
    # ax.set_xlim3d(0, 50)
    # ax.set_ylim3d(-50, 200)
    # ax.set_zlim3d(-50, 50)
    
    c=['green','blue','red','green','blue','red']
    for idx in range(nblades*nrotor):
        ax.plot_wireframe(rings[0, idx*(nspan-1):(idx+1)*(nspan-1), :ntheta+1],
                      rings[1, idx*(nspan-1):(idx+1)*(nspan-1), :ntheta+1],
                      rings[2, idx*(nspan-1):(idx+1)*(nspan-1), :ntheta+1],
                      color=c[idx], cstride=0)

plotDoubleRotor(100,100,100)

# =============================================================================
# Rotor performance plots
# =============================================================================
def plottingFunction(solver,prop_geo,data):
    plt.plot(prop_geo.centerPoints,data[-1][:len(prop_geo.centerPoints)])
    plt.show()
    return

# plottingFunction(solver,prop_geo,data)

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
# fig = plt.figure(dpi=150)
# ax1 = fig.add_subplot(131, projection="3d")
# ax2 = fig.add_subplot(132, projection="3d")
# ax3 = fig.add_subplot(133, projection="3d")
# axes = [ax1, ax2, ax3]

# split blade panels
# blades_split = np.split(blade, 3)

# for blade_i, ax in zip(blades_split, axes):
#     x_le = blade_i[:, 0]
#     y_le = blade_i[:, 1]
#     z_le = blade_i[:, 2]
#     x_te = blade_i[:, 9]
#     y_te = blade_i[:, 10]
#     z_te = blade_i[:, 11]
#     x = np.vstack((x_le, x_te))
#     y = np.vstack((y_le, y_te))
#     z = np.vstack((z_le, z_te))
#     ax.view_init(10, 180)
#     ax.set_xlabel("x")
#     ax.set_ylabel("y")
#     ax.set_zlabel("z")
#     ax.plot_wireframe(x, y, z)
# ax1.plot(solver.filaments[0,0,:],solver.filaments[1,0,:],solver.filaments[2,0,:])

# =============================================================================
# plotjes
# =============================================================================

# fig = plt.figure()
# ax = Axes3D(fig)
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')

# ax.plot_wireframe(rings[0, :nspan-1, :ntheta+1],
#                   rings[1, :nspan-1, :ntheta+1],
#                   rings[2, :nspan-1, :ntheta+1],
#                   color='green', cstride=0)
# ax.plot_wireframe(rings[0, nspan-1:2*nspan-2, :ntheta+1],
#                   rings[1, nspan-1:2*nspan-2, :ntheta+1],
#                   rings[2, nspan-1:2*nspan-2, :ntheta+1],
#                   color='blue', cstride=0)
# ax.plot_wireframe(rings[0, (2*nspan-2):, :ntheta+1],
#                   rings[1, (2*nspan-2):, :ntheta+1],
#                   rings[2, (2*nspan-2):, :ntheta+1],
#                   color='red', cstride=0)