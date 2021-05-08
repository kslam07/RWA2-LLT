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


# =============================================================================
# BLADE VIZ
# =============================================================================
fig = plt.figure(dpi=150)
ax1 = fig.add_subplot(111, projection="3d")
# ax2 = fig.add_subplot(132, projection="3d")
# ax3 = fig.add_subplot(133, projection="3d")
# axes = [ax1, ax2, ax3]

# split blade panels
blades_split = np.split(blade, 3)

for blade_i in (blades_split):
    x = blade_i[:, [0, 9]]
    y = blade_i[:, [1, 10]]
    z = blade_i[:, [2, 11]]

    ax1.plot_surface(x, y, z, color='k')

# =============================================================================
# Draw rotor hub
# =============================================================================
theta = np.linspace(0, 2 * np.pi)
y = solver.geo.span_arr[0]*np.cos(theta)*solver.r_rotor
z = solver.geo.span_arr[0]*np.sin(theta)*solver.r_rotor

ax1.plot(y*0, y, z, color='k')

# =============================================================================
# WAKE VIZ
# =============================================================================

ax1.plot_wireframe(rings[0, :nspan-1, :ntheta+1],
                  rings[1, :nspan-1, :ntheta+1],
                  rings[2, :nspan-1, :ntheta+1],
                  color='green', cstride=0)
ax1.plot_wireframe(rings[0, nspan-1:2*nspan-2, :ntheta+1],
                  rings[1, nspan-1:2*nspan-2, :ntheta+1],
                  rings[2, nspan-1:2*nspan-2, :ntheta+1],
                  color='blue', cstride=0)
ax1.plot_wireframe(rings[0, (2*nspan-2):, :ntheta+1],
                  rings[1, (2*nspan-2):, :ntheta+1],
                  rings[2, (2*nspan-2):, :ntheta+1],
                  color='red', cstride=0)

ax1.view_init(0, 180)
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_zlabel("z")
ax1.set_xlim(-5, 5)
ax1.set_ylim(-50, 50)
ax1.set_zlim(-50, 50)
ax1.axis("equal")