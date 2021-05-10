"""
Plotting function
"""
from create_geometry import BladeGeometry, doubleRotor
from lifting_line_solver import LiftingLineSolver
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from read_BEMdata_into_Python import read_matlab_data
import numpy as np

nspan = 100
ntheta = 100
nblades = 3
nrotor = 1

solver = BladeGeometry(50, 8, 10, 3, nspan, ntheta,spacing = 'equal')
cp = solver.cp
solver.compute_ring()
solver.discretize_blade()
blade = solver.bladepanels
rings = solver.filaments
prop_geo = BladeGeometry(radius=50.0, tsr=8, v_inf=10.0, n_blades=3, n_span=nspan, n_theta=ntheta, spacing='equal')
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
ax1.axis("auto")

# =============================================================================
# BEM COMPARISON
# =============================================================================
plt.close('All')

[BEM_rR, BEM_alpha, BEM_phi, BEM_rho, BEM_Ax, BEM_Az, BEM_Gamma , BEM_CT, BEM_CP] = read_matlab_data()

# Radial distribution alpha and phi

plt.figure()
plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_alpha, 0), BEM_rR.shape)[0, :]*180/np.pi, '-r', label=r'$\alpha$ BEM')
plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_phi, 0), BEM_rR.shape)[0, :]*180/np.pi, '-b', label=r'$\phi$ BEM')
plt.plot(data[2][:nspan-1, 0], np.degrees(np.resize(data[6], data[2].shape)[:nspan-1, 0]), '--r', label=r'$\alpha$ LLM')
plt.plot(data[2][:nspan-1, 0], np.degrees(np.resize(data[7], data[2].shape)[:nspan-1, 0]), '--b', label=r'$\phi$ LLM')
plt.xlabel('r/R (-)')
plt.ylabel('angle (deg)')
plt.legend()
plt.grid(True)

# Radial distribution F_tan en F_ax

plt.figure()
plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_Ax, 0), BEM_rR.shape)[0, :]*BEM_rho[0], '-r', label=r'$F_{ax}$ BEM')
plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_Az, 0), BEM_rR.shape)[0, :]*BEM_rho[0], '-b', label=r'$F_{tan}$ BEM')
# Plot one blade of LLM
plt.plot(data[2][:nspan-1, 0], np.resize(data[3], data[2].shape)[:nspan-1, 0]*BEM_rho[0], '--r', label=r'$F_{ax}$ LLM')
plt.plot(data[2][:nspan-1, 0], np.resize(data[4], data[2].shape)[:nspan-1, 0]*BEM_rho[0], '--b', label=r'$F_{tan}$ LLM')
plt.xlabel('r/R (-)')
plt.ylabel('F (N)')
plt.legend()
plt.grid(True)

# Radial distribution circulation

# made non-dimensional with (np.pi * Uinf**2) / (NBlades*Omega)
circ_nondim = (np.pi*solver.geo.v_inf**2)/(solver.geo.tsr*solver.geo.v_inf/solver.geo.radius)

plt.figure()
plt.plot(BEM_rR[0, :], BEM_Gamma[0, :], label=r'$\Gamma$ BEM')
plt.plot(data[2][:nspan-1, 0], np.resize(data[5], data[2].shape)[:nspan-1, 0]/circ_nondim, label=r'$\Gamma$ LLM')
plt.xlabel('r/R (-)')
plt.ylabel(r'$\Gamma$ (-)')
plt.legend()
plt.grid(True)

# Radial distribution CT

## STILL NEED TO DIVIDE BY SEGMENT AREA
area = 1
CT_LLM = np.resize(data[3], data[2].shape)[:, 0]/(0.5*area*solver.geo.v_inf**2)

plt.figure()
plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_CT, 0), BEM_rR.shape)[0, :], '-r', label=r'$C_T$ BEM')
plt.plot(data[2][:nspan-1, 0], CT_LLM[:nspan-1], '--r', label=r'$C_T$ LLM')
plt.xlabel('r/R (-)')
plt.ylabel(r'$C_T$ (-)')
plt.legend()
plt.grid(True)

# Radial distribution CP

CP_LLM = np.resize(data[3], data[2].shape)[:, 0]*np.resize(data[0], data[2].shape)[:, 0]\
         *data[2][:, 0]*solver.geo.radius*(solver.geo.tsr*solver.geo.v_inf/solver.geo.radius)\
         /(0.5*(solver.geo.v_inf**3)*np.pi*solver.geo.radius**2)

plt.figure()
plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_CP, 0), BEM_rR.shape)[0, :], '-r', label=r'$C_P$ BEM')
plt.plot(data[2][:nspan-1, 0], CP_LLM[:nspan-1], '--r', label=r'$C_P$ LLM')
plt.xlabel('r/R (-)')
plt.ylabel('$C_P$ (-)')
plt.legend()
plt.grid(True)

plt.show()