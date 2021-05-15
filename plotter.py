"""
Plotting function
"""
from create_geometry import BladeGeometry
from lifting_line_solver import LiftingLineSolver
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from read_BEMdata_into_Python import read_matlab_data
import numpy as np

nspan = 20
ntheta = 200
nblades = 3
spacing = 'equal'
nrotor = 2

prop_geo = BladeGeometry(radius=50.0, tsr=8, v_inf=10.0, n_blades=3, n_span=nspan, n_theta=ntheta, spacing=spacing, a=0)

blade = prop_geo.bladepanels
rings = prop_geo.filaments
solver = LiftingLineSolver(geo=prop_geo, r_rotor=50, weight=0.5, tol=1e-8, n_iter=100, double_rotor=False)

data = solver.run_solver()
omega = solver.geo.tsr * solver.geo.v_inf / solver.geo.radius
[CP_LLM, CT_LLM] = solver.CP_and_CT(np.resize(data[0], data[2].shape), np.resize(data[1], data[2].shape), data[2],
                                    np.resize(data[3], data[2].shape), np.resize(data[4], data[2].shape),
                                    solver.geo.v_inf, omega, solver.geo.radius, nblades)


# =============================================================================
# Rotor performance plots
# =============================================================================
def plottingFunction(prop_geo, data):
    plt.plot(prop_geo.centerPoints, data[-1][:len(prop_geo.centerPoints)])
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
y = solver.geo.span_arr[0] * np.cos(theta) * solver.r_rotor
z = solver.geo.span_arr[0] * np.sin(theta) * solver.r_rotor

ax1.plot(y * 0, y, z, color='k')

# =============================================================================
# WAKE VIZ
# =============================================================================

ax1.plot_wireframe(rings[0, :nspan - 1, :ntheta + 1],
                   rings[1, :nspan - 1, :ntheta + 1],
                   rings[2, :nspan - 1, :ntheta + 1],
                   color='green', cstride=0)
ax1.plot_wireframe(rings[0, nspan - 1:2 * nspan - 2, :ntheta + 1],
                   rings[1, nspan - 1:2 * nspan - 2, :ntheta + 1],
                   rings[2, nspan - 1:2 * nspan - 2, :ntheta + 1],
                   color='blue', cstride=0)
ax1.plot_wireframe(rings[0, (2 * nspan - 2):, :ntheta + 1],
                   rings[1, (2 * nspan - 2):, :ntheta + 1],
                   rings[2, (2 * nspan - 2):, :ntheta + 1],
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

[BEM_rR, BEM_alpha, BEM_phi, BEM_rho, BEM_Ax, BEM_Az, BEM_Gamma, BEM_CT, BEM_CP, BEM_a, BEM_aline] = read_matlab_data()

# Radial distribution alpha and phi

plt.figure()
plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_alpha, 0), BEM_rR.shape)[0, :] * 180 / np.pi, '-r', label=r'$\alpha$ BEM')
plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_phi, 0), BEM_rR.shape)[0, :] * 180 / np.pi, '-b', label=r'$\phi$ BEM')
plt.plot(data[2][:nspan - 1, 0], np.degrees(np.resize(data[6], data[2].shape)[:nspan - 1, 0]), '--r',
         label=r'$\alpha$ LLM')
plt.plot(data[2][:nspan - 1, 0], np.degrees(np.resize(data[7], data[2].shape)[:nspan - 1, 0]), '--b',
         label=r'$\phi$ LLM')
plt.xlabel('r/R (-)')
plt.ylabel('angle (deg)')
plt.legend()
plt.grid(True)

# Radial distribution F_tan en F_ax

plt.figure()
plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_Ax, 0), BEM_rR.shape)[0, :] * BEM_rho[0], '-r', label=r'$F_{ax}$ BEM')
plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_Az, 0), BEM_rR.shape)[0, :] * BEM_rho[0], '-b', label=r'$F_{tan}$ BEM')
# Plot one blade of LLM
plt.plot(data[2][:nspan - 1, 0], np.resize(data[3], data[2].shape)[:nspan - 1, 0] * BEM_rho[0], '--r',
         label=r'$F_{ax}$ LLM')
plt.plot(data[2][:nspan - 1, 0], np.resize(data[4], data[2].shape)[:nspan - 1, 0] * BEM_rho[0], '--b',
         label=r'$F_{tan}$ LLM')
plt.xlabel('r/R (-)')
plt.ylabel('F (N)')
plt.legend()
plt.grid(True)

# Radial distribution circulation

# made non-dimensional with (np.pi * Uinf**2) / (NBlades*Omega)
circ_nondim = (np.pi * solver.geo.v_inf ** 2) / (nblades * omega)

plt.figure()
plt.plot(BEM_rR[0, :], BEM_Gamma[0, :], label=r'$\Gamma$ BEM')
plt.plot(data[2][:nspan - 1, 0], np.resize(data[5], data[2].shape)[:nspan - 1, 0] / circ_nondim, label=r'$\Gamma$ LLM')
plt.xlabel('r/R (-)')
plt.ylabel(r'$\Gamma$ (-)')
plt.legend()
plt.grid(True)

# INDUCTION FACTORS
fig, ax = plt.subplots(1, 2, dpi=150)
ax[0].plot(BEM_rR[0, :], BEM_a[0, :], label="BEM")
ax[0].plot(data[2][:nspan - 1], data[0][:nspan - 1], label="Lifting Line")
ax[1].plot(BEM_rR[0, :], BEM_aline[0, :], label="BEM")
ax[1].plot(data[2][:nspan - 1], data[1][:nspan - 1], label="Lifting Line")
ax[0].set_xlabel("r/R [-]")
ax[0].set_ylabel("a [-]")
ax[1].set_xlabel("r/R [-]")
ax[1].set_ylabel("a' [-]")
ax[0].grid()
ax[1].grid()
ax[0].legend()

# RADIAL DISTRIBUTION CT

CT_LLM2 = np.resize(data[3], data[2].shape)[:, 0] / (0.5 * np.pi * (solver.geo.radius ** 2) * solver.geo.v_inf ** 2)
print('CT Carlos:', np.sum(CT_LLM))
print('CT:', np.sum(CT_LLM2))

# plt.figure()
# plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_CT, 0), BEM_rR.shape)[0, :], '-r', label=r'$C_T$ BEM')
# plt.plot(data[2][:nspan-1, 0], CT_LLM[:nspan-1], '--r', label=r'$C_T$ LLM Carlos')
# plt.plot(data[2][:nspan-1, 0], CT_LLM2[:nspan-1], '--g', label=r'$C_T$ LLM 2')
# plt.xlabel('r/R (-)')
# plt.ylabel(r'$C_T$ (-)')
# plt.legend()
# plt.grid(True)

# Radial distribution CP

CP_LLM2 = np.resize(data[3], data[2].shape)[:, 0]*np.resize(data[0], data[2].shape)[:, 0]\
          *data[2][:, 0]*solver.geo.radius*omega/(0.5*(solver.geo.v_inf**3)*np.pi*solver.geo.radius**2)

print('CP Carlos:', np.sum(CP_LLM))
print('CP:', np.sum(CP_LLM2))

# plt.figure()
# plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_CP, 0), BEM_rR.shape)[0, :], '-r', label=r'$C_P$ BEM')
# plt.plot(data[2][:nspan-2, 0], CP_LLM[:nspan-2], '--r', label=r'$C_P$ LLM Carlos')
# plt.plot(data[2][:nspan-1, 0], CP_LLM2[:nspan-1], '--g', label=r'$C_P$ LLM 2')
# plt.xlabel('r/R (-)')
# plt.ylabel('$C_P$ (-)')
# plt.legend()
# plt.grid(True)
