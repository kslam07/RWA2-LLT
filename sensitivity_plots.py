"""
Sensitivity Plots
"""
from create_geometry import BladeGeometry
from create_geometry_varwake import BladeGeometry2
from lifting_line_solver import LiftingLineSolver
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from read_BEMdata_into_Python import read_matlab_data
import numpy as np

plot_Spacing = False
plot_DiscrAzimuthal = False
plot_ConvecSpeed = False
plot_WakeLength = False

# =============================================================================
#  Wake Discretization
# =============================================================================
if plot_DiscrAzimuthal:

    nspan = 15
    nthetaList = [50,100,400]
    nthetaList = [25,50,75]
    nblades = 3
    spacing = 'equal'
    
    plt.figure(figsize=(8, 6), dpi=150)
    for ntheta in nthetaList:
        prop_geo = BladeGeometry(radius=50.0, tsr=8, v_inf=10.0, n_blades=3, n_span=nspan, n_theta=ntheta,
            spacing=spacing, a=0, xshift=0, yshift=100, zshift=0, phase_diff=0, double_rotor=False)
        # prop_geo.doubleRotor()
        solver = LiftingLineSolver(geo=prop_geo, r_rotor=50, weight=0.5, tol=1e-6,
                                   n_iter=200)
        
        data = solver.run_solver()
        omega = solver.geo.tsr * solver.geo.v_inf / solver.geo.radius
        
        circ_nondim = (np.pi * solver.geo.v_inf ** 2) / (nblades * omega)
        
        
        # plt.plot(BEM_rR[0, :], BEM_Gamma[0, :], '-b', label=r'$\Gamma$ BEM')
        plt.plot(data[2][:nspan - 1, 0], np.resize(data[5], data[2].shape)[:nspan - 1, 0] / circ_nondim, '-', label='$N_{theta}$= '+str(ntheta))

    plt.xlabel('Radial location r/R (-)', fontsize=15)
    plt.ylabel(r'Circulation $\Gamma$ (-)', fontsize=15)
    plt.title(r'Radial distribution of $\Gamma$', fontsize=16)
    plt.legend(fontsize=15)
    plt.grid(True)

if plot_ConvecSpeed:

    nspan = 20
    ntheta = 200
    nblades = 3
    spacing = 'equal'

    tsr_list = [4, 6, 8]
    circ_list = []
    Fax_list = []
    Faz_list = []

    for tsr in tsr_list:

        prop_geo = BladeGeometry(radius=50.0, tsr=tsr, v_inf=10.0, n_blades=3, n_span=nspan, n_theta=ntheta,
                                 spacing=spacing, a=0, xshift=0, yshift=100, zshift=0, phase_diff=0, double_rotor=False)
        # prop_geo.doubleRotor()
        solver = LiftingLineSolver(geo=prop_geo, r_rotor=50, weight=0.5, tol=1e-6, n_iter=200)

        data = solver.run_solver()
        omega = solver.geo.tsr * solver.geo.v_inf / solver.geo.radius

        circ_nondim = (np.pi * solver.geo.v_inf ** 2) / (nblades * omega)
        circ_list.append(np.resize(data[5], data[2].shape)[:nspan - 1, 0] / circ_nondim)

        F_nondim_LLM = 0.5 * (solver.geo.v_inf ** 2) * solver.geo.radius
        Fax_list.append(np.resize(data[3], data[2].shape)[:nspan - 1, 0] / F_nondim_LLM)
        Faz_list.append(np.resize(data[4], data[2].shape)[:nspan - 1, 0] / F_nondim_LLM)

    plt.figure(figsize=(8, 6), dpi=150)
    plt.plot(data[2][:nspan - 1, 0], circ_list[0], '-r', label='$\lambda$ = ' + str(tsr_list[0]))
    plt.plot(data[2][:nspan - 1, 0], circ_list[1], '-b', label='$\lambda$ = ' + str(tsr_list[1]))
    plt.plot(data[2][:nspan - 1, 0], circ_list[2], '-g', label='$\lambda$ = ' + str(tsr_list[2]))
    plt.xlabel('Radial location r/R (-)', fontsize=15)
    plt.ylabel(r'Circulation $\Gamma$ (-)', fontsize=15)
    plt.title(r'Radial distribution of $\Gamma$', fontsize=16)
    plt.legend(fontsize=15)
    plt.grid(True)
    plt.savefig('Sensitivity_figures/convec_circ.eps', format='eps')

    plt.figure(figsize=(8, 6), dpi=150)
    plt.plot(data[2][:nspan - 1, 0], Fax_list[0], '-r', label=r'$F_{ax}$ | $\lambda$ = ' + str(tsr_list[0]))
    plt.plot(data[2][:nspan - 1, 0], Fax_list[1], '-b', label=r'$F_{ax}$ | $\lambda$ = ' + str(tsr_list[1]))
    plt.plot(data[2][:nspan - 1, 0], Fax_list[2], '-g', label=r'$F_{ax}$ | $\lambda$ = ' + str(tsr_list[2]))
    plt.plot(data[2][:nspan - 1, 0], Faz_list[0], '--r', label=r'$F_{az}$ | $\lambda$ = ' + str(tsr_list[0]))
    plt.plot(data[2][:nspan - 1, 0], Faz_list[1], '--b', label=r'$F_{az}$ | $\lambda$ = ' + str(tsr_list[1]))
    plt.plot(data[2][:nspan - 1, 0], Faz_list[2], '--g', label=r'$F_{az}$ | $\lambda$ = ' + str(tsr_list[2]))
    plt.xlabel('Radial location r/R (-)', fontsize=15)
    plt.ylabel('Force $F$ (-)', fontsize=15)
    plt.title(r'Radial distribution of $F_{ax}$ and $F_{az}$', fontsize=16)
    plt.legend(fontsize=15)
    plt.grid(True)
    plt.savefig('Sensitivity_figures/convec_forces.eps', format='eps')

if plot_WakeLength:

    nspan = 20
    ntheta = 200
    nblades = 3
    spacing = 'equal'

    N_rotations = [5, 10, 15]
    circ_list = []
    Fax_list = []
    Faz_list = []
    error_list = []

    fig = plt.figure(dpi=150)
    i = 0

    for rotation in N_rotations:

        prop_geo = BladeGeometry2(radius=50.0, tsr=8, v_inf=10.0, n_blades=3, n_span=nspan, n_theta=ntheta,
                                  spacing=spacing, a=0, n_rotations=rotation, xshift=0, yshift=100, zshift=0,
                                  phase_diff=0, double_rotor=False)

        blade = prop_geo.bladepanels
        rings = prop_geo.filaments

        solver = LiftingLineSolver(geo=prop_geo, r_rotor=50, weight=0.5, tol=1e-6, n_iter=200)

        data = solver.run_solver()
        omega = solver.geo.tsr * solver.geo.v_inf / solver.geo.radius

        circ_nondim = (np.pi * solver.geo.v_inf ** 2) / (nblades * omega)
        circ_list.append(np.resize(data[5], data[2].shape)[:nspan - 1, 0] / circ_nondim)

        F_nondim_LLM = 0.5 * (solver.geo.v_inf ** 2) * solver.geo.radius
        Fax_list.append(np.resize(data[3], data[2].shape)[:nspan - 1, 0] / F_nondim_LLM)
        Faz_list.append(np.resize(data[4], data[2].shape)[:nspan - 1, 0] / F_nondim_LLM)
        error_list.append(data[8])

        plot_wake = False

        if plot_wake:

            i += 1

            ax = fig.add_subplot(1, 3, i, projection="3d")

            # =============================================================================
            # WAKE VIZ
            # =============================================================================

            ax.plot_wireframe(rings[0, :nspan - 1, :ntheta + 1]/(40),
                               rings[1, :nspan - 1, :ntheta + 1],
                               rings[2, :nspan - 1, :ntheta + 1],
                               color='green', ccount=0, rcount=5)
            ax.plot_wireframe(rings[0, nspan - 1:2 * nspan - 2, :ntheta + 1]/(40),
                               rings[1, nspan - 1:2 * nspan - 2, :ntheta + 1],
                               rings[2, nspan - 1:2 * nspan - 2, :ntheta + 1],
                               color='blue', ccount=0, rcount=5)
            ax.plot_wireframe(rings[0, (2 * nspan - 2):, :ntheta + 1]/(40),
                               rings[1, (2 * nspan - 2):, :ntheta + 1],
                               rings[2, (2 * nspan - 2):, :ntheta + 1],
                               color='red', ccount=0, rcount=5)

            # split blade panels
            blades_split = np.split(blade, 3)

            for blade_i in (blades_split):
                x = blade_i[:, [0, 9]]
                y = blade_i[:, [1, 10]]
                z = blade_i[:, [2, 11]]

                ax.plot_surface(x, y, z, color='k')

            # =============================================================================
            # Draw rotor hub
            # =============================================================================
            theta = np.linspace(0, 2 * np.pi)
            y = solver.geo.span_arr[0] * np.cos(theta) * solver.r_rotor
            z = solver.geo.span_arr[0] * np.sin(theta) * solver.r_rotor

            ax.plot(y * 0, y, z, color='k')

            ax.view_init(0, 180)
            ax.set_yticklabels([])
            ax.set_zticklabels([])
            ax.set_xlim(0, N_rotations[-1])
            ax.set_xlabel(r"Wakelength $\frac{L}{2 \pi}$ (-)")
            ax.set_title('{} Rotations'.format(rotation))
            ax.view_init(elev=10, azim=-450)

    plt.figure(figsize=(8, 6), dpi=150)
    plt.plot(data[2][:nspan - 1, 0], circ_list[0], '-r', label='# Rotations = ' + str(N_rotations[0]))
    plt.plot(data[2][:nspan - 1, 0], circ_list[1], '-b', label='# Rotations = ' + str(N_rotations[1]))
    plt.plot(data[2][:nspan - 1, 0], circ_list[2], '-g', label='# Rotations = ' + str(N_rotations[2]))
    plt.xlabel('Radial location r/R (-)', fontsize=15)
    plt.ylabel(r'Circulation $\Gamma$ (-)', fontsize=15)
    plt.title(r'Radial distribution of $\Gamma$', fontsize=16)
    plt.legend(fontsize=15)
    plt.grid(True)
    plt.savefig('Sensitivity_figures/varwake_circ.eps', format='eps')

    plt.figure(figsize=(8, 6), dpi=150)
    plt.plot(data[2][:nspan - 1, 0], Fax_list[0], '-r', label=r'$F_{ax}$ | # Rotations = ' + str(N_rotations[0]))
    plt.plot(data[2][:nspan - 1, 0], Fax_list[1], '-b', label=r'$F_{ax}$ | # Rotations = ' + str(N_rotations[1]))
    plt.plot(data[2][:nspan - 1, 0], Fax_list[2], '-g', label=r'$F_{ax}$ | # Rotations = ' + str(N_rotations[2]))
    plt.plot(data[2][:nspan - 1, 0], Faz_list[0], '--r', label=r'$F_{az}$ | # Rotations = ' + str(N_rotations[0]))
    plt.plot(data[2][:nspan - 1, 0], Faz_list[1], '--b', label=r'$F_{az}$ | # Rotations = ' + str(N_rotations[1]))
    plt.plot(data[2][:nspan - 1, 0], Faz_list[2], '--g', label=r'$F_{az}$ | # Rotations = ' + str(N_rotations[2]))
    plt.xlabel('Radial location r/R (-)', fontsize=15)
    plt.ylabel('Force $F$ (-)', fontsize=15)
    plt.title(r'Radial distribution of $F_{ax}$ and $F_{az}$', fontsize=16)
    plt.legend(fontsize=15)
    plt.grid(True)
    plt.savefig('Sensitivity_figures/varwake_forces.eps', format='eps')

    plt.figure(figsize=(8, 6), dpi=150)
    plt.semilogy(np.arange(1, len(error_list[0]) + 1), error_list[0], '-r', label='# Rotations = ' + str(N_rotations[0]))
    plt.semilogy(np.arange(1, len(error_list[1]) + 1), error_list[1], '-b', label='# Rotations = ' + str(N_rotations[1]))
    plt.semilogy(np.arange(1, len(error_list[2]) + 1), error_list[2], '-g', label='# Rotations = ' + str(N_rotations[2]))
    # plt.plot(np.arange(1, len(error_list[0]) + 1), error_list[0], '-r', label='# Rotations = ' + str(N_rotations[0]))
    # plt.plot(np.arange(1, len(error_list[1]) + 1), error_list[1], '-b', label='# Rotations = ' + str(N_rotations[1]))
    # plt.plot(np.arange(1, len(error_list[2]) + 1), error_list[2], '-g', label='# Rotations = ' + str(N_rotations[2]))
    plt.xlabel('Iterations N', fontsize=15)
    plt.ylabel(r'Log error $\epsilon$ (-)', fontsize=15)
    plt.title('Error', fontsize=16)
    plt.legend(fontsize=15)
    plt.grid(True)
    plt.savefig('Sensitivity_figures/varwake_convergence.eps', format='eps')

plt.show()