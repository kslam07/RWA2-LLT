"""
Sensitivity Plots
"""
from create_geometry import BladeGeometry
from lifting_line_solver import LiftingLineSolver
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from read_BEMdata_into_Python import read_matlab_data
import numpy as np

plot_Spacing = False
plot_DiscrAzimuthal = False

# =============================================================================
#  Wake Discretization
# =============================================================================
if plot_DiscrAzimuthal:
    nspan = 15
    nthetaList = [50,100,400]
    nthetaList = [25,50,75]
    nblades = 3
    spacing = 'equal'
    nrotor = 2
    
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
        plt.plot(data[2][:nspan - 1, 0], np.resize(data[5], data[2].shape)[:nspan - 1, 0] / circ_nondim, '--', label='$N_{theta}$= '+str(ntheta))
    plt.xlabel('Radial location r/R (-)', fontsize=15)
    plt.ylabel(r'Circulation $\Gamma$ (-)', fontsize=15)
    plt.title(r'Radial distribution of $\Gamma$', fontsize=16)
    plt.legend(fontsize=15)
    plt.grid(True)