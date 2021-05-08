
import matplotlib.pyplot as plt
import numpy as np
import h5py

# Read .mat file
f = h5py.File('BEM_data.mat','r')
# print(f.keys())

# Load the BEM variables that you need
BEM_rR = np.asarray(f.get('rR'))
BEM_alpha = np.asarray(f.get('alpha'))
BEM_phi = np.asarray(f.get('phi'))
BEM_rho = np.asarray(f.get('rho'))
BEM_Ax = np.asarray(f.get('Ax'))
BEM_Az = np.asarray(f.get('Az'))
BEM_Gamma = np.asarray(f.get('gamma'))
BEM_CT = np.asarray(f.get('CT'))
BEM_CP = np.asarray(f.get('CP'))

# Close .mat file
f.close()

plt.close('All')

# Radial distribution alpha and phi

plt.figure()
plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_alpha, 0), BEM_rR.shape)[0, :]*180/np.pi, label=r'$\alpha$')
plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_phi, 0), BEM_rR.shape)[0, :]*180/np.pi, '--', label=r'$\phi$')
plt.xlabel('r/R (-)')
plt.ylabel('angle (deg)')
plt.legend()
plt.grid(True)

# Radial distribution F_tan en F_ax

plt.figure()
plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_Ax, 0), BEM_rR.shape)[0, :]*BEM_rho[0], label=r'F_{ax}')
plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_Az, 0), BEM_rR.shape)[0, :]*BEM_rho[0],'--', label=r'F_{tan}')
plt.xlabel('r/R (-)')
plt.ylabel('F (N)')
plt.legend()
plt.grid(True)

# Radial distribution circulation

plt.figure()
# made non-dimensional with (np.pi * Uinf**2) / (NBlades*Omega)
plt.plot(BEM_rR[0, :], BEM_Gamma[0, :], label=r'$\Gamma$')
plt.xlabel('r/R (-)')
plt.ylabel(r'$\Gamma$ (-)')
plt.legend()
plt.grid(True)

# Radial distribution CT

plt.figure()
plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_CT, 0), BEM_rR.shape)[0, :], label=r'$C_T$')
plt.xlabel('r/R (-)')
plt.ylabel(r'$C_T$ (-)')
plt.legend()
plt.grid(True)

# Radial distribution CP

plt.figure()
plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_CP, 0), BEM_rR.shape)[0, :], label=r'$C_P$')
plt.xlabel('r/R (-)')
plt.ylabel('$C_P$ (-)')
plt.legend()
plt.grid(True)

# Polar plot

# minVal = min([solverS0.alpha(:); solverS15.alpha(:); solverS30.alpha(:)]);
# maxVal = max([solverS0.alpha(:); solverS15.alpha(:); solverS30.alpha(:)]);
# t1 = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
# nexttile
# pplot = contourf(x, y, rad2deg(solverS0.alpha));
# h=colorbar;
# caxis manual
# caxis([rad2deg(minVal), rad2deg(maxVal)])
# xlabel(h,'\alpha','Rotation',0,'FontSize',18)
# xlabel('x/R (-)', 'FontSize',17)
# ylabel('y/R (-)', 'FontSize',17)
# colormap('default')

plt.show()
