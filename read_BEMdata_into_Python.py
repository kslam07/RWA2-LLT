
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

# Radial distribution alpha and phi

plt.figure()
plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_alpha, 0), BEM_rR.shape)[0, :]*180/np.pi, label=r'$\alpha$')
plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_phi, 0), BEM_rR.shape)[0, :]*180/np.pi, '--', label=r'$\phi$')
plt.xlabel('r/R (-)')
plt.ylabel('angle (deg)')
plt.legend()

# Radial distribution F_tan en F_ax

plt.figure()
plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_Ax[:,0]*BEM_rho, 0), BEM_rR.shape)[0, :],'--', label=r'F_{ax}')
plt.plot(BEM_rR[0, :], np.resize(np.mean(BEM_Az[:,0]*BEM_rho, 0), BEM_rR.shape)[0, :],'--', label=r'F_{tan}')
plt.xlabel('r/R (-)')
plt.ylabel('F (N)')
plt.legend()

# Radial distribution circulation

# figure()
# # made nondimensional with (np.pi*Uinf**2/(NBlades*Omega)
# plot(solver.rR,solver.gamma)
# xlabel('r/R (-)')
# ylabel('\Gamma (-)')

# Radial distribution CT

# plot(solver6.rR, mean(solver6.CT,2), solver6.rR, mean(solver6.Cq,2),'--', solver6.rR, mean(solver6.CN,2),'-.')
# plot(solver8.rR, mean(solver8.CT,2), solver8.rR, mean(solver8.Cq,2),'--', solver8.rR, mean(solver8.CN,2),'-.')
# plot(solver10.rR, mean(solver10.CT,2), solver10.rR, mean(solver10.Cq,2),'--', solver10.rR, mean(solver10.CN,2),'-.')
# xlabel('r/R (-)')
# ylabel('C (-)')
# legend('C_T_{,TSR=6}', 'C_Q_{,TSR=6}', 'C_N_{,TSR=6}','C_T_{,TSR=8}', 'C_Q_{,TSR=8}', 'C_N_{,TSR=8}','C_T_{,TSR=10}', 'C_Q_{,TSR=10}', 'C_N_{,TSR=10}','Location','eastoutside')

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
print('Done')