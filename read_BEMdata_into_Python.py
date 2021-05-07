
import numpy as np
import h5py

# Read .mat file
f = h5py.File('BEM_data.mat','r')
print(f.keys())

# Load BEM variables
BEM_rR = np.asarray(f.get('rR'))
BEM_alpha = np.asarray(f.get('alpha'))
BEM_phi = np.asarray(f.get('alpha'))

# Close .mat file
f.close()

# Radial distribution alpha

# figure()
# plot(solver6.rR, mean(solver6.alpha,2)*180/pi, solver6.rR, mean(solver6.phi,2)*180/pi,'--');
# plot(solver8.rR, mean(solver8.alpha,2)*180/pi, solver8.rR, mean(solver8.phi,2)*180/pi,'--');
# plot(solver10.rR, mean(solver10.alpha,2)*180/pi, solver10.rR, mean(solver10.phi,2)*180/pi,'--');
# xlabel('r/R (-)')
# ylabel('angle (deg)')
# legend('\alpha_{TSR=6}', "\phi_{TSR=6}",'\alpha_{TSR=8}','\phi_{TSR=8}','\alpha_{TSR=10}','\phi_{TSR=10}')

print('Done')