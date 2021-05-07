
import numpy as np
import h5py

f = h5py.File('BEM_data.mat','r')
print(f.keys())

BEM_rR = np.asarray(f.get('rR'))
BEM_alpha = np.asarray(f.get('alpha'))


print('Done')