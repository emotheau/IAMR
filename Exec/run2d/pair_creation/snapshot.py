import numpy as np
import matplotlib.pyplot as plt
import h5py


with h5py.File('/global/cscratch1/sd/jpathak/iamr_training_pairs/kolmogorov/2048/training_pairs.h5', 'r') as f:
  print(list(f.keys()))
  fields_hr = f['fields'][0:3,:,:,:]
  fields_tilde_upsampled = f['fields_tilde'][0:3,:,:,:]

energy_hr = fields_hr[:,0:1,:,:]**2 + fields_hr[:,1:2,:,:]**2
energy_upsampled = fields_tilde_upsampled[:,0:1,:,:]**2 + fields_tilde_upsampled[:,1:2,:,:]**2
energy_residual = energy_hr[0,0,:,:] - energy_upsampled[0,0,:,:]
arr = np.concatenate((energy_residual[0:-1,0:-1],energy_residual[0:-1,0:-1],energy_residual[0:-1,0:-1]), axis = 1) 

plt.rcParams["figure.figsize"] = (150,50)
plt.figure()
im = plt.imshow(arr, cmap = 'RdBu', vmax = 0.002, vmin =-0.002)

plt.savefig('snapshot.eps')
