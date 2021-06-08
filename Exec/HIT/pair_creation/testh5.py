import os
import yt
import h5py
import numpy as np

filepath='/global/cscratch1/sd/jpathak/iamr_hit/128/'

truth = 'plots/plots00482'
tilde_plot = 'refined_plot/plots_00480'

#plotfiles = ['plots00475', 'plots00393' , 'plots00000'] 

outdir = '/global/cscratch1/sd/jpathak/iamr_hit/128/h5/'

if not os.path.isdir(outdir):
    os.makedirs(outdir)


ds = yt.load(os.path.join(filepath, truth))
ad = ds.all_data()
res = 128
X = np.reshape(ad["x_velocity"], (res, res, res))
print(X.shape)

ds = yt.load(os.path.join(filepath, tilde_plot))
ad = ds.all_data()
res = 128
X_tilde = np.reshape(ad["x_velocity"], (res, res, res))
print(X.shape)

plotfile = 'residual_test.h5'

with h5py.File(os.path.join(outdir, plotfile), 'w') as f:
    f.create_dataset("x_velocity_residual", data  = X-X_tilde)
    f.flush()
