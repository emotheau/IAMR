import os
import yt
import h5py
import numpy as np

filepath='/global/cscratch1/sd/jpathak/iamr_hit/128/plots/'

plotfiles = ['plots00475', 'plots00393' , 'plots00000'] 

outdir = '/global/cscratch1/sd/jpathak/iamr_hit/128/h5/'

if not os.path.isdir(outdir):
    os.makedirs(outdir)

for plotfile in plotfiles:
    print(plotfile)

    ds = yt.load(os.path.join(filepath, plotfile))
    ad = ds.all_data()
    res = 128
    X = np.reshape(ad["x_velocity"], (res, res, res))
    print(X.shape)
    
    with h5py.File(os.path.join(outdir, plotfile), 'w') as f:
        f.create_dataset("x_velocity", data  = X)
        f.flush()
