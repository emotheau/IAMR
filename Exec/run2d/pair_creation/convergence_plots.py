import yt
import h5py
import numpy as np


filename = '/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/convergence_test/512/plots69696/'
res = 512
ds = yt.load(filename)
ad = ds.all_data()
X = np.reshape( ad['x_velocity'], (res, res))

with h5py.File('/global/cscratch1/sd/jpathak/convergence_test/512/plots69696.h5', 'w') as f:
    f.create_dataset('fields' ,data = X)
    f.flush()

