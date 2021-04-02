import h5py
import numpy as np
import os

path = '/global/cscratch1/sd/jpathak/iamr_training_pairs/kolmogorov/256/flt32' 
outputpath = '/global/cscratch1/sd/jpathak/iamr_training_pairs/kolmogorov/256/flt32' 
if not os.path.isdir(outputpath):
    os.makedirs(outputpath)


for filename in os.listdir(path):
    if filename.endswith("_5.h5"):
        print(os.path.join(path, filename))
        filepath = os.path.join(path, filename)
        with h5py.File(filepath, 'r') as f:
            fields = f['fields'][:,:,:,:].astype(np.float32)
            fields_tilde = f['fields_tilde'][:,:,:,:].astype(np.float32)
        output_filepath = os.path.join(outputpath, filename)
        with h5py.File(output_filepath, 'w') as g:
            g.create_dataset('fields', data = fields)
            g.create_dataset('fields_tilde', data = fields_tilde)
    else:
        continue



