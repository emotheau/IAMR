import numpy as np
import h5py
import os.path
import glob

mypath = '/global/cscratch1/sd/jpathak/iamr_training_pairs/kolmogorov/2048/ur16/h5_history/'
filearray = glob.glob(mypath + '*.h5')

outfilepath = '/global/cscratch1/sd/jpathak/iamr_training_pairs/kolmogorov/2048/ur16/training_pairs_full.h5'

initialize = True

for filepath in filearray:
  print(filepath, flush = True)
  if 'train' in filepath:
    with h5py.File(filepath, 'r') as f:
      print(f['fields'].shape, flush = True)
      print(f['fields_tilde'].shape, flush = True)
      fields = f['fields'][:,:,:,:]
      fields_tilde = f['fields_tilde'][:,0:2,:,:]
      if initialize == True:
        with h5py.File(outfilepath , 'w') as g:
          g.create_dataset('fields', data = fields,  maxshape=(None, fields.shape[1], fields.shape[2], fields.shape[3]))
          g.create_dataset('fields_tilde', data = fields_tilde,  maxshape=(None, fields.shape[1], fields.shape[2], fields.shape[3]))
          g.flush()
          initialize = False
      else:
        with h5py.File(outfilepath , 'a') as g:
          g['fields'].resize((g['fields'].shape[0] + fields.shape[0]), axis = 0)
          g['fields'][-fields.shape[0]:, :, :, :] = fields
          g['fields_tilde'].resize((g['fields_tilde'].shape[0] + fields_tilde.shape[0]), axis = 0)
          g['fields_tilde'][-fields_tilde.shape[0]:, :, :, :] = fields_tilde
          g.flush()




