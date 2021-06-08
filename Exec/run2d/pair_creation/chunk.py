import h5py
import numpy as np
import os

filepath = '/global/cscratch1/sd/jpathak/iamr_training_pairs/kolmogorov/2048/ur16/training_pairs_full.h5'

chunk_path = '/global/cscratch1/sd/jpathak/iamr_training_pairs/kolmogorov/2048/ur16/chunks'

if not os.path.isdir(chunk_path):
    os.makedirs(chunk_path)
 
chunk_size = 200

with h5py.File(filepath, 'r') as f:

    num_images = f['fields'].shape[0]
    print(num_images)

    num_chunks = num_images // chunk_size

    print(num_chunks)

    for ii in range(num_chunks):

        outfilepath = os.path.join(chunk_path, 'training_pairs_' + str(ii) + '.h5')
        print(outfilepath)

        with h5py.File(os.path.join(chunk_path, 'training_pairs_' + str(ii) + '.h5'), 'w' ) as g:

            g.create_dataset('fields', data = f['fields'][chunk_size*ii:chunk_size*(ii+1),:,:,:])
            g.create_dataset('fields_tilde', data = f['fields_tilde'][chunk_size*ii:chunk_size*(ii+1),:,:,:])
            g.flush()


