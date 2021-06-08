import h5py
import numpy as np
import os
import glob

mypath = '/global/cscratch1/sd/jpathak/iamr_training_pairs/kolmogorov/2048/ur16/h5_history/'
filearray = glob.glob(mypath + '*.h5')

im_in_file = []
chunk_size = 200
num_res = 0
total_im = 0
for filename in filearray:
    if 'train' in filename:
        with h5py.File(filename, 'r') as f:
            num_im = f['fields'].shape[0]
            print(num_im)
            im_in_file.append(num_im)
            total_im+=num_im

print(total_im)
print(im_in_file)
num_chunks = total_im // chunk_size

        



    
