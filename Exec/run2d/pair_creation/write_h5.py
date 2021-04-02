import yt
import numpy as np
import h5py
import os
from collections import OrderedDict
from yt.funcs import mylog
mylog.setLevel(40) #set yt log level to ERROR
import json
# function to get a dictonary of snapshots in a directory ordered by code time
def getTimestampedDict(dirpath):

    dirlist = os.listdir(dirpath)
    num_plotfiles = len(dirlist)
    timestamps = np.zeros( (num_plotfiles,))
    timestamped_dict = OrderedDict()
    for ii, plotfile in enumerate(dirlist):
        if ii % 1000 ==0:
            print(ii)
        plotfilepath = os.path.join(dirpath, plotfile)
        ds = yt.load(plotfilepath)
        time = ds.current_time
        timestamped_dict[ plotfile ] = float(time)
    sorted_d = sorted(timestamped_dict.items(), key = lambda x:x[1])
    timestamped_dict = OrderedDict(sorted_d)
    return(timestamped_dict)


basedir ='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/256' #'/project/projectdirs/dasrepo/jpathak/iamr_data/ldc/case1a'

X_dir = basedir + '/plots'
X_tilde_dir = basedir + '/refined_plots'
h5dir = '/global/cscratch1/sd/jpathak/iamr_training_pairs/kolmogorov/256/'

if not os.path.exists(h5dir):
    os.makedirs(h5dir)
train_file = os.path.join(h5dir, 'training_pairs.h5')
res = 256

X_dict= getTimestampedDict(X_dir)
X_tilde_dict = getTimestampedDict(X_tilde_dir)


with open('odict_fields.json', 'w') as f:
    f.write(json.dumps(X_dict))

with open('odict_fields_tilde.json', 'w') as f:
    f.write(json.dumps(X_tilde_dir))

fields_list = ['x_velocity', 'y_velocity']

num_X = len(list(X_dict.keys()))

X_array = np.zeros((num_X-1, 2, res, res))
X_tilde_array = np.zeros((num_X-1, 2, res, res))
initialize = True
for ii, (filename, timestamp) in enumerate(X_tilde_dict.items()):
    if ii < num_X -1:
        X_tilde_path = os.path.join(X_tilde_dir, filename)
        ds_X_tilde = yt.load(X_tilde_path)
        ad_X_tilde = ds_X_tilde.all_data()
        X_tilde_array[ii,:, :,:] = np.stack( (np.reshape( ad_X_tilde['x_velocity'] , (res, res)) , np.reshape(  ad_X_tilde['y_velocity'], (res, res) ) ) , axis = 0  )
        time_tilde = ds_X_tilde.current_time
        X_filename = list(X_dict.keys())[ii+1]
        X_path = os.path.join(X_dir, X_filename)
        ds_X = yt.load(X_path)
        ad_X = ds_X.all_data()
        X_array[ii,:,:,:] = np.stack( (np.reshape( ad_X['x_velocity'] , (res, res)) , np.reshape(  ad_X['y_velocity'], (res, res) ) ) , axis = 0  )
        time = ds_X.current_time
        print(ii)
        print(float(time_tilde) - float(time))
        print(time_tilde, 'time_tilde')
        print(time, 'time')
        if ii % 5000 == 0:
            if initialize == True:
                with h5py.File(train_file, 'w') as f:
                    f.create_dataset('fields', data = X_array, maxshape=(None, X_array.shape[1], X_array.shape[2], X_array.shape[3]))
                    f.create_dataset('fields_tilde', data = X_tilde_array, maxshape=(None, X_array.shape[1], X_array.shape[2], X_array.shape[3]))
                    f.flush()
                    initialize = False
            else:
                with h5py.File(train_file,'a') as f:
                    f['fields'].resize((f['fields'].shape[0] + X_array.shape[0]), axis = 0)
                    f['fields'][-X_array.shape[0]:, :, :, :] = X_array
                    f['fields_tilde'].resize((f['fields_tilde'].shape[0] + X_tilde_array.shape[0]), axis = 0)
                    f['fields_tilde'][-X_tilde_array.shape[0]:, :, :, :] = X_tilde_array
                    f.flush()    
                   
with h5py.File(train_file, 'w') as f:
    f.create_dataset('fields', data = X_array)
    f.create_dataset('fields_tilde', data = X_tilde_array)
    f.flush()
        
    

