import yt
import numpy as np
import h5py
import os
from collections import OrderedDict
from yt.funcs import mylog
mylog.setLevel(40) #set yt log level to ERROR

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


basedir ='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/2048' #'/project/projectdirs/dasrepo/jpathak/iamr_data/ldc/case1a'

X_dir = basedir + '/plots'
X_tilde_dir = basedir + '/refined_plots'
h5dir = '/global/cscratch1/sd/jpathak/iamr_training_pairs/kolmogorov/2048/'

if not os.path.exists(h5dir):
    os.makedirs(h5dir)
train_file = os.path.join(h5dir, 'training_pairs.h5')
res = 2048

X_dict= getTimestampedDict(X_dir)
X_tilde_dict = getTimestampedDict(X_tilde_dir)
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
#        if ii % 5000 == 0:
#            if initialize == True:
#                with h5py.File(train_file, 'w') as f:
#                    f.create_dataset('fields', data = X_array, maxshape=(None, X_array.shape[1], X_array.shape[2], X_array.shape[3]))
#                    f.create_dataset('fields_tilde', data = X_tilde_array, maxshape=(None, X_array.shape[1], X_array.shape[2], X_array.shape[3]))
#                    f.flush()
#                    initialize = False
#            else:
#                with h5py.File(train_file,'a') as f:
#                    f['fields'].resize((f['fields'].shape[0] + X_array.shape[0]), axis = 0)
#                    f['fields'][-X_array.shape[0]:, :, :, :] = X_array
#                    f['fields_tilde'].resize((f['fields_tilde'].shape[0] + X_tilde_array.shape[0]), axis = 0)
#                    f['fields_tilde'][-X_tilde_array.shape[0]:, :, :, :] = X_tilde_array
#                    f.flush()    
                   
with h5py.File(train_file, 'w') as f:
    f.create_dataset('fields', data = X_array)
    f.create_dataset('fields_tilde', data = X_tilde_array)
    f.flush()
        
    
#if initialize == True:
#  with h5py.File('/global/cscratch1/sd/jpathak/rbc2d/clean_pairsRa9ds4cfl1/training_pairs.h5' , 'w') as g:
#    g.create_dataset('fields', data = fields,  maxshape=(None, fields.shape[1], fields.shape[2], fields.shape[3]))
#    g.create_dataset('fields_tilde', data = fields_tilde,  maxshape=(None, fields.shape[1], fields.shape[2], fields.shape[3]))
#    g.flush()
#    initialize = False
#else:
#  with h5py.File('/global/cscratch1/sd/jpathak/rbc2d/clean_pairsRa9ds4cfl1/training_pairs.h5' , 'a') as g:
#    g['fields'].resize((g['fields'].shape[0] + fields.shape[0]), axis = 0)
#    g['fields'][-fields.shape[0]:, :, :, :] = fields
#    g['fields_tilde'].resize((g['fields_tilde'].shape[0] + fields_tilde.shape[0]), axis = 0)
#    g['fields_tilde'][-fields_tilde.shape[0]:, :, :, :] = fields_tilde
#    g.flush()    
    

#    x_velocity_array = np.reshape(ad['x_velocity'], (res,res))
#    y_velocity_array = np.reshape(ad['y_velocity'], (res,res))
#    fields = np.stack( (x_velocity_array, y_velocity_array), axis = 0)
#    if initialize_array == True:
#        with h5py.File(output_filename, 'w') as f:
#            f.create_dataset('fields', data = fields, maxshape =(None, fields.shape[1], fields.shape[2]))
#            f.flush()
#    else:
#        with h5py.File(output_filenmae, 'a') as f:
#            f['fields'].resize((f['fields'].shape[0] + fields.shape[0]), axis = 0)
#            f['fields'][-fields.shape[0]:, :, :, :] = fields

#ds = yt.load('/project/projectdirs/dasrepo/jpathak/iamr_scripts/plotfiles/plt400000')
#
#
#field = 'x_velocity'
#time = ds.current_time
#ad = ds.all_data()
#field_array = np.reshape(ad[field], (256,256))
#with h5py.File('testout.h5', 'w') as f:
#    f['x_velocity'] = field_array
#    f['time'] = time

