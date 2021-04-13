import yt
import sys
import numpy as np
import h5py
import os
from collections import OrderedDict
import csv
from yt.funcs import mylog
mylog.setLevel(40) #set yt log level to ERROR
import json

def dict_to_csv(odict_object, filename):
    
    keys, values = [], []
    
    for key, value in odict_object.items():
        keys.append(key)
        values.append(value)
    
    with open(filename + ".csv", "w") as outfile:
        csvwriter = csv.writer(outfile)
        csvwriter.writerow(keys)
        csvwriter.writerow(values)

def csv_to_dict(filename):
    file = open(filename, mode='r')
    odict = OrderedDict()
    csvReader = csv.reader(file)
    keys = next(csvReader)
    values = next(csvReader)
    for key, value in zip(keys, values):
        odict[key] = value
    return(odict)

chunk_size = 1000

ii = int(sys.argv[1])

startidx = chunk_size*ii
endidx = chunk_size*(ii+1)

idx = int(sys.argv[2])

X_dict = csv_to_dict("X_dict_ur8"+ str(idx) + ".csv")
X_tilde_dict = csv_to_dict("X_tilde_dict_ur8" +str(idx) +  ".csv")



basedir ='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/256/ur8/' +  str(idx) #'/project/projectdirs/dasrepo/jpathak/iamr_data/ldc/case1a'

X_dir = basedir + '/plots'
X_tilde_dir = basedir + '/refined_plots'
h5dir = '/global/cscratch1/sd/jpathak/iamr_training_pairs/kolmogorov/256/ur8/h5_with_history/'

if ii == 0:
    if not os.path.exists(h5dir):
        os.makedirs(h5dir)

train_file = os.path.join(h5dir, 'training_pairs_h1_' + str(ii)+ '_' + str(idx) +  '.h5')

res = 256

fields_list = ['x_velocity', 'y_velocity']

num_X = len(list(X_dict.keys()))

X_array = np.zeros((chunk_size, 2, res, res))
X_tilde_array = np.zeros((chunk_size, 4, res, res))
arr_iter = 0
for dictidx in range(startidx, endidx):
    tilde_filename = list(X_tilde_dict.keys())[dictidx]
    tilde_time = X_tilde_dict[tilde_filename]
    t_minus_time = str(round(float(tilde_time) - 0.1,1))
    X_tilde_path = os.path.join(X_tilde_dir, tilde_filename)
    ds_X_tilde = yt.load(X_tilde_path)
    ad_X_tilde = ds_X_tilde.all_data()
    if (tilde_time in X_dict) and (t_minus_time in X_dict):
        X_tilde_array[arr_iter,0:2, :,:] = np.stack( (np.reshape( ad_X_tilde['x_velocity'] , (res, res)) , np.reshape(  ad_X_tilde['y_velocity'], (res, res) ) ) , axis = 0  ).astype(np.float32)
        X_filename = X_dict[tilde_time]    
        tilde_time_from_file = round(float(ds_X_tilde.current_time), 1)
#        print(tilde_time_from_file, tilde_time)
        X_path = os.path.join(X_dir, X_filename)
        ds_X = yt.load(X_path)
        ad_X = ds_X.all_data()
        X_array[arr_iter,:,:,:] = np.stack( (np.reshape( ad_X['x_velocity'] , (res, res)) , np.reshape(  ad_X['y_velocity'], (res, res) ) ) , axis = 0  ).astype(np.float32)
#        print( ds_X.current_time)
        X_minus_filename = X_dict[t_minus_time]
        X_minus_path = os.path.join(X_dir, X_minus_filename)
        ds_X_minus = yt.load(X_minus_path)
        ad_X_minus = ds_X_minus.all_data()
        X_tilde_array[arr_iter,2:4,:,:] = np.stack( (np.reshape( ad_X_minus['x_velocity'] , (res, res)) , np.reshape(  ad_X_minus['y_velocity'], (res, res) ) ) , axis = 0  ).astype(np.float32)

        arr_iter += 1
print("valid_pairs=", arr_iter)
with h5py.File(train_file, 'w') as f:
    f.create_dataset('fields', data = X_array[0:arr_iter,:,:,:].astype(np.float32))
    f.create_dataset('fields_tilde', data = X_tilde_array[0:arr_iter,:,:,:].astype(np.float32))
    f.flush()


