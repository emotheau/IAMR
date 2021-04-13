import yt
import csv
import numpy as np
import h5py
import os
from collections import OrderedDict
from yt.funcs import mylog
import sys
mylog.setLevel(40) #set yt log level to ERROR
import pandas as pd
# function to get a dictonary of snapshots in a directory ordered by code time

def getTimestampedDictXtilde(dirpath):
    dirlist = os.listdir(dirpath)
    num_plotfiles = len(dirlist)
    timestamps = np.zeros( (num_plotfiles,))
    timestamped_dict = OrderedDict()
    for ii, plotfile in enumerate(dirlist):
        if ii % 1000 ==0:
            print(ii)
        plotfilepath = os.path.join(dirpath, plotfile)
        ds = yt.load(plotfilepath)
        time = round(float(ds.current_time), 1)
        timestamped_dict[ plotfile ] = time
    sorted_d = sorted(timestamped_dict.items(), key = lambda x:x[1])
    timestamped_dict = OrderedDict(sorted_d)
    return(timestamped_dict)


def getTimestampedDictX(dirpath):
    dirlist = os.listdir(dirpath)
    num_plotfiles = len(dirlist)
    timestamps = np.zeros( (num_plotfiles,))
    timestamped_dict = OrderedDict()
    for ii, plotfile in enumerate(dirlist):
        if ii % 1000 ==0:
            print(ii)
        plotfilepath = os.path.join(dirpath, plotfile)
        ds = yt.load(plotfilepath)
        time = round(float(ds.current_time), 1)
        timestamped_dict[ time ] = plotfile
    sorted_d = sorted(timestamped_dict.items(), key = lambda x:x[0])
    timestamped_dict = OrderedDict(sorted_d)
    return(timestamped_dict)


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




basedir = sys.argv[1] #'/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/inference/3'

ml_dir = basedir + '/plotsHR'
baseline_dir = basedir + '/baselineHR'
h5dir = sys.argv[2]#'/global/cscratch1/sd/jpathak/iamr_inference/kolmogorov/3' + '/h5' 
truth_dir = '/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/256/3/plots'
truth_csv_name = "/project/projectdirs/dasrepo/jpathak/IAMR-EM/Exec/run2d/analysis/X_dict_3.csv"

if not os.path.exists(h5dir):
    os.makedirs(h5dir)
inference_file = os.path.join(h5dir, 'inference_data.h5')

res = 256

ml_dict= getTimestampedDictXtilde(ml_dir)

baseline_dict = getTimestampedDictX(baseline_dir)

truth_dict = csv_to_dict(truth_csv_name)

fields_list = ['x_velocity', 'y_velocity']

num_ml = len(list(ml_dict.keys()))

ml_array = np.zeros((num_ml, 2, res, res))
truth_array = np.zeros((num_ml, 2, res, res))
baseline_array = np.zeros((num_ml, 2, res, res))
arr_iter = 0

for idx, ml_filename in enumerate(list(ml_dict.keys())):

    print(idx)
    ml_time = ml_dict[ml_filename]
    truth_filename = truth_dict[str(ml_time)]
    print("ml time ", ml_time)
    baseline_filename = baseline_dict[ml_time]
    ml_path = os.path.join(ml_dir, ml_filename)
    truth_path = os.path.join(truth_dir, truth_filename)
    baseline_path = os.path.join(baseline_dir, baseline_filename)

    ds_ml = yt.load(ml_path)
    ad_ml = ds_ml.all_data()
    ds_truth = yt.load(truth_path)
    ad_truth = ds_truth.all_data()
    ds_baseline = yt.load(baseline_path)
    ad_baseline = ds_baseline.all_data()
    
    ml_array[idx,:, :,:] = np.stack( (np.reshape( ad_ml['x_velocity'] , (res, res)) , np.reshape(  ad_ml['y_velocity'], (res, res) ) ) , axis = 0  )
    truth_array[idx,:, :,:] = np.stack( (np.reshape( ad_truth['x_velocity'] , (res, res)) , np.reshape(  ad_truth['y_velocity'], (res, res) ) ) , axis = 0  )
    baseline_array[idx,:, :,:] = np.stack( (np.reshape( ad_baseline['x_velocity'] , (res, res)) , np.reshape(  ad_baseline['y_velocity'], (res, res) ) ) , axis = 0  )

with h5py.File(inference_file, 'w') as f:
    f.create_dataset('ml', data = ml_array[0:num_ml,:,:,:])
    f.create_dataset('truth', data = truth_array[0:num_ml,:,:,:])
    f.create_dataset('baseline', data = baseline_array[0:num_ml,:,:,:])
    f.flush()







