import yt
import csv
import numpy as np
import h5py
import os
from collections import OrderedDict
from yt.funcs import mylog
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
   


idx = '5'

basedir ='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/256/' + idx #'/project/projectdirs/dasrepo/jpathak/iamr_data/ldc/case1a'

X_dir = basedir + '/plots'
X_tilde_dir = basedir + '/refined_plots'
h5dir = '/global/cscratch1/sd/jpathak/iamr_training_pairs/kolmogorov/256/'

#if not os.path.exists(h5dir):
#    os.makedirs(h5dir)
#train_file = os.path.join(h5dir, 'training_pairs.h5')
res = 256

X_dict= getTimestampedDictX(X_dir)
X_tilde_dict = getTimestampedDictXtilde(X_tilde_dir)
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


dict_to_csv(X_dict, "X_dict_" + idx )
dict_to_csv(X_tilde_dict, "X_tilde_dict_" + idx)

#X_dict = csv_to_dict("X_dict.csv")
#X_tilde_dict = csv_to_dict("X_tilde_dict.csv")


