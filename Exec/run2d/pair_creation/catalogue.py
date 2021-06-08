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
   
idx_list = [0, 1, 2,3,4,5,6,7,8, 9, 10, 11, 12, 13, 14, 15]

#idx = str(sys.argv[1])
res = 2048
ur = 16
for idx in idx_list:
    
    basedir ='/project/projectdirs/dasrepo/jpathak/iamr_expts/kolmogorov/processed_files/'+ str(res) + '/ur' + str(ur) + '/' + str(idx) 
    
    X_dir = basedir + '/plots'
    X_tilde_dir = basedir + '/refined_plots'
#    X_tilde_lr_dir = basedir + '/tilde_plots'
    catalogdir = basedir + '/csv/'
    
    if not os.path.exists(catalogdir):
        os.makedirs(catalogdir)
    
    X_dict= getTimestampedDictX(X_dir)
    X_tilde_dict = getTimestampedDictXtilde(X_tilde_dir)
#    X_tilde_lr_dict = getTimestampedDictX(X_tilde_lr_dir)
    
    X_dict_name = catalogdir + "X_dict_" + str(res) + "ur" + str(ur) + "_" + str(idx)
    X_tilde_dict_name = catalogdir + "X_tilde_dict_" + str(res) + "ur" + str(ur) + "_" + str(idx)
#    X_tilde_lr_dict_name = catalogdir + "X_tilde_dict_lr_" + str(res) + "ur" + str(ur) + "_" + str(idx)
    
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
    
    
    dict_to_csv(X_dict, X_dict_name)
    dict_to_csv(X_tilde_dict, X_tilde_dict_name)
#    dict_to_csv(X_tilde_lr_dict, X_tilde_lr_dict_name)
    
    
    
