import yt
import csv
import numpy as np
import h5py
import os
from collections import OrderedDict
from yt.funcs import mylog


def csv_to_dict(filename):
    file = open(filename, mode='r')
    odict = OrderedDict()
    csvReader = csv.reader(file)
    keys = next(csvReader)
    values = next(csvReader)
    for key, value in zip(keys, values):
        odict[key] = value
    return(odict)



truth_csv_name = "X_dict_3.csv"


truth_dict = csv_to_dict(truth_csv_name)
i=0
for key, val in truth_dict.items():
    print(key, val)
    i+=1
    if i > 10:
        break

print(truth_dict[str(3081.9)])
